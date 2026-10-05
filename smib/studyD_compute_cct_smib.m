function out = studyD_compute_cct_smib(cluster_pre, faultModel, opts, mode)
% ias_compute_cct_smib for a chosen limiter mode. A clearing angle beyond
% the post-fault UEP counts as unstable whatever the energy: under
% q-priority the potential drops again past the UEP, and without this the
% verdict would not be monotone in the clearing time.

if ~isfield(faultModel,'kV'), faultModel.kV = 0.2; end
if ~isfield(faultModel,'kX'), faultModel.kX = 1.0; end
if nargin < 3 || isempty(opts), opts = struct; end
if ~isfield(opts,'t_hi'), opts.t_hi = 2.0; end
if ~isfield(opts,'tol'),  opts.tol  = 1e-3; end
if ~isfield(opts,'f0'),   opts.f0   = 60;   end
if nargin < 4 || isempty(mode), mode = 'circular'; end

if isfield(opts,'delta_h')
    delta_h = opts.delta_h(:);
else
    delta_h = linspace(-pi, +pi, 4000).';
end

cluster_post = cluster_pre;

Pe_pre  = ias_compute_pe_lim(cluster_pre,  delta_h, mode);
[delta_SEP_pre, ~, ~] = ias_find_equilibria(delta_h, Pe_pre, cluster_pre.Pref_h);

Pe_post = ias_compute_pe_lim(cluster_post, delta_h, mode);
[delta_SEP_post, delta_UEP_post, ~] = ias_find_equilibria(delta_h, Pe_post, cluster_post.Pref_h);

if isnan(delta_SEP_pre) || isnan(delta_SEP_post) || isnan(delta_UEP_post)
    out.CCT  = 0.0;
    out.note = 'no SEP/UEP pair, CCT set to 0';
    out.delta_SEP_pre  = delta_SEP_pre;
    out.delta_SEP_post = delta_SEP_post;
    out.delta_UEP_post = delta_UEP_post;
    out.Vcrit          = NaN;
    out.limiter_mode   = mode;
    return;
end

[U_post, Vcrit] = ias_energy_barrier(delta_h, Pe_post, cluster_post.Pref_h, ...
                                     delta_SEP_post, delta_UEP_post);

cluster_fault         = cluster_pre;
cluster_fault.Vth_mag = cluster_pre.Vth_mag * faultModel.kV;
cluster_fault.Xnet    = cluster_pre.Xnet    * faultModel.kX;

Pe_fault = ias_compute_pe_lim(cluster_fault, delta_h, mode);

PeF = griddedInterpolant(delta_h, Pe_fault, 'linear', 'nearest');
Up  = griddedInterpolant(delta_h, U_post,  'linear', 'nearest');

f0 = opts.f0; ws = 2*pi*f0;
H = cluster_pre.Hh; D = cluster_pre.Dh; Pm = cluster_pre.Pref_h;

odefun = @(t,x) [ x(2); (ws/(2*H))*(Pm - PeF(x(1))) - (D/(2*H))*x(2) ];

M = 2*H/ws;
Eclear = @(xcl) 0.5*M*(xcl(2)^2) + Up(xcl(1));

isStable = @(tc) local_stable(tc, odefun, delta_SEP_pre, Eclear, Vcrit, delta_UEP_post);

if isStable(opts.t_hi)
    out.CCT = opts.t_hi;
else
    lo = 0.0; hi = opts.t_hi;
    while (hi - lo) > opts.tol
        mid = 0.5*(lo + hi);
        if isStable(mid), lo = mid; else, hi = mid; end
    end
    out.CCT = lo;
end

out.delta_SEP_pre  = delta_SEP_pre;
out.delta_SEP_post = delta_SEP_post;
out.delta_UEP_post = delta_UEP_post;
out.Vcrit          = Vcrit;
out.limiter_mode   = mode;
out.Pe_pre         = Pe_pre;
out.Pe_post        = Pe_post;
out.Pe_fault       = Pe_fault;
end

function ok = local_stable(tc, odefun, delta0, Eclear, Vcrit, delta_UEP)
    x0 = [delta0; 0.0];
    sol = ode45(odefun, [0 tc], x0);
    xcl = sol.y(:,end);
    ok  = (xcl(1) < delta_UEP) && (Eclear(xcl) < Vcrit);
end
