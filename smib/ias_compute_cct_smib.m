function out = ias_compute_cct_smib(cluster_pre, faultModel, opts)
% Critical clearing time of a single cluster, eqs. (33)-(37). The swing
% equation is integrated along the fault-on curve up to the clearing time,
% the post-fault energy of the clearing state is compared with the barrier,
% and the clearing time is found by bisection.
%   faultModel.kV   fault-on Thevenin voltage relative to pre-fault (0.2)
%   faultModel.kX   scaling of Xnet during the fault (1)
%   opts.t_hi       upper end of the search window, s (2)
%   opts.tol        bisection tolerance, s (1e-3)
%   opts.delta_h    angle grid, rad
% The post-fault network is the pre-fault one. A cluster that is still
% stable at t_hi gets CCT = t_hi.

if ~isfield(faultModel,'kV'), faultModel.kV = 0.2; end
if ~isfield(faultModel,'kX'), faultModel.kX = 1.0; end

if nargin < 3, opts = struct; end
if ~isfield(opts,'t_hi'),  opts.t_hi = 2.0; end
if ~isfield(opts,'tol'),   opts.tol  = 1e-3; end
if ~isfield(opts,'f0'),    opts.f0   = 60; end

if isfield(opts,'delta_h')
    delta_h = opts.delta_h(:);
else
    delta_h = linspace(-pi, +pi, 4000).';
end

cluster_post = cluster_pre;

Pe_pre = ias_compute_pe_lim(cluster_pre, delta_h);
[delta_SEP_pre, ~, ~] = ias_find_equilibria(delta_h, Pe_pre, cluster_pre.Pref_h);

Pe_post = ias_compute_pe_lim(cluster_post, delta_h);
[delta_SEP_post, delta_UEP_post, ~] = ias_find_equilibria(delta_h, Pe_post, cluster_post.Pref_h);

if isnan(delta_SEP_pre) || isnan(delta_SEP_post) || isnan(delta_UEP_post)
    out.CCT  = 0.0;
    out.note = 'no SEP/UEP pair, CCT set to 0';
    return;
end

[U_post, Vcrit] = ias_energy_barrier(delta_h, Pe_post, cluster_post.Pref_h, delta_SEP_post, delta_UEP_post);

cluster_fault = cluster_pre;
cluster_fault.Vth_mag = cluster_pre.Vth_mag * faultModel.kV;
cluster_fault.Xnet    = cluster_pre.Xnet    * faultModel.kX;

Pe_fault = ias_compute_pe_lim(cluster_fault, delta_h);

PeF = griddedInterpolant(delta_h, Pe_fault, 'linear', 'nearest');
Up  = griddedInterpolant(delta_h, U_post,  'linear', 'nearest');

f0 = opts.f0;
ws = 2*pi*f0;

H = cluster_pre.Hh;
D = cluster_pre.Dh;
Pm = cluster_pre.Pref_h;

% x = [delta; d(delta)/dt] in rad and rad/s
odefun = @(t,x) [ ...
    x(2); ...
    (ws/(2*H))*(Pm - PeF(x(1))) - (D/(2*H))*x(2) ];

% kinetic plus potential energy of the clearing state, eq. (34), with the
% speed term written for d(delta)/dt in rad/s
M = 2*H/ws;
Eclear = @(xcl) 0.5*M*(xcl(2)^2) + Up(xcl(1));

isStable = @(tc) local_stable(tc, odefun, delta_SEP_pre, Eclear, Vcrit);

if isStable(opts.t_hi)
    out.CCT = opts.t_hi;
else
    lo = 0.0; hi = opts.t_hi;
    while (hi - lo) > opts.tol
        mid = 0.5*(lo + hi);
        if isStable(mid)
            lo = mid;
        else
            hi = mid;
        end
    end
    out.CCT = lo;
end

out.delta_SEP_pre  = delta_SEP_pre;
out.delta_SEP_post = delta_SEP_post;
out.delta_UEP_post = delta_UEP_post;
out.Vcrit          = Vcrit;

end

function ok = local_stable(tc, odefun, delta0, Eclear, Vcrit)
x0 = [delta0; 0.0];
sol = ode45(odefun, [0 tc], x0);
xcl = sol.y(:,end);
ok  = (Eclear(xcl) < Vcrit);
end
