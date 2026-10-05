function out = studyE_compute_cct_smib(cluster_pre, faultModel, opts, mode)
% studyD_compute_cct_smib with a series resistance and a choice of frame
% (both through ias_compute_pe_lim_R). The fault scales Rnet together with
% Xnet, and the equilibria come from ias_find_equilibria_fwd. Extra fields:
%   feasible     a SEP with a UEP above it exists before and after the fault
%   censored     still stable at opts.t_hi, so CCT = t_hi is a lower bound
%   Aacc         accelerating area on the fault-on curve, SEP to UEP
%   fracNegIP0   share of the SEP-UEP interval on which the unlimited
%                in-phase current is negative

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

Pe_pre = ias_compute_pe_lim_R(cluster_pre, delta_h, mode);
[delta_SEP_pre, ~, ~, ok_pre] = ias_find_equilibria_fwd(delta_h, Pe_pre, cluster_pre.Pref_h);

Pe_post = ias_compute_pe_lim_R(cluster_post, delta_h, mode);
[delta_SEP_post, delta_UEP_post, ~, ok_post] = ...
    ias_find_equilibria_fwd(delta_h, Pe_post, cluster_post.Pref_h);

out.limiter_mode   = mode;
out.delta_SEP_pre  = delta_SEP_pre;
out.delta_SEP_post = delta_SEP_post;
out.delta_UEP_post = delta_UEP_post;

if ~ok_pre || ~ok_post
    out.CCT        = NaN;
    out.feasible   = false;
    out.censored   = false;
    out.Vcrit      = NaN;
    out.Aacc       = NaN;
    out.fracNegIP0 = NaN;
    out.note = 'no SEP with a UEP above it';
    return;
end
out.feasible = true;

[U_post, Vcrit] = ias_energy_barrier(delta_h, Pe_post, cluster_post.Pref_h, ...
                                     delta_SEP_post, delta_UEP_post);

cluster_fault         = cluster_pre;
cluster_fault.Vth_mag = cluster_pre.Vth_mag * faultModel.kV;
cluster_fault.Xnet    = cluster_pre.Xnet    * faultModel.kX;
if isfield(cluster_pre,'Rnet') && ~isempty(cluster_pre.Rnet)
    cluster_fault.Rnet = cluster_pre.Rnet * faultModel.kX;
end

Pe_fault = ias_compute_pe_lim_R(cluster_fault, delta_h, mode);

PeF = griddedInterpolant(delta_h, Pe_fault, 'linear', 'nearest');
Up  = griddedInterpolant(delta_h, U_post,  'linear', 'nearest');

f0 = opts.f0; ws = 2*pi*f0;
H = cluster_pre.Hh; D = cluster_pre.Dh; Pm = cluster_pre.Pref_h;

odefun = @(t,x) [ x(2); (ws/(2*H))*(Pm - PeF(x(1))) - (D/(2*H))*x(2) ];

M = 2*H/ws;
Eclear = @(xcl) 0.5*M*(xcl(2)^2) + Up(xcl(1));

isStable = @(tc) local_stable(tc, odefun, delta_SEP_pre, Eclear, Vcrit, delta_UEP_post);

if isStable(opts.t_hi)
    out.CCT      = opts.t_hi;
    out.censored = true;
else
    out.censored = false;
    lo = 0.0; hi = opts.t_hi;
    while (hi - lo) > opts.tol
        mid = 0.5*(lo + hi);
        if isStable(mid), lo = mid; else, hi = mid; end
    end
    out.CCT = lo;
end

selA = delta_h >= delta_SEP_pre & delta_h <= delta_UEP_post;
if nnz(selA) > 1
    out.Aacc = trapz(delta_h(selA), Pm - Pe_fault(selA));
else
    out.Aacc = NaN;
end

Zn = cluster_post.Xnet;
if isfield(cluster_post,'Rnet') && ~isempty(cluster_post.Rnet)
    Zc = cluster_post.Rnet + 1j*Zn;
else
    Zc = 1j*Zn;
end
I0 = (cluster_post.Eh_mag .* exp(1j*delta_h) - ...
      cluster_post.Vth_mag .* exp(1j*cluster_post.theta_th)) ./ Zc;
IP0 = real(I0 .* exp(-1j*cluster_post.theta_th));
selB = delta_h >= delta_SEP_post & delta_h <= delta_UEP_post;
if nnz(selB) > 1
    out.fracNegIP0 = nnz(IP0(selB) < 0) / nnz(selB);
else
    out.fracNegIP0 = NaN;
end

out.Vcrit    = Vcrit;
out.Pe_pre   = Pe_pre;
out.Pe_post  = Pe_post;
out.Pe_fault = Pe_fault;
end

function ok = local_stable(tc, odefun, delta0, Eclear, Vcrit, delta_UEP)
    x0 = [delta0; 0.0];
    sol = ode45(odefun, [0 tc], x0);
    xcl = sol.y(:,end);
    ok  = (xcl(1) < delta_UEP) && (Eclear(xcl) < Vcrit);
end
