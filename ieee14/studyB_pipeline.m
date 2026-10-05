function out = studyB_pipeline(S, opts)
% Steps 2 to 4 as one function, for a state struct S that the caller may
% have perturbed: stage Thevenin equivalents, post-fault equilibria and
% margin, energy barrier and CCT of each cluster.
%   opts.tc_max, .rtol, .atol, .tol_tc, .max_it   as in step4
%   opts.t_clear_threshold   clearing time the system CCT is compared
%                            with, s (0.1)

define_constants;
if nargin < 2, opts = struct; end
if ~isfield(opts,'tc_max'),  opts.tc_max  = 10.0; end
if ~isfield(opts,'rtol'),    opts.rtol    = 1e-7; end
if ~isfield(opts,'atol'),    opts.atol    = 1e-9; end
if ~isfield(opts,'tol_tc'),  opts.tol_tc  = 1e-6; end
if ~isfield(opts,'max_it'),  opts.max_it  = 70;   end
if ~isfield(opts,'t_clear_threshold'), opts.t_clear_threshold = 0.100; end

Nscan       = 8001;
delta_grid  = linspace(-2*pi, 2*pi, Nscan).';
fmin_tol    = 1e-4;
root_dedup  = 1e-6;
wrap_dedup  = 1e-6;
min_alpha   = 1e-7;
eq_tol      = 1e-6;
Iref        = 1e6;

wrapPi  = @(x) atan2(sin(x), cos(x));
wrap2pi = @(x) mod(x, 2*pi);
circDist = @(a,b) abs(atan2(sin(a-b), cos(a-b)));

t_start = tic;

% rebuilt here because Vop, Sload or the reactances may have been changed
Yload = conj(S.Sload) ./ max(1e-12, abs(S.Vop).^2);
Yload_mat = spdiags(Yload, 0, S.nb, S.nb);

Zs_sg  = 1j*S.Xdprime_sg;  Ys_sg  = 1/Zs_sg;
Zf_gfm = 1j*S.Xf_gfm;      Yf_gfm = 1/Zf_gfm;

Ypre   = S.Ybus_pre;
Yfault = Ypre;
kf = S.fault_bus;
Yfault(kf,kf) = Yfault(kf,kf) + 1/S.Zfault_pu;

mpc_post = S.mpc0;
br = mpc_post.branch;
for kk = 1:size(S.trip_pairs,1)
    a = S.trip_pairs(kk,1); b = S.trip_pairs(kk,2);
    idx = find(((br(:,F_BUS)==a & br(:,T_BUS)==b) | (br(:,F_BUS)==b & br(:,T_BUS)==a)) & br(:,BR_STATUS)==1, 1);
    % line already out before the fault: tripping it changes nothing
    if isempty(idx)
        idx = find((br(:,F_BUS)==a & br(:,T_BUS)==b) | (br(:,F_BUS)==b & br(:,T_BUS)==a), 1);
    end
    if isempty(idx)
        error('studyB_pipeline: post-fault trip line %d-%d not found.', a, b);
    end
    br(idx, BR_STATUS) = 0;
end
mpc_post.branch = br;
[Ypost, ~, ~] = makeYbus(S.baseMVA, mpc_post.bus, mpc_post.branch);

gfm_buses = S.gfm_buses;
Vth_pre   = containers.Map('KeyType','double','ValueType','any');
Zth_pre   = containers.Map('KeyType','double','ValueType','any');
Vth_fault = containers.Map('KeyType','double','ValueType','any');
Zth_fault = containers.Map('KeyType','double','ValueType','any');
Vth_post  = containers.Map('KeyType','double','ValueType','any');
Zth_post  = containers.Map('KeyType','double','ValueType','any');

for tbus = gfm_buses
    [vp, zp] = stage_thevenin(Ypre,   Yload_mat, S.nb, S.sg_bus, Zs_sg, Ys_sg, S.gfl_buses, gfm_buses, S.Sgen_bus, S.Vop, Zf_gfm, Yf_gfm, tbus);
    [vf, zf] = stage_thevenin(Yfault, Yload_mat, S.nb, S.sg_bus, Zs_sg, Ys_sg, S.gfl_buses, gfm_buses, S.Sgen_bus, S.Vop, Zf_gfm, Yf_gfm, tbus);
    [vo, zo] = stage_thevenin(Ypost,  Yload_mat, S.nb, S.sg_bus, Zs_sg, Ys_sg, S.gfl_buses, gfm_buses, S.Sgen_bus, S.Vop, Zf_gfm, Yf_gfm, tbus);
    Vth_pre(tbus)   = vp; Zth_pre(tbus)   = zp;
    Vth_fault(tbus) = vf; Zth_fault(tbus) = zf;
    Vth_post(tbus)  = vo; Zth_post(tbus)  = zo;
end

delta_g = containers.Map('KeyType','double','ValueType','double');
for b = gfm_buses
    Iop = conj(S.Sgen_bus(b) ./ S.Vop(b));
    Eop = S.Vop(b) + Zf_gfm * Iop;
    delta_g(b) = wrapPi(angle(Eop) - angle(Vth_pre(b)));
end

nb_g = numel(gfm_buses);
MhCL = nan(nb_g,1);
Delta = nan(nb_g,1);
DeltaRef = nan(nb_g,1);
CCT_h = nan(nb_g,1);
Ecrit_h = nan(nb_g,1);
stable_to_tcmax_h = false(nb_g,1);
CCT_is_gt_h = false(nb_g,1);

for ii = 1:nb_g
    b = gfm_buses(ii);
    Pm   = S.Pm_map(b);
    Eabs = S.Eabs_map(b);
    Imax = S.Imax_map(b);

    Vp = Vth_pre(b);   Zp = Zth_pre(b);
    Vf = Vth_fault(b); Zf = Zth_fault(b);
    Vo = Vth_post(b);  Zo = Zth_post(b);

    Pe_pre      = pe_limited(delta_grid, Vp, Zp, Eabs, Zf_gfm, Imax);
    Pe_post     = pe_limited(delta_grid, Vo, Zo, Eabs, Zf_gfm, Imax);
    Pe_post_ref = pe_limited(delta_grid, Vo, Zo, Eabs, Zf_gfm, Iref);

    eq_pre   = find_equilibria_grid(delta_grid, Pe_pre,   Pm, fmin_tol, root_dedup);
    eq_post  = find_equilibria_grid(delta_grid, Pe_post,  Pm, fmin_tol, root_dedup);
    eq_postR = find_equilibria_grid(delta_grid, Pe_post_ref, Pm, fmin_tol, root_dedup);

    [delta0, ~, ~] = select_sep(eq_pre, delta_g(b), wrap2pi, circDist, wrap_dedup);
    if isnan(delta0), continue; end

    Pe_fault_0 = pe_limited(delta0, Vf, Zf, Eabs, Zf_gfm, Imax);
    acc0 = Pm - Pe_fault_0; dir = sign(acc0); if dir==0, dir=+1; end

    Pe_post_max = max(Pe_post); Pe_post_min = min(Pe_post);
    if ~(Pm <= Pe_post_max + eq_tol && Pm >= Pe_post_min - eq_tol)
        MhCL(ii) = 0; Delta(ii) = 0; DeltaRef(ii) = NaN; CCT_h(ii) = 0;
        continue;
    end

    [sepL, uepL, ~, ~, ~] = select_sep_uep(eq_post,  delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);
    [sepR, uepR, ~, ~, ~] = select_sep_uep(eq_postR, delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);
    if isnan(sepL) || isnan(uepL), continue; end

    Delta(ii) = abs(dir*(uepL - sepL));
    if ~isnan(sepR) && ~isnan(uepR)
        DeltaRef(ii) = abs(dir*(uepR - sepR));
        MhCL(ii) = Delta(ii) / max(1e-12, DeltaRef(ii));
    end

    alpha_u = dir*(uepL - sepL);
    while alpha_u <= 0, alpha_u = alpha_u + 2*pi; end
    alpha_grid = linspace(0, alpha_u, 10001).';
    delta_path = sepL + dir*alpha_grid;
    Pe_post_path = pe_limited(delta_path, Vo, Zo, Eabs, Zf_gfm, Imax);
    g = (Pe_post_path - Pm) * dir;
    U = cumtrapz(alpha_grid, g); U = U - U(1);
    Ecrit = U(end);
    Ecrit_h(ii) = Ecrit;

    odefun = @(tt,xx) fault_ode(tt, xx, Pm, S.H_gfm, S.D_gfm, S.omega_s, Vf, Zf, Eabs, Zf_gfm, Imax);
    odeopts = odeset('RelTol',opts.rtol,'AbsTol',opts.atol);
    [tt, X] = ode45(odefun, [0 opts.tc_max], [delta0; 0], odeopts);

    verdict_ok = @(tc) verdict_energy(tc, tt, X, sepL, dir, alpha_u, alpha_grid, U, Ecrit, S.H_gfm, S.omega_s);

    if verdict_ok(opts.tc_max)
        CCT_h(ii) = opts.tc_max;
        stable_to_tcmax_h(ii) = true;
        CCT_is_gt_h(ii) = true;
    else
        lo = 0.0; hi = opts.tc_max;
        if ~verdict_ok(lo)
            CCT_h(ii) = 0.0;
        else
            for it = 1:opts.max_it
                if (hi-lo) <= opts.tol_tc, break; end
                mid = 0.5*(lo+hi);
                if verdict_ok(mid), lo = mid; else, hi = mid; end
            end
            CCT_h(ii) = lo;
        end
    end
end

CCT_sys = min(CCT_h, [], 'omitnan');
Msys    = min(MhCL,  [], 'omitnan');

out.Vth_pre_map = Vth_pre; out.Zth_pre_map = Zth_pre;
out.Vth_fault_map = Vth_fault; out.Zth_fault_map = Zth_fault;
out.Vth_post_map = Vth_post; out.Zth_post_map = Zth_post;
out.MhCL = MhCL; out.Msys = Msys;
out.Delta = Delta; out.DeltaRef = DeltaRef;
out.CCT_h = CCT_h; out.CCT_sys = CCT_sys;
out.Ecrit_h = Ecrit_h;
out.stable_to_tcmax_h = stable_to_tcmax_h;
out.CCT_is_gt_h = CCT_is_gt_h;
out.gfm_buses = gfm_buses;
out.verdict_stable_to_breaker = (CCT_sys >= opts.t_clear_threshold);
out.verdict_post_eq_exists    = (CCT_sys > 0);
out.t_pipeline = toc(t_start);
end
