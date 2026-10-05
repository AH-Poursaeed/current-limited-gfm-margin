% Step 4: post-fault energy barrier and CCT of each cluster, eqs. (33)-(38)
% (Table X). The fault-on trajectory is integrated once up to tc_max, the
% verdict for a clearing time is read from it, and the CCT follows by
% bisection. A cluster still stable at tc_max is reported as "> tc_max".
% Writes STEP4_state.mat.

clear; close all;

load('STEP3_state.mat','S');

tc_max  = 10.0;
rtol    = 1e-7;
atol    = 1e-9;
tol_tc  = 1e-6;
max_it  = 70;

Nscan     = 8001;
delta_grid= linspace(-2*pi, 2*pi, Nscan).';
fmin_tol  = 1e-4;
root_dedup= 1e-6;
wrap_dedup= 1e-6;
min_alpha = 1e-7;
eq_tol    = 1e-6;

wrap2pi = @(x) mod(x, 2*pi);
circDist = @(a,b) abs(atan2(sin(a-b), cos(a-b)));

fprintf('\nStep 4: energy barrier and CCT\n');

CCT_h = nan(numel(S.gfm_buses),1);
Ecrit_h = nan(numel(S.gfm_buses),1);

stable_to_tcmax_h = false(numel(S.gfm_buses),1);
CCT_is_gt_h       = false(numel(S.gfm_buses),1);
CCT_report_h      = strings(numel(S.gfm_buses),1);

for ii = 1:numel(S.gfm_buses)
    b = S.gfm_buses(ii);

    Pm   = S.Pm_map(b);
    Eabs = S.Eabs_map(b);

    Vth_pre   = S.Vth_pre_map(b);    Zth_pre   = S.Zth_pre_map(b);
    Vth_fault = S.Vth_fault_map(b);  Zth_fault = S.Zth_fault_map(b);
    Vth_post  = S.Vth_post_map(b);   Zth_post  = S.Zth_post_map(b);

    Imax = S.Imax_map(b);

    Pe_pre  = pe_limited(delta_grid, Vth_pre,  Zth_pre,  Eabs, S.Zf_gfm, Imax);
    Pe_post = pe_limited(delta_grid, Vth_post, Zth_post, Eabs, S.Zf_gfm, Imax);

    eq_pre  = find_equilibria_grid(delta_grid, Pe_pre,  Pm, fmin_tol, root_dedup);
    eq_post = find_equilibria_grid(delta_grid, Pe_post, Pm, fmin_tol, root_dedup);

    [delta0, ~, pre_reason] = select_sep(eq_pre, S.delta_g_map(b), wrap2pi, circDist, wrap_dedup);
    if isnan(delta0)
        fprintf('\n-- bus %d --\n', b);
        fprintf('no pre-fault SEP (%s)\n', pre_reason);
        CCT_h(ii) = NaN;
        continue;
    end

    Pe_fault_0 = pe_limited(delta0, Vth_fault, Zth_fault, Eabs, S.Zf_gfm, Imax);
    acc0 = Pm - Pe_fault_0;
    dir = sign(acc0); if dir==0, dir=+1; end

    Pe_post_max = max(Pe_post);  Pe_post_min = min(Pe_post);
    if ~(Pm <= Pe_post_max + eq_tol && Pm >= Pe_post_min - eq_tol)
        fprintf('\n-- bus %d --\n', b);
        fprintf('no post-fault equilibrium: Pm=%.6f not in [%.6f, %.6f], CCT = 0\n', Pm, Pe_post_min, Pe_post_max);
        CCT_h(ii) = 0.0;
        continue;
    end

    [sep, uep_unwrapped, ~, ~, rPost] = select_sep_uep(eq_post, delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);
    if isnan(sep) || isnan(uep_unwrapped)
        fprintf('\n-- bus %d --\n', b);
        fprintf('post-fault SEP/UEP not found (%s)\n', rPost);
        CCT_h(ii) = NaN;
        continue;
    end

    % potential energy along alpha, the angle from the SEP in the swing
    % direction; its value at the UEP is the barrier, eqs. (33) and (35)
    alpha_u = dir*(uep_unwrapped - sep);
    while alpha_u <= 0
        alpha_u = alpha_u + 2*pi;
    end

    Nalpha = 10001;
    alpha_grid = linspace(0, alpha_u, Nalpha).';
    delta_path = sep + dir*alpha_grid;

    Pe_post_path = pe_limited(delta_path, Vth_post, Zth_post, Eabs, S.Zf_gfm, Imax);

    g = (Pe_post_path - Pm) * dir;
    U = cumtrapz(alpha_grid, g);
    U = U - U(1);
    Ecrit = U(end);

    Ecrit_h(ii) = Ecrit;

    % fault-on trajectory from the pre-fault SEP
    odefun = @(tt,xx) fault_ode(tt, xx, Pm, S.H_gfm, S.D_gfm, S.omega_s, Vth_fault, Zth_fault, Eabs, S.Zf_gfm, Imax);
    opts = odeset('RelTol',rtol,'AbsTol',atol);
    [t, X] = ode45(odefun, [0 tc_max], [delta0; 0], opts);

    verdict_ok = @(tc) verdict_energy(tc, t, X, sep, dir, alpha_u, alpha_grid, U, Ecrit, S.H_gfm, S.omega_s);

    ok_max = verdict_ok(tc_max);
    if ok_max
        CCT = tc_max;
        stable_to_max = 1;
    else
        stable_to_max = 0;
        lo = 0.0; hi = tc_max;
        if ~verdict_ok(lo)
            CCT = 0.0;
        else
            for it = 1:max_it
                if (hi - lo) <= tol_tc, break; end
                mid = 0.5*(lo+hi);
                if verdict_ok(mid)
                    lo = mid;
                else
                    hi = mid;
                end
            end
            CCT = lo;
        end
    end

    CCT_h(ii) = CCT;
    CCT_is_gt = stable_to_max;
    if CCT_is_gt
        CCT_report = ">" + string(tc_max);
    else
        CCT_report = string(CCT);
    end

    stable_to_tcmax_h(ii) = stable_to_max;
    CCT_is_gt_h(ii)       = CCT_is_gt;
    CCT_report_h(ii)      = CCT_report;

    fprintf('\n-- bus %d --\n', b);
    fprintf('Pm=%.6f | dir=%+d | delta0=%.6f rad\n', Pm, dir, delta0);
    fprintf('SEP_post=%.6f rad | UEP_post=%.6f rad | alpha_u=%.6f rad\n', sep, uep_unwrapped, alpha_u);
    fprintf('Ecrit=%.6e | CCT=%.6f s | stable_to_tcmax=%d\n', Ecrit, CCT, stable_to_max);
end

CCT_sys = min(CCT_h,[],'omitnan');
[~,idx_lim] = min(CCT_h,[],'omitnan');
lim_bus = S.gfm_buses(idx_lim);

fprintf('\nCCT_sys = %.6f s (bus %d)\n', CCT_sys, lim_bus);

save('STEP4_state.mat','S','CCT_h','CCT_sys','Ecrit_h', ...
     'stable_to_tcmax_h','CCT_is_gt_h','CCT_report_h','tc_max');
fprintf('saved STEP4_state.mat\n');
