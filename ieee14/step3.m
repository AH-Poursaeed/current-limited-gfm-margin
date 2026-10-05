% Step 3: current-limited power-angle curves of each GFM cluster, its
% post-fault SEP and UEP, and the margin M_CL = Delta / Delta_ref, eqs.
% (30)-(32) (Table IX). Writes STEP3_state.mat.

clear; close all;
load('STEP2_state.mat','S');

Nscan = 8001;
delta_grid = linspace(-2*pi, 2*pi, Nscan).';
fmin_tol   = 1e-4;
root_dedup = 1e-6;
wrap_dedup = 1e-6;
min_alpha  = 1e-7;
eq_tol     = 1e-6;

wrap2pi = @(x) mod(x, 2*pi);
circDist = @(a,b) abs(atan2(sin(a-b), cos(a-b)));

fprintf('\nStep 3: equilibria and margins\n');

MhCL = nan(numel(S.gfm_buses),1);
Delta = nan(numel(S.gfm_buses),1);
DeltaRef = nan(numel(S.gfm_buses),1);

for ii = 1:numel(S.gfm_buses)
    b = S.gfm_buses(ii);

    Pm   = S.Pm_map(b);
    Eabs = S.Eabs_map(b);

    Vth_pre   = S.Vth_pre_map(b);    Zth_pre   = S.Zth_pre_map(b);
    Vth_fault = S.Vth_fault_map(b);  Zth_fault = S.Zth_fault_map(b);
    Vth_post  = S.Vth_post_map(b);   Zth_post  = S.Zth_post_map(b);

    Imax  = S.Imax_map(b);
    Iref  = 1e6;    % no current limit, for Delta_ref

    Pe_pre   = pe_limited(delta_grid, Vth_pre,   Zth_pre,   Eabs, S.Zf_gfm, Imax);
    Pe_post  = pe_limited(delta_grid, Vth_post,  Zth_post,  Eabs, S.Zf_gfm, Imax);

    Pe_post_ref = pe_limited(delta_grid, Vth_post, Zth_post, Eabs, S.Zf_gfm, Iref);

    eq_pre   = find_equilibria_grid(delta_grid, Pe_pre,  Pm, fmin_tol, root_dedup);
    eq_post  = find_equilibria_grid(delta_grid, Pe_post, Pm, fmin_tol, root_dedup);
    eq_postR = find_equilibria_grid(delta_grid, Pe_post_ref, Pm, fmin_tol, root_dedup);

    [delta0, ~, pre_reason] = select_sep(eq_pre, S.delta_g_map(b), wrap2pi, circDist, wrap_dedup);
    if isnan(delta0)
        fprintf('\n-- bus %d --\n', b);
        fprintf('no pre-fault SEP (%s)\n', pre_reason);
        continue;
    end

    % direction of the swing from the acceleration at fault inception
    Pe_fault_0 = pe_limited(delta0, Vth_fault, Zth_fault, Eabs, S.Zf_gfm, Imax);
    acc0 = Pm - Pe_fault_0;
    dir = sign(acc0); if dir==0, dir=+1; end

    Pe_post_max = max(Pe_post);  Pe_post_min = min(Pe_post);
    if ~(Pm <= Pe_post_max + eq_tol && Pm >= Pe_post_min - eq_tol)
        fprintf('\n-- bus %d --\n', b);
        fprintf('no post-fault equilibrium: Pm=%.6f not in [%.6f, %.6f]\n', Pm, Pe_post_min, Pe_post_max);
        MhCL(ii)=0; Delta(ii)=0; DeltaRef(ii)=0;
        continue;
    end

    [sepL, uepL, ~, ~, rL] = select_sep_uep(eq_post,  delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);
    [sepR, uepR, ~, ~, rR] = select_sep_uep(eq_postR, delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);

    fprintf('\n-- bus %d --\n', b);
    fprintf('Pm=%.6f | dir=%+d | delta0=%.6f rad (%.3f deg)\n', Pm, dir, delta0, rad2deg(delta0));

    if isnan(sepL) || isnan(uepL)
        fprintf('limited curve: SEP/UEP not found (%s)\n', rL);
        continue;
    end
    if isnan(sepR) || isnan(uepR)
        fprintf('reference curve: SEP/UEP not found (%s)\n', rR);
        continue;
    end

    Delta(ii)    = dir*(uepL - sepL);
    DeltaRef(ii) = dir*(uepR - sepR);

    if Delta(ii) <= 0,    Delta(ii) = abs(Delta(ii)); end
    if DeltaRef(ii) <= 0, DeltaRef(ii) = abs(DeltaRef(ii)); end

    MhCL(ii) = Delta(ii) / max(1e-12, DeltaRef(ii));

    fprintf('limited  : sep=%.6f | uep=%.6f | Delta=%.6f rad (%.3f deg)\n', ...
        sepL, uepL, Delta(ii), rad2deg(Delta(ii)));
    fprintf('no limit : sep=%.6f | uep=%.6f | Delta_ref=%.6f rad (%.3f deg)\n', ...
        sepR, uepR, DeltaRef(ii), rad2deg(DeltaRef(ii)));
    fprintf('M_CL = %.6f\n', MhCL(ii));
end

Msys = min(MhCL,[],'omitnan');
fprintf('\nM_sys = %.6f\n', Msys);

save('STEP3_state.mat','S','MhCL','Delta','DeltaRef','Msys');
fprintf('saved STEP3_state.mat\n');
