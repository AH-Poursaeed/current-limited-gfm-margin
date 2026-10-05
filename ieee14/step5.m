% Step 5: one table with the results of all clusters (Tables IX and X):
% equilibria, margin, energy barrier, CCT, and the point on the post-fault
% SEP-UEP path where the unlimited current first reaches the limit.
% Writes ieee14_results.csv and ieee14_results.mat.

clear; close all;

load('STEP2_state.mat','S');
load('STEP4_state.mat','CCT_h','Ecrit_h','CCT_is_gt_h','CCT_report_h');

Nscan      = 8001;
delta_grid = linspace(-2*pi, 2*pi, Nscan).';
fmin_tol   = 1e-4;
root_dedup = 1e-6;
wrap_dedup = 1e-6;
min_alpha  = 1e-7;
eq_tol     = 1e-6;

wrap2pi = @(x) mod(x, 2*pi);
circDist = @(a,b) abs(atan2(sin(a-b), cos(a-b)));

Iref = 1e6;

rows = [];

fprintf('\nStep 5: results per cluster\n');

for ii = 1:numel(S.gfm_buses)
    b = S.gfm_buses(ii);

    Pm   = S.Pm_map(b);
    Eabs = S.Eabs_map(b);

    Vth_pre   = S.Vth_pre_map(b);    Zth_pre   = S.Zth_pre_map(b);
    Vth_post  = S.Vth_post_map(b);   Zth_post  = S.Zth_post_map(b);

    Imax = S.Imax_map(b);

    Pe_pre      = pe_limited(delta_grid, Vth_pre,  Zth_pre,  Eabs, S.Zf_gfm, Imax);
    Pe_post_L   = pe_limited(delta_grid, Vth_post, Zth_post, Eabs, S.Zf_gfm, Imax);
    Pe_post_ref = pe_limited(delta_grid, Vth_post, Zth_post, Eabs, S.Zf_gfm, Iref);

    eq_pre    = find_equilibria_grid(delta_grid, Pe_pre,    Pm, fmin_tol, root_dedup);
    eq_post_L = find_equilibria_grid(delta_grid, Pe_post_L, Pm, fmin_tol, root_dedup);
    eq_post_R = find_equilibria_grid(delta_grid, Pe_post_ref, Pm, fmin_tol, root_dedup);

    [delta0, ~, pre_reason] = select_sep(eq_pre, S.delta_g_map(b), wrap2pi, circDist, wrap_dedup);
    if isnan(delta0)
        warning('bus %d: no pre-fault SEP (%s), skipped', b, pre_reason);
        continue;
    end

    Pe_fault_0 = pe_limited(delta0, S.Vth_fault_map(b), S.Zth_fault_map(b), Eabs, S.Zf_gfm, Imax);
    dir = sign(Pm - Pe_fault_0); if dir==0, dir=+1; end

    Pe_post_max = max(Pe_post_L);  Pe_post_min = min(Pe_post_L);
    eq_exists = (Pm <= Pe_post_max + eq_tol) && (Pm >= Pe_post_min - eq_tol);

    sepL=NaN; uepL=NaN; sepR=NaN; uepR=NaN; reasonL=""; reasonR="";
    if eq_exists
        [sepL, uepL, ~, ~, reasonL] = select_sep_uep(eq_post_L, delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);
        [sepR, uepR, ~, ~, reasonR] = select_sep_uep(eq_post_R, delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);
    end
    uep_wrapped = NaN;
    if ~isnan(uepL)
        uep_wrapped = mod(uepL, 2*pi);
    end
    Delta    = NaN; DeltaRef = NaN; MhCL = NaN;
    alpha_u  = NaN;

    if eq_exists && ~isnan(sepL) && ~isnan(uepL) && ~isnan(sepR) && ~isnan(uepR)
        alpha_u  = dir*(uepL - sepL); while alpha_u <= 0, alpha_u = alpha_u + 2*pi; end
        Delta    = alpha_u;

        DeltaRef = dir*(uepR - sepR); while DeltaRef <= 0, DeltaRef = DeltaRef + 2*pi; end
        MhCL     = Delta / max(1e-12, DeltaRef);
    end

    % first angle on the SEP-UEP path at which the unlimited current
    % exceeds Imax (alphaCL_on), its share of the path (MI), and the
    % terminal voltage there (Vcr)
    alphaCL_on = NaN; MI = NaN; Vcr = NaN;

    if ~isnan(alpha_u)
        Nalpha = 8001;
        alpha_grid = linspace(0, alpha_u, Nalpha).';
        delta_path = sepL + dir*alpha_grid;

        [Iu_abs, Vbus_abs] = unlimited_current_voltage(delta_path, S.Vth_post_map(b), S.Zth_post_map(b), Eabs, S.Zf_gfm);

        g = Iu_abs - Imax;
        idx = find(g(1:end-1) <= 0 & g(2:end) > 0, 1, 'first');

        if isempty(idx)
            alphaCL_on = alpha_u;
            MI = 1.0;
            Vcr = NaN;
        else
            a1 = alpha_grid(idx); a2 = alpha_grid(idx+1);
            g1 = g(idx);          g2 = g(idx+1);
            if abs(g2-g1) < 1e-15
                alphaCL_on = 0.5*(a1+a2);
            else
                alphaCL_on = a1 - g1*(a2-a1)/(g2-g1);
            end
            alphaCL_on = max(0, min(alpha_u, alphaCL_on));
            MI = alphaCL_on / max(1e-12, alpha_u);

            Vcr = interp1(alpha_grid, Vbus_abs, alphaCL_on, 'linear');
        end
    end

    Ecrit = Ecrit_h(ii);
    CCT   = CCT_h(ii);

    r = table();
    r.bus = b;
    r.Pm_pu = Pm;
    r.dir = dir;
    r.delta0 = delta0;

    r.sep_post = sepL;
    r.uep_post_unwrapped = uepL;
    r.Delta_rad = Delta;
    r.DeltaRef_rad = DeltaRef;
    r.MhCL = MhCL;

    r.alpha_u = alpha_u;
    r.alphaCL_on = alphaCL_on;
    r.MI = MI;
    r.Vcr_pu = Vcr;

    r.Ecrit = Ecrit;
    r.CCT_s = CCT;

    r.uep_post_wrapped = uep_wrapped;
    r.CCT_is_gt  = CCT_is_gt_h(ii);
    r.CCT_report = CCT_report_h(ii);

    r.eq_exists = eq_exists;
    r.reason_postL = string(reasonL);
    r.reason_postR = string(reasonR);

    rows = [rows; r];
end

Msys = min(rows.MhCL,[],'omitnan');
[CCTsys_calc, idxLim] = min(rows.CCT_s,[],'omitnan');
limBus = rows.bus(idxLim);

fprintf('\nM_sys   = %.6f\n', Msys);
fprintf('CCT_sys = %.6f s (bus %d)\n', CCTsys_calc, limBus);

disp(rows);

writetable(rows, 'ieee14_results.csv');
save('ieee14_results.mat','rows','Msys','CCTsys_calc','limBus','S');

fprintf('\nsaved ieee14_results.csv and ieee14_results.mat\n');

function [Iu_abs, Vpcc_abs] = unlimited_current_voltage(delta, Vth, Zth, Eabs, Zf)
% magnitudes of the current and terminal voltage without the limiter
    theta_th = angle(Vth);
    delta_col = delta(:);

    E  = Eabs .* exp(1j*(theta_th + delta_col));
    Iu = (E - Vth) ./ (Zth + Zf);

    Vpcc = Vth + Zth .* Iu;

    Iu_abs   = abs(Iu);
    Vpcc_abs = abs(Vpcc);
end
