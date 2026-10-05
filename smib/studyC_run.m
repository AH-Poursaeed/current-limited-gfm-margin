function studyC_run()
% Aggregation error on a two-converter hub (Table XII). For every case the
% CCT of the aggregated cluster is compared with the CCT of the two
% converters simulated one by one, both judged by the same time-domain
% criterion. Case families:
%   C1  electrical separation, parameters matched
%   C2  inertia ratio          C3  damping ratio
%   C4  current-limit ratio    (each at three separations)
%   C5  inertia and current-limit ratios together
%   C6  the cases of C1 once more
% Writes studyC_aggregation_error.csv.

t_total = tic;

optsCCT_agg = struct('t_hi',2.0,'tol',1e-3,'f0',60);
optsCCT_per = struct('t_hi',2.0,'tol',1e-3,'T_post',10.0, ...
                     'instab_angle_rad', pi, ...
                     'rtol',1e-7,'atol',1e-9);

fprintf('\nAggregation error on the two-converter hub\n');

% matched pair at the hub itself. Aggregation is exact here, so the two
% time-domain CCTs must agree and the energy-based CCT must be that of
% case 1 in Table II.
P0 = studyC_make_2GFM('Xa', 0, 'Xb', 0, 'Xc', 0.7);
[agg0, CCT_agg0_EAC, M_agg0] = studyC_aggregate(P0, optsCCT_agg);
CCT_per_TD_0 = studyC_per_converter_cct(P0, optsCCT_per);
P0_aggTD = build_agg_TD_P(P0, agg0);
CCT_agg_TD_0 = studyC_per_converter_cct(P0_aggTD, optsCCT_per);

d0_rel = (CCT_per_TD_0 - CCT_agg_TD_0) / max(1e-12, CCT_agg_TD_0);
EAC_vs_TD_pct = 100*(CCT_per_TD_0 - CCT_agg0_EAC)/max(1e-12, CCT_agg0_EAC);

fprintf('matched pair at the hub:\n');
fprintf('  CCT, aggregate, energy criterion      : %.6f s\n', CCT_agg0_EAC);
fprintf('  CCT, aggregate, time-domain criterion : %.6f s\n', CCT_agg_TD_0);
fprintf('  CCT, two converters, time domain      : %.6f s\n', CCT_per_TD_0);
fprintf('  difference between the last two       : %+.4f %%\n', 100*d0_rel);
fprintf('  time domain against energy criterion  : %+.2f %%\n', EAC_vs_TD_pct);
fprintf('  aggregate margin %.6f, cluster EMF %.6f pu\n', M_agg0, agg0.E_h);

if abs(CCT_agg0_EAC - 0.40625) >= 5e-4 || abs(M_agg0 - 0.7940) >= 5e-4
    error('studyC_run:check', 'aggregate of the matched pair does not reproduce case 1');
end
if abs(d0_rel) > 5e-3
    error('studyC_run:check', 'matched pair: the two CCTs differ by %.4f %%', 100*d0_rel);
end

cases = build_cases();
N = numel(cases);
fprintf('\n%d cases:\n', N);

rows = repmat(struct(), N, 1);

for k = 1:N
    c = cases(k);
    P = c.makeP();
    try
        [agg, CCT_agg_EAC, M_agg] = studyC_aggregate(P, optsCCT_agg);
        [CCT_per_TD, info_per]    = studyC_per_converter_cct(P, optsCCT_per);
        P_aggTD = build_agg_TD_P(P, agg);
        CCT_agg_TD = studyC_per_converter_cct(P_aggTD, optsCCT_per);
    catch ME
        % cases without an operating point end up here
        fprintf('  %s: %s\n', c.label, ME.message);
        CCT_agg_EAC = NaN; M_agg = NaN;
        CCT_per_TD = NaN;  CCT_agg_TD = NaN;
        info_per = struct('I_coh',NaN,'peak_dw_diff',NaN,'N_ode45',NaN,'N_bisect',NaN);
    end
    % negative: the aggregate gives the longer clearing time
    if isnan(CCT_per_TD) || isnan(CCT_agg_TD) || abs(CCT_agg_TD) < 1e-12
        dCCT_rel = NaN;
    else
        dCCT_rel = (CCT_per_TD - CCT_agg_TD) / CCT_agg_TD;
    end
    rows(k).case        = c.label;
    rows(k).family      = c.family;
    rows(k).description = c.description;
    rows(k).Xa          = P.Xa;
    rows(k).Xb          = P.Xb;
    rows(k).Xsep        = P.Xa + P.Xb;
    rows(k).Xc          = P.Xc;
    rows(k).H1          = P.H1; rows(k).H2 = P.H2;
    rows(k).D1          = P.D1; rows(k).D2 = P.D2;
    rows(k).Pref1       = P.Pref1; rows(k).Pref2 = P.Pref2;
    rows(k).Imax1       = P.Imax1; rows(k).Imax2 = P.Imax2;
    rows(k).ratio_H     = P.H1/P.H2;
    rows(k).ratio_D     = P.D1/P.D2;
    rows(k).ratio_Imax  = P.Imax1/P.Imax2;
    rows(k).CCT_agg_EAC = CCT_agg_EAC;
    rows(k).CCT_agg_TD  = CCT_agg_TD;
    rows(k).CCT_per_TD  = CCT_per_TD;
    rows(k).dCCT_rel    = dCCT_rel;
    rows(k).M_agg       = M_agg;
    rows(k).I_coh       = info_per.I_coh;
    rows(k).peak_dw_diff= info_per.peak_dw_diff;
    rows(k).N_ode45     = info_per.N_ode45;
    if mod(k, 10) == 0
        fprintf('  %3d/%d  %-22s  error %+.3f %%\n', k, N, c.label, 100*dCCT_rel);
    end
end

T = struct2table(rows);
writetable(T, 'studyC_aggregation_error.csv');

feasible = ~isnan(T.dCCT_rel);
T_feas = T(feasible,:);
abs_dCCT_pct = 100*abs(T_feas.dCCT_rel);
[mx_abs_dCCT, idx_mx_local] = max(abs_dCCT_pct);
idx_mx = find(feasible);
idx_mx = idx_mx(idx_mx_local);
n_infeasible = sum(~feasible);

within_1   = abs_dCCT_pct <= 1;
within_5   = abs_dCCT_pct <= 5;
within_10  = abs_dCCT_pct <= 10;

fprintf('\n%d cases: %d with an operating point, %d without (%.1f s)\n', ...
    height(T), height(T_feas), n_infeasible, toc(t_total));
fprintf('largest |error| %.3f %% in %s\n', mx_abs_dCCT, T.case{idx_mx});
fprintf('within 1 %%: %d, within 5 %%: %d, within 10 %%: %d\n', ...
    sum(within_1), sum(within_5), sum(within_10));

fams = unique(T.family, 'stable');
for ii = 1:numel(fams)
    f = fams(ii);
    idx_f = strcmp(T.family, f);
    n_inf_f = sum(~feasible & idx_f);
    idx = idx_f & feasible;
    if ~any(idx)
        fprintf('%s : %2d cases, none with an operating point\n', f, sum(idx_f));
        continue;
    end
    raw_f   = 100*abs(T.dCCT_rel(idx));
    fprintf('%s : %2d cases (%d without operating point), max |error| %5.2f %%\n', ...
        f, sum(idx_f), n_inf_f, max(raw_f));
end

end

function P_agg = build_agg_TD_P(P, agg)
% the aggregated cluster written as two identical half-size units at the
% hub, so that studyC_per_converter_cct can be applied to it
P_agg = studyC_make_2GFM( ...
    'H1', agg.H_h/2, 'H2', agg.H_h/2, ...
    'D1', agg.D_h/2, 'D2', agg.D_h/2, ...
    'Pref1', agg.Pref_h/2, 'Pref2', agg.Pref_h/2, ...
    'Qref1', agg.Qref_h/2, 'Qref2', agg.Qref_h/2, ...
    'Imax1', agg.Imax_h/2, 'Imax2', agg.Imax_h/2, ...
    'E1', agg.E_h, 'E2', agg.E_h, ...
    'Vth', P.Vth, 'theta_th', P.theta_th, ...
    'Xa', 0, 'Xb', 0, 'Xc', P.Xc, ...
    'kV', P.kV, 'kX', P.kX, 'f0', P.f0);
end

function cases = build_cases()
% X is the total reactance between the two converters; ratios are unit 1
% over unit 2 with the cluster total held constant
X_grid = [0.01 0.02 0.05 0.10 0.20 0.40 0.70 1.00];
ratio_H = [0.1 0.25 0.5 1 2 4 10];
ratio_D = ratio_H;
ratio_I = ratio_H;

cases = struct([]);

for k = 1:numel(X_grid)
    X = X_grid(k);
    cases(end+1).label = sprintf('C1_X_%05.3f', X);
    cases(end).family = "C1";
    cases(end).description = sprintf('separation X_a=X_b=%.4f, homogeneous params', X/2);
    cases(end).makeP = @() studyC_make_2GFM('Xa', X/2, 'Xb', X/2);
end

X_subset = [0.05 0.20 0.70];
for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(ratio_H)
        r = ratio_H(k);
        Hsum = 5;
        H1 = Hsum * r/(1+r);   H2 = Hsum - H1;
        cases(end+1).label = sprintf('C2_H_%05.3f_X_%05.3f', r, X);
        cases(end).family = "C2";
        cases(end).description = sprintf('H_1/H_2 = %.3f, X = %.3f', r, X);
        cases(end).makeP = @() studyC_make_2GFM('Xa', X/2, 'Xb', X/2, 'H1', H1, 'H2', H2);
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(ratio_D)
        r = ratio_D(k);
        Dsum = 10;
        D1 = Dsum * r/(1+r);   D2 = Dsum - D1;
        cases(end+1).label = sprintf('C3_D_%05.3f_X_%05.3f', r, X);
        cases(end).family = "C3";
        cases(end).description = sprintf('D_1/D_2 = %.3f, X = %.3f', r, X);
        cases(end).makeP = @() studyC_make_2GFM('Xa', X/2, 'Xb', X/2, 'D1', D1, 'D2', D2);
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(ratio_I)
        r = ratio_I(k);
        Isum = 1.2;
        I1 = Isum * r/(1+r);   I2 = Isum - I1;
        cases(end+1).label = sprintf('C4_I_%05.3f_X_%05.3f', r, X);
        cases(end).family = "C4";
        cases(end).description = sprintf('Imax_1/Imax_2 = %.3f, X = %.3f', r, X);
        cases(end).makeP = @() studyC_make_2GFM('Xa', X/2, 'Xb', X/2, 'Imax1', I1, 'Imax2', I2);
    end
end

joint_combos = struct( ...
    'r_H',  num2cell([10, 0.1, 4, 0.25]), ...
    'r_I',  num2cell([4, 0.25, 10, 0.1]));
for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(joint_combos)
        rH = joint_combos(k).r_H;
        rI = joint_combos(k).r_I;
        Hsum = 5; H1 = Hsum*rH/(1+rH); H2 = Hsum - H1;
        Isum = 1.2; I1 = Isum*rI/(1+rI); I2 = Isum - I1;
        cases(end+1).label = sprintf('C5_jH%05.3f_jI%05.3f_X%05.3f', rH, rI, X);
        cases(end).family = "C5";
        cases(end).description = sprintf('joint: H ratio %.2f, Imax ratio %.2f, X=%.3f', rH, rI, X);
        cases(end).makeP = @() studyC_make_2GFM('Xa', X/2, 'Xb', X/2, ...
            'H1', H1, 'H2', H2, 'Imax1', I1, 'Imax2', I2);
    end
end

for k = 1:numel(X_grid)
    X = X_grid(k);
    cases(end+1).label = sprintf('C6_homog_X_%05.3f', X);
    cases(end).family = "C6";
    cases(end).description = sprintf('homogeneous baseline at X=%.4f (reference)', X);
    cases(end).makeP = @() studyC_make_2GFM('Xa', X/2, 'Xb', X/2);
end
end
