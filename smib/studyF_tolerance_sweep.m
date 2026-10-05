function studyF_tolerance_sweep()
% Finer mismatch sweep on the two-converter hub, behind Table XIII. Model
% and stability criterion are those of studyC_run, with
%   - 24 spread levels from 0.1 % to 164 % of the cluster mean,
%   - dispatch (F5) and internal EMF (F7) next to separation (F1), inertia
%     (F2), damping (F3) and current limit (F4),
%   - all four of H, D, Imax, Pref mismatched together, the larger dispatch
%     on the unit with the larger current limit (F6) or with the smaller
%     one (F8),
%   - a bisection tolerance of 1e-4 s on the time-domain CCT.
% The spread m of a setting p is (max p_k - min p_k) / mean p_k with the
% cluster total held constant: p_1,2 = (total/2)(1 +- m/2). A ratio r of
% studyC_run corresponds to m = 2(r-1)/(r+1), so r = 2, 4, 10 are on the
% grid.
% One row per case is appended to studyF_tolerance_sweep.csv, so a run
% that was interrupted continues where it stopped. studyF_envelope is
% called at the end.

t_total = tic;
CSV = 'studyF_tolerance_sweep.csv';

optsAgg = struct('t_hi',2.0,'tol',1e-3,'f0',60);
optsPer = struct('t_hi',2.0,'tol',1e-4,'T_post',10.0, ...
                 'instab_angle_rad', pi, 'rtol',1e-7,'atol',1e-9);
optsPer_ref = optsPer; optsPer_ref.tol = 1e-3;

fprintf('\nAggregation error on the two-converter hub, fine mismatch grid\n');
fprintf('bisection tolerance %.0e s (%.0e s in studyC_run)\n', ...
    optsPer.tol, optsPer_ref.tol);

% eight cases of studyC_run, at its tolerance, have to give the CCTs
% stored with the Table XII data
fprintf('\neight Table XII cases at the 1e-3 s tolerance:\n');
chk = replay_table12_cases(optsAgg, optsPer_ref);
if ~chk.pass
    error('studyF_tolerance_sweep:check', 'Table XII cases are not reproduced');
end
fprintf('  %d of %d reproduced exactly\n', chk.n_exact, chk.n);

% matched pair at the hub: no aggregation error by construction
fprintf('\nmatched pair at the hub:\n');
P0 = studyC_make_2GFM('Xa',0,'Xb',0,'Xc',0.7);
[agg0, CCT_agg0_EAC, M_agg0] = studyC_aggregate(P0, optsAgg);
[floor_d, floor_per, floor_agg] = one_case(P0, agg0, optsPer);
fprintf('  aggregate, energy criterion : CCT %.6f s, margin %.6f\n', ...
    CCT_agg0_EAC, M_agg0);
fprintf('  time domain : two converters %.9f s, aggregate %.9f s\n', floor_per, floor_agg);
fprintf('  error %.4f %%\n', 100*abs(floor_d));
if abs(CCT_agg0_EAC - 0.40625) >= 5e-4 || abs(M_agg0 - 0.7940) >= 5e-4
    error('studyF_tolerance_sweep:check', 'aggregate of the matched pair does not reproduce case 1');
end
if abs(floor_d) > 5e-3
    error('studyF_tolerance_sweep:check', 'matched pair: error above 0.5 %%');
end

cases = build_cases();
N = numel(cases);
fprintf('\n%d cases\n', N);

done = strings(0,1);
if isfile(CSV)
    Tprev = readtable(CSV, 'VariableNamingRule','preserve');
    done = string(Tprev.("case"));
    fprintf('  %d cases already in %s\n', numel(done), CSV);
end
n_new = 0;

for k = 1:N
    c = cases(k);
    if any(done == string(c.label)), continue; end

    P = c.makeP();
    feasible = true; note = 'ok';
    try
        [agg, CCT_agg_EAC, M_agg] = studyC_aggregate(P, optsAgg);
        [d, CCT_per, CCT_agg_TD, info] = one_case(P, agg, optsPer);
        E_h = agg.E_h; Imax_h = agg.Imax_h;
    catch ME
        d = NaN; CCT_per = NaN; CCT_agg_TD = NaN; CCT_agg_EAC = NaN;
        M_agg = NaN; E_h = NaN; Imax_h = NaN;
        info = struct('I_coh',NaN,'N_ode45',0);
        feasible = false; note = ME.message;
    end
    % d is also NaN when the aggregated pair has no operating point in the
    % time-domain model
    if isnan(d)
        feasible = false;
        if strcmp(note,'ok'), note = 'no operating point for the aggregated pair'; end
    end

    r.case        = c.label;
    r.family      = c.family;
    r.param       = c.param;
    r.mismatch    = c.mismatch;
    r.ratio       = c.ratio;
    r.Xsep        = P.Xa + P.Xb;
    r.Xc          = P.Xc;
    r.H1 = P.H1;   r.H2 = P.H2;
    r.D1 = P.D1;   r.D2 = P.D2;
    r.Pref1 = P.Pref1; r.Pref2 = P.Pref2;
    r.Imax1 = P.Imax1; r.Imax2 = P.Imax2;
    r.E_h         = E_h;
    r.Imax_h      = Imax_h;
    r.M_agg       = M_agg;
    r.CCT_agg_EAC = CCT_agg_EAC;
    r.CCT_agg_TD  = CCT_agg_TD;
    r.CCT_per_TD  = CCT_per;
    r.dCCT_rel    = d;
    r.feasible    = double(feasible);
    r.I_coh       = info.I_coh;
    r.N_ode45     = info.N_ode45;
    r.note        = note;

    Trow = struct2table(r, 'AsArray', true);
    if isfile(CSV)
        writetable(Trow, CSV, 'WriteMode','append', 'WriteVariableNames', false);
    else
        writetable(Trow, CSV);
    end
    n_new = n_new + 1;

    if mod(n_new, 10) == 0 || k == N
        fprintf('  %3d/%3d  %-26s  m = %.4f  error %+8.3f %%  (%.0f s)\n', ...
            k, N, c.label, c.mismatch, 100*d, toc(t_total));
    end
end

T = readtable(CSV, 'VariableNamingRule','preserve');
fprintf('\n%d cases in %s, %.1f s\n', height(T), CSV, toc(t_total));

studyF_envelope();
end

function chk = replay_table12_cases(optsAgg, optsPer_ref)
Tref = readtable(studyF_reference_csv(), 'VariableNamingRule','preserve');
ref_case = string(Tref.("case"));

probe = { ...
    'C1_X_0.010',          @() studyC_make_2GFM('Xa',0.005,'Xb',0.005); ...
    'C1_X_0.700',          @() studyC_make_2GFM('Xa',0.35,'Xb',0.35); ...
    'C2_H_2.000_X_0.200',  @() mk_pair('H', 2.0,  0.20); ...
    'C2_H_10.000_X_0.700', @() mk_pair('H', 10.0, 0.70); ...
    'C3_D_0.250_X_0.700',  @() mk_pair('D', 0.25, 0.70); ...
    'C3_D_4.000_X_0.050',  @() mk_pair('D', 4.0,  0.05); ...
    'C4_I_2.000_X_0.050',  @() mk_pair('I', 2.0,  0.05); ...
    'C4_I_0.500_X_0.050',  @() mk_pair('I', 0.5,  0.05)};

n = size(probe,1); n_exact = 0;
fprintf('%-22s %16s %16s %16s %16s  %s\n', ...
    'case','CCT_per','CCT_per ref','CCT_agg','CCT_agg ref','same');
det = struct([]);
for k = 1:n
    lbl = probe{k,1};
    P = probe{k,2}();
    [agg, ~, ~] = studyC_aggregate(P, optsAgg);
    [~, CCT_per, CCT_agg] = one_case(P, agg, optsPer_ref);
    j = find(ref_case == string(lbl), 1);
    per_ref = Tref.("CCT_per_TD")(j);
    agg_ref = Tref.("CCT_agg_TD")(j);
    ok = isequaln(CCT_per, per_ref) && isequaln(CCT_agg, agg_ref);
    n_exact = n_exact + ok;
    fprintf('%-22s %16.10f %16.10f %16.10f %16.10f  %s\n', ...
        lbl, CCT_per, per_ref, CCT_agg, agg_ref, string(ok));
    d.case = lbl; d.CCT_per = CCT_per; d.CCT_per_ref = per_ref;
    d.CCT_agg = CCT_agg; d.CCT_agg_ref = agg_ref; d.exact = ok;
    if isempty(det), det = d; else, det(end+1) = d; end
end
chk.n = n; chk.n_exact = n_exact; chk.pass = (n_exact == n); chk.detail = det;
end

function P = mk_pair(which, r, X)
% parameterisation of studyC_run: cluster total fixed, ratio r = p1/p2
switch which
    case 'H', s = 5;   v1 = s*r/(1+r); P = studyC_make_2GFM('Xa',X/2,'Xb',X/2,'H1',v1,'H2',s-v1);
    case 'D', s = 10;  v1 = s*r/(1+r); P = studyC_make_2GFM('Xa',X/2,'Xb',X/2,'D1',v1,'D2',s-v1);
    case 'I', s = 1.2; v1 = s*r/(1+r); P = studyC_make_2GFM('Xa',X/2,'Xb',X/2,'Imax1',v1,'Imax2',s-v1);
end
end

function P_agg = mk_aggTD(P, agg)
% aggregated cluster as two identical half-size units at the hub
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

function [d, CCT_per, CCT_agg, info] = one_case(P, agg, optsPer)
% relative CCT error of the aggregate, both sides by the time-domain
% criterion; negative when the aggregate gives the longer clearing time
[CCT_per, info] = studyC_per_converter_cct(P, optsPer);
CCT_agg         = studyC_per_converter_cct(mk_aggTD(P, agg), optsPer);
if isnan(CCT_per) || isnan(CCT_agg) || abs(CCT_agg) < 1e-12
    d = NaN;
else
    d = (CCT_per - CCT_agg) / CCT_agg;
end
end

function cases = build_cases()
m_grid = [0.001 0.002 0.005 0.01 0.02 0.03 0.05 0.075 0.10 0.125 0.15 ...
          0.20 0.25 0.30 0.35 0.40 0.50 0.60 2/3 0.80 1.00 1.20 1.40 18/11];
% EMF spread up to 0.30 only (1.15 pu against 0.85 pu)
m_grid_E = m_grid(m_grid <= 0.30);
X_subset = [0.05 0.20 0.70];
X_grid = [0 0.005 0.01 0.02 0.03 0.05 0.075 0.10 0.15 0.20 0.30 0.40 ...
          0.55 0.70 0.85 1.00];
m_joint = [0.01 0.02 0.03 0.05 0.10 0.15 0.20 0.25 0.30 0.40 0.50 2/3 1.00];

SUM_H = 5; SUM_D = 10; SUM_I = 1.2; SUM_P = 0.5; SUM_E = 2.0;

cases = struct([]);

for k = 1:numel(X_grid)
    X = X_grid(k);
    cases = add(cases, sprintf('F1_X_%06.4f', X), 'F1', 'Xsep', 0, NaN, ...
        @() studyC_make_2GFM('Xa', X/2, 'Xb', X/2));
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(m_grid)
        m = m_grid(k); [v1,v2] = split(SUM_H, m);
        cases = add(cases, sprintf('F2_H_m%06.4f_X%05.3f', m, X), 'F2', 'H', m, v1/v2, ...
            @() studyC_make_2GFM('Xa',X/2,'Xb',X/2,'H1',v1,'H2',v2));
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(m_grid)
        m = m_grid(k); [v1,v2] = split(SUM_D, m);
        cases = add(cases, sprintf('F3_D_m%06.4f_X%05.3f', m, X), 'F3', 'D', m, v1/v2, ...
            @() studyC_make_2GFM('Xa',X/2,'Xb',X/2,'D1',v1,'D2',v2));
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(m_grid)
        m = m_grid(k); [v1,v2] = split(SUM_I, m);
        cases = add(cases, sprintf('F4_I_m%06.4f_X%05.3f', m, X), 'F4', 'Imax', m, v1/v2, ...
            @() studyC_make_2GFM('Xa',X/2,'Xb',X/2,'Imax1',v1,'Imax2',v2));
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(m_grid)
        m = m_grid(k); [v1,v2] = split(SUM_P, m);
        cases = add(cases, sprintf('F5_P_m%06.4f_X%05.3f', m, X), 'F5', 'Pref', m, v1/v2, ...
            @() studyC_make_2GFM('Xa',X/2,'Xb',X/2,'Pref1',v1,'Pref2',v2));
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(m_joint)
        m = m_joint(k);
        [h1,h2] = split(SUM_H, m); [d1,d2] = split(SUM_D, m);
        [i1,i2] = split(SUM_I, m); [p1,p2] = split(SUM_P, m);
        cases = add(cases, sprintf('F6_all_m%06.4f_X%05.3f', m, X), 'F6', 'joint', m, NaN, ...
            @() studyC_make_2GFM('Xa',X/2,'Xb',X/2, ...
                'H1',h1,'H2',h2,'D1',d1,'D2',d2, ...
                'Imax1',i1,'Imax2',i2,'Pref1',p1,'Pref2',p2));
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(m_grid_E)
        m = m_grid_E(k); [v1,v2] = split(SUM_E, m);
        cases = add(cases, sprintf('F7_E_m%06.4f_X%05.3f', m, X), 'F7', 'E', m, v1/v2, ...
            @() studyC_make_2GFM('Xa',X/2,'Xb',X/2,'E1',v1,'E2',v2));
    end
end

for j = 1:numel(X_subset)
    X = X_subset(j);
    for k = 1:numel(m_joint)
        m = m_joint(k);
        [h1,h2] = split(SUM_H, m); [d1,d2] = split(SUM_D, m);
        [i1,i2] = split(SUM_I, m); [p1,p2] = split(SUM_P, m);
        cases = add(cases, sprintf('F8_adv_m%06.4f_X%05.3f', m, X), 'F8', 'joint_adverse', m, NaN, ...
            @() studyC_make_2GFM('Xa',X/2,'Xb',X/2, ...
                'H1',h1,'H2',h2,'D1',d1,'D2',d2, ...
                'Imax1',i2,'Imax2',i1, ...   % current limits swapped
                'Pref1',p1,'Pref2',p2));
    end
end
end

function [v1, v2] = split(s, m)
v1 = (s/2)*(1 + m/2);
v2 = (s/2)*(1 - m/2);
end

function cases = add(cases, label, family, param, m, r, fn)
c.label = label; c.family = family; c.param = param;
c.mismatch = m; c.ratio = r; c.makeP = fn;
if isempty(cases), cases = c; else, cases(end+1) = c; end
end
