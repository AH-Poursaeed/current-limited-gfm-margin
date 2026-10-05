function varargout = check_paper_values(src)
% Compares computed results with the numbers printed in the paper.
%   check_paper_values          results written by the scripts (smib/, ieee14/)
%   check_paper_values data     the result files in data/
%   n = check_paper_values(...) also returns the number of failed checks
% One line per check with the largest deviation and the tolerance, which is
% half a unit of the last printed digit unless noted.

if nargin < 1, src = 'results'; end
root = fileparts(mfilename('fullpath'));
switch lower(src)
    case 'results'
        dS = fullfile(root, 'smib');
        dN = fullfile(root, 'ieee14');
    case 'data'
        dS = fullfile(root, 'data', 'smib');
        dN = fullfile(root, 'data', 'ieee14');
    otherwise
        error('check_paper_values: use ''results'' or ''data''.');
end
rd = @(d, f) readtable(fullfile(d, f), 'VariableNamingRule', 'preserve');

nfail = 0;
ncheck = 0;
fprintf('\n%-58s %12s %10s\n', 'check', 'deviation', 'tolerance');

%% Table I: SEP-UEP separation and margin, cases 1-4
T = rd(dS, 'ias_ccpac_margins.csv');
num('Table I    Delta_ref (deg)', T.DeltaRef_deg, [139.03; 151.10; 153.08; 139.03], 5e-3);
num('Table I    Delta (deg)',     T.DeltaLim_deg, [110.40; 116.34; 125.88; 130.57], 5e-3);
num('Table I    M_CL',            T.M_CL, [0.79404; 0.76996; 0.82235; 0.93912], 5e-6);

%% Table II: energy barrier and CCT
% V_crit is given to five significant digits (1.2433 for case 4)
T = rd(dS, 'ias_ccpac_cct.csv');
num('Table II   V_crit (pu)', T.Vcrit,   [0.62515; 0.74355; 0.98089; 1.24330], 5e-5);
num('Table II   CCT (s)',     T.CCT_sec, [0.40625; 0.40625; 0.44434; 0.53516], 5e-6);

%% Table III: CCT under the three limiter modes
T = rd(dS, 'studyD_main_sweep.csv');
modes = {'circular', 'd-priority', 'q-priority'};
ref3 = [0.40625 0.40625 0.44434 0.53516
        0.56445 0.63379 0.64453 0.58008
        0.29199 0.21582 0.27441 0.45605];
cct3 = nan(3, 4); m3 = nan(3, 4);
for m = 1:3
    for c = 1:4
        k = strcmp(T.Mode, modes{m}) & T.Case == c;
        cct3(m, c) = T.CCT_s(k);
        m3(m, c)   = T.MhCL(k);
    end
end
num('Table III  CCT (s), 12 entries', cct3, ref3, 5e-6);
shift = 100 * (mean(cct3(2:3, :), 2) - mean(cct3(1, :))) / mean(cct3(1, :));
num('Sec. IV.A  mean CCT shift d, q (%)', shift, [35.2; -30.9], 0.05);
istrue('Sec. IV.A  margin ordered d > circular > q, 4 cases', ...
    all(m3(2, :) > m3(1, :) & m3(1, :) > m3(3, :)));
istrue('Sec. IV.A  margin equals one under d-priority', all(m3(2, :) == 1));
% the margin does not depend on the fault depth; what the sweep can show
% is the CCT ordering (the three modes tie where all are stable up to 2 s)
T = rd(dS, 'studyD_kV_sweep.csv');
istrue('Sec. IV.A  CCT ordered d >= circular >= q, kV sweep', ...
    all(T.CCT_d_priority >= T.CCT_circular & T.CCT_circular >= T.CCT_q_priority));

%% Section IV.A: ordering against R/X
T = rd(dS, 'studyE_resistive_sweep.csv');
key = strcat(T.param, '|', string(T.rho), '|', string(T.case_id));
u = unique(key(T.rho <= 0.25));
nScored = 0; nHold = 0;
for i = 1:numel(u)
    s = key == u(i);
    if any(T.feasible(s) == 0) || any(T.censored(s) == 1) || any(T.limiter_active(s) == 0)
        continue;
    end
    nScored = nScored + 1;
    cc = T.CCT_s(s & strcmp(T.mode, 'circular'));
    cd_ = T.CCT_s(s & strcmp(T.mode, 'd-priority'));
    cq = T.CCT_s(s & strcmp(T.mode, 'q-priority'));
    nHold = nHold + (cd_ > cc && cc > cq);
end
num('Sec. IV.A  sweep points with R/X <= 0.25', numel(u), 64, 0);
num('Sec. IV.A  of which resolved', nScored, 60, 0);
num('Sec. IV.A  ordering holds in', nHold, 60, 0);
% smallest R/X at which q-priority is not below d-priority
uu = unique(key);
rhoRev = inf;
for i = 1:numel(uu)
    s = key == uu(i);
    if any(T.feasible(s) == 0) || any(T.censored(s) == 1) || any(T.limiter_active(s) == 0)
        continue;
    end
    if T.CCT_s(s & strcmp(T.mode, 'q-priority')) >= T.CCT_s(s & strcmp(T.mode, 'd-priority'))
        rhoRev = min(rhoRev, T.rho(find(s, 1)));
    end
end
num('Sec. IV.A  X/R at the first reversal ("about 1.4")', 1 / rhoRev, 1.4, 0.05);

%% Table XIV: CCT under the two frames of the priority split
T = rd(dS, 'studyE_frame_check.csv');
T = T(T.rho == 0, :);
frames = {'thevenin', 'converter'};
ref14 = [0.4065 0.5651 0.2927; 0.4065 0.2501 0.5217
         0.4070 0.6339 0.2166; 0.4070 0.2367 0.5116
         0.4447 0.6446 0.2744; 0.4447 0.2679 0.5115
         0.5359 0.5801 0.4565; 0.5359 0.4870 0.6240];
got = nan(8, 3);
for c = 1:4
    for f = 1:2
        for m = 1:3
            k = T.case_id == c & strcmp(T.frame, frames{f}) & strcmp(T.mode, modes{m});
            got(2*(c-1) + f, m) = T.CCT_s(k);
        end
    end
end
num('Table XIV  CCT (s), 24 entries', got, ref14, 5e-5);
num('Table XIV  Thevenin rows against Table III (< 1e-3 s)', got(1:2:end, :), ref3.', 1e-3);

%% Table XII: aggregation error by source of heterogeneity
T = rd(dS, 'studyC_aggregation_error.csv');
fam = string(T.family);
e = 100 * T.dCCT_rel;
num('Table XII  separation, max error (%)',    max(abs(e(fam == "C1"))), 8.33, 5e-3);
num('Table XII  damping, max error (%)',       max(abs(e(fam == "C3"))), 6.93, 5e-3);
num('Table XII  inertia, max error (%)',       max(abs(e(fam == "C2"))), 20.82, 5e-3);
num('Table XII  current limit, max error (%)', max(abs(e(fam == "C4"))), 33.41, 5e-3);
istrue('Table XII  inertia + current limit: no operating point', all(isnan(e(fam == "C5"))));
num('Table XII  smallest error, matched pair (%)', min(abs(e(~isnan(e)))), 0.23, 5e-3);
istrue('Sec. IV.C  aggregate CCT longer in every case with a value', all(e(~isnan(e)) < 0));

%% Table XIII: admissible spread at 0.05 pu separation, 5 % bound
T = rd(dS, 'studyF_admissibility_envelope.csv');
k = T.bound_pct == 5 & abs(T.Xsep - 0.05) < 1e-9;
f13  = {'F4', 'F5', 'F2', 'F3', 'F8'};
nm13 = {'current limit', 'dispatch', 'inertia', 'damping', 'all four'};
ref13 = [0.075 4.86; 0.10 4.88; 0.15 4.45; 18/11 4.83; 0.03 4.52];
for i = 1:5
    j = k & strcmp(T.family, f13{i});
    num(sprintf('Table XIII %s, spread', nm13{i}), T.max_admissible(j), ref13(i, 1), 5e-4);
    num(sprintf('Table XIII %s, error (%%)', nm13{i}), T.err_at_max(j), ref13(i, 2), 5e-3);
end
T = rd(dS, 'studyF_tolerance_sweep.csv');
num('Sec. IV.C  cases in the fine sweep', height(T), 424, 0);
e = 100 * T.dCCT_rel;
f1 = strcmp(T.family, 'F1');
num('Sec. II.A  separation 0.15 pu alone (%)', abs(e(f1 & abs(T.Xsep - 0.15) < 1e-9)), 0.83, 5e-3);
num('Sec. II.A  separation 0.20 pu alone (%)', abs(e(f1 & abs(T.Xsep - 0.20) < 1e-9)), 1.19, 5e-3);
num('Sec. IV.C  largest case with shorter aggregate CCT (%)', max(e(T.feasible == 1)), 1.19, 5e-3);
% the two sweeps at the grid points they share (ratios 2, 4 and 10)
Tc = rd(dS, 'studyC_aggregation_error.csv');
cname = string(Tc.("case")); fname = string(T.("case"));
map = {'F2', 'C2', 'H'; 'F3', 'C3', 'D'; 'F4', 'C4', 'I'};
gap = 0;
for q = 1:3
    for x = [0.05 0.20 0.70]
        for r = [2 4 10]
            jF = find(fname == sprintf('%s_%s_m%06.4f_X%05.3f', map{q, 1}, map{q, 3}, 2*(r-1)/(r+1), x), 1);
            if isempty(jF) || T.feasible(jF) ~= 1, continue; end
            for rr = [r 1/r]
                jC = find(cname == sprintf('%s_%s_%05.3f_X_%05.3f', map{q, 2}, map{q, 3}, rr, x), 1);
                if ~isempty(jC) && ~isnan(Tc.dCCT_rel(jC))
                    gap = max(gap, abs(abs(e(jF)) - 100 * abs(Tc.dCCT_rel(jC))));
                end
            end
        end
    end
end
num('Sec. IV.C  two sweeps at shared points (pct. points)', gap, 0.20, 5e-3);
T = rd(dS, 'studyF_emf_check.csv');
s = T(abs(T.Xsep - 0.05) < 1e-9, :);
hit = s.m_E(s.limited_prefault == 1 | s.SEP_exists == 0);
num('Table XIII note: EMF spread on the limit at 0.05 pu', min(hit), 0.03, 5e-4);

%% Table VII: Thevenin equivalents per stage
T = rd(dN, 'thevenin_equivalents.csv');
ref7 = [1.000858  -6.147 0.023189 0.107511
        0.304477  -4.076 0.015312 0.054343
        0.999135  -4.935 0.024283 0.116390
        0.975805 -11.237 0.057279 0.149635
        0.228281 -15.215 0.039241 0.093678
        0.968533 -11.898 0.058730 0.151580
        1.024299 -11.787 0.049425 0.400366
        0.044739  -7.832 0.009949 0.312824
        1.005108 -13.555 0.062332 0.415100];
num('Table VII  |Vth| (pu)',       T.Vth_pu,  ref7(:, 1), 5e-7);
num('Table VII  angle of Vth (deg)', T.Vth_deg, ref7(:, 2), 5e-4);
num('Table VII  Rth (pu)',         T.Rth_pu,  ref7(:, 3), 5e-7);
num('Table VII  Xth (pu)',         T.Xth_pu,  ref7(:, 4), 5e-7);

%% Table VIII: reduction against the full network
% round-off level numbers; they need not agree digit for digit on every
% machine, so only their size is checked
T = rd(dN, 'studyA_validation_residuals.csv');
st = {'pre', 'fault', 'post'};
res = nan(3, 2);
for i = 1:3
    k = strcmp(T.stage, st{i});
    res(i, :) = [max(T.maxAbsV_pu(k)), max(T.maxAbsPe_pu(k))];
end
fprintf('%-58s  %s\n', 'Table VIII max|dV|  pre/fault/post (pu)', sprintf('%.2e  ', res(:, 1)));
fprintf('%-58s  %s\n', '           max|dPe| pre/fault/post (pu)', sprintf('%.2e  ', res(:, 2)));
istrue('Table VIII residuals at round-off level (< 1e-14 pu)', all(res(:) < 1e-14));

%% Tables IX and X: equilibria, margin, barrier and CCT per cluster
% the tables print four decimals, most of them cut rather than rounded,
% hence 1e-4
T = rd(dN, 'ieee14_results.csv');
num('Table IX   P_m (pu)',         T.Pm_pu,              [0.4; 0.1; 0.1], 5e-7);
num('Table IX   delta_0 (rad)',    T.delta0,             [0.0757; 0.0154; 0.0383], 1e-4);
num('Table IX   delta_SEP (rad)',  T.sep_post,           [0.0789; 0.0140; 0.0366], 1e-4);
num('Table IX   delta_UEP (rad)',  T.uep_post_unwrapped, [2.3818; -3.5672; 2.9165], 1e-4);
num('Table IX   Delta (rad)',      T.Delta_rad,          [2.3029; 3.5812; 2.8799], 1e-4);
num('Table IX   Delta_ref (rad)',  T.DeltaRef_rad,       [3.1885; 2.7436; 3.3002], 1e-4);
num('Table IX   M_CL',             T.MhCL,               [0.7222; 1.3052; 0.8726], 1e-4);
num('Table X    Delta_CL (rad)',   T.alphaCL_on,         [0.1946; 0.3569; 0.6047], 1e-4);
num('Table X    M_I',              T.MI,                 [0.0845; 0.0996; 0.2099], 1e-4);
num('Table X    V_cr (pu)',        T.Vcr_pu,             [1.0473; 0.9605; 1.0663], 1e-4);
num('Table X    E_crit',           T.Ecrit,              [1.2215; 2.7973; 1.9066], 1e-4);
num('Table X    CCT, buses 2 and 8 (s)', T.CCT_s([1 3]), [0.233869; 0.453235], 5e-7);
istrue('Table X    bus 3 stable up to the 10 s horizon', T.CCT_is_gt(2) == 1);
num('Sec. IV.B  system margin',    min(T.MhCL),  0.7222, 1e-4);
num('Sec. IV.B  system CCT (s)',   min(T.CCT_s), 0.233869, 5e-7);

%% Section IV.C: 22 perturbations of the external network
T = rd(dN, 'studyB_robustness_table.csv');
P = T(2:end, :);
num('Sec. IV.C  number of cases', height(P), 22, 0);
num('Sec. IV.C  margin range',  [min(P.Msys_CL); max(P.Msys_CL)], [0.682; 0.728], 5e-4);
num('Sec. IV.C  nominal margin', T.Msys_CL(1), 0.722, 5e-4);
[dev, i] = max(abs(P.Msys_relDev));
num('Sec. IV.C  largest margin deviation (%)', 100 * dev, 5.62, 5e-3);
istrue('Sec. IV.C  ... at the pre-fault loss of line 1-2', strcmp(P.("case"){i}, 'B2_pretrip_1_2'));
istrue('Sec. IV.C  post-fault equilibrium in every case', all(P.post_eq_exists == 1));
istrue('Sec. IV.C  verdict at 100 ms unchanged in every case', ...
    all(P.stable_to_breaker == T.stable_to_breaker(1)));
num('Sec. IV.C  smallest system CCT (s)', min(P.CCT_sys_s), 0.169, 5e-4);
num('Sec. IV.C  Vth shift, line 1-2: post, pre (pu)', ...
    [P.maxVthShift_post_pu(i); P.maxVthShift_pre_pu(i)], [0.46; 0.43], 5e-3);
istrue('Sec. IV.C  ... the largest shift of the set', ...
    P.maxVthShift_post_pu(i) == max(P.maxVthShift_post_pu));

fprintf('\n%d checks, %d failed\n', ncheck, nfail);
if nargout > 0
    varargout{1} = nfail;
end

    function num(name, val, ref, tol)
        % a missing or NaN value fails
        if numel(val) ~= numel(ref)
            d = NaN;
        else
            d = max(abs(val(:) - ref(:)), [], 'includenan');
        end
        ok = d <= tol + 1e-12;
        report(name, sprintf('%12.3g %10.1g', d, tol), ok);
    end

    function istrue(name, cond)
        report(name, sprintf('%12s %10s', '', ''), ~isempty(cond) && all(cond(:)));
    end

    function report(name, mid, ok)
        ncheck = ncheck + 1;
        if ok
            tag = 'ok';
        else
            tag = 'FAILED';
            nfail = nfail + 1;
        end
        fprintf('%-58s %s  %s\n', name, mid, tag);
    end
end
