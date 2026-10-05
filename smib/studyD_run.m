function studyD_run()
% CCT and margin of cases 1-4 under circular, d-priority and q-priority
% current limiting (Table III), and the same for case 1 over a range of
% fault depths kV. Writes studyD_main_sweep.csv and studyD_kV_sweep.csv.

t_total = tic;
Npts = 2000;
deltaGrid = linspace(-pi, +pi, Npts).';
optsCCT = struct('t_hi',2.0,'tol',1e-3,'f0',60,'delta_h',deltaGrid);
faultModel = struct('kV',0.2,'kX',1.0);

modes = {'circular','d-priority','q-priority'};
numCases = 4;

fprintf('\nPriority current limiting, single-converter cases 1-4\n');

% with the circular limiter the studyD_* routines have to return the
% values of Table II; checked before anything else is computed
target_CCT = [0.40625 0.40625 0.44434 0.53516];
target_M_case1 = 0.79404;

CCT_chk = nan(numCases,1);
M_chk   = nan(numCases,1);
for h = 1:numCases
    p = ias_cluster_params_ccpac(h);
    c = studyD_compute_cluster_margin(p, deltaGrid, 'circular');
    o = studyD_compute_cct_smib(c, faultModel, optsCCT, 'circular');
    CCT_chk(h) = o.CCT;
    M_chk(h)   = c.M_CL;
end
dCCT_chk  = CCT_chk - target_CCT(:);
dM1_chk   = M_chk(1) - target_M_case1;

fprintf('\ncircular limiter against Table II:\n');
for h = 1:numCases
    fprintf('  Case %d : CCT = %.6f s (target %.5f, dev %.2e) | M = %.6f\n', ...
        h, CCT_chk(h), target_CCT(h), abs(dCCT_chk(h)), M_chk(h));
end
fprintf('  Case 1 margin, deviation from 0.79404 : %.2e\n', abs(dM1_chk));
chk_ok = max(abs(dCCT_chk)) < 5e-4 && abs(dM1_chk) < 5e-4;
if ~chk_ok
    error('studyD_run:check', 'circular limiter does not reproduce Table II');
end

CCT_mat = nan(numCases, numel(modes));
M_mat   = nan(numCases, numel(modes));
fprintf('\nfour cases, three limiter modes, kV = %.2f:\n', faultModel.kV);
fprintf('  %-6s %-12s | %-9s | %-9s\n', 'Case', 'Mode', 'CCT (s)', 'M_CL');
fprintf('  %s\n', repmat('-', 1, 50));
for h = 1:numCases
    p = ias_cluster_params_ccpac(h);
    for m = 1:numel(modes)
        mode = modes{m};
        c = studyD_compute_cluster_margin(p, deltaGrid, mode);
        o = studyD_compute_cct_smib(c, faultModel, optsCCT, mode);
        CCT_mat(h,m) = o.CCT;
        M_mat(h,m)   = c.M_CL;
        fprintf('  Case %-1d %-12s | %9.6f | %9.6f\n', h, mode, o.CCT, c.M_CL);
    end
end

% differences are taken against the circular limiter of the same case
Case = repmat((1:numCases).', numel(modes), 1);
Mode = repmat(modes, numCases, 1); Mode = Mode(:);
CCT  = CCT_mat(:);
MhCL = M_mat(:);
CCT_circular = repmat(CCT_mat(:,1), numel(modes), 1);
M_circular   = repmat(M_mat(:,1),   numel(modes), 1);
dCCT_vs_circ = CCT - CCT_circular;
dCCT_rel_pct = 100 * dCCT_vs_circ ./ max(1e-12, CCT_circular);
dM_vs_circ   = MhCL - M_circular;
dM_rel_pct   = 100 * dM_vs_circ ./ max(1e-12, M_circular);
T_main = table(Case, Mode, CCT, MhCL, dCCT_vs_circ, dCCT_rel_pct, dM_vs_circ, dM_rel_pct, ...
    'VariableNames', {'Case','Mode','CCT_s','MhCL','dCCT_vs_circ_s','dCCT_rel_pct','dM_vs_circ','dM_rel_pct'});

writetable(T_main, 'studyD_main_sweep.csv');

% case 1 for different fault depths (kV = fault-on voltage / pre-fault)
kV_grid = [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90];
CCT_kV = nan(numel(kV_grid), numel(modes));
M_kV   = nan(numel(kV_grid), numel(modes));

fprintf('\ncase 1 against kV:\n');
fprintf('  %-6s | %-9s | %-9s | %-9s\n', 'kV', 'circular', 'd-prior', 'q-prior');
fprintf('  %s\n', repmat('-', 1, 48));
for k = 1:numel(kV_grid)
    fm = faultModel; fm.kV = kV_grid(k);
    p  = ias_cluster_params_ccpac(1);
    for m = 1:numel(modes)
        mode = modes{m};
        c = studyD_compute_cluster_margin(p, deltaGrid, mode);
        o = studyD_compute_cct_smib(c, fm, optsCCT, mode);
        CCT_kV(k,m) = o.CCT;
        M_kV(k,m)   = c.M_CL;
    end
    fprintf('  %5.2f  |  %7.4f  |  %7.4f  |  %7.4f\n', kV_grid(k), CCT_kV(k,1), CCT_kV(k,2), CCT_kV(k,3));
end

T_kV = table(kV_grid(:), CCT_kV(:,1), CCT_kV(:,2), CCT_kV(:,3), M_kV(:,1), M_kV(:,2), M_kV(:,3), ...
    'VariableNames', {'kV','CCT_circular','CCT_d_priority','CCT_q_priority','M_circular','M_d_priority','M_q_priority'});
writetable(T_kV, 'studyD_kV_sweep.csv');

% mean shift over the four cases, relative to the mean circular CCT
gap_d = mean(CCT_mat(:,2) - CCT_mat(:,1));
gap_q = mean(CCT_mat(:,3) - CCT_mat(:,1));
fprintf('\nmean CCT shift against circular:\n');
fprintf('  d-priority : %+.4f s (%+.2f %%)\n', gap_d, 100*gap_d/mean(CCT_mat(:,1)));
fprintf('  q-priority : %+.4f s (%+.2f %%)\n', gap_q, 100*gap_q/mean(CCT_mat(:,1)));
fprintf('elapsed %.1f s\n', toc(t_total));
end
