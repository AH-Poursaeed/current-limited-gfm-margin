function studyB_run()
% Sensitivity of the margin and the CCT to what is frozen inside the stage
% Thevenin equivalents (Section IV.C). The clusters keep their settings;
% 22 cases change the rest of the system:
%   B1  operating point: generator and GFL power, voltage magnitude and
%       angle at the non-cluster buses, and the two joint extremes
%   B2  one line out before the fault, with the power flow solved again
%   B3  generator reactance and all branch reactances
% Recorded per case: system margin, system CCT, whether the CCT stays above
% a 100 ms clearing time, and the shift of the Thevenin voltages.
% Needs STEP1_state.mat. Writes studyB_robustness_table.csv.

load('STEP1_state.mat','S');
S0 = S;

opts = struct('tc_max',10.0,'rtol',1e-7,'atol',1e-9,'tol_tc',1e-6, ...
              'max_it',70,'t_clear_threshold',0.100);

t_total = tic;

fprintf('\nMargin and CCT against changes in the frozen external network\n');
fprintf('clearing-time reference %.0f ms | tc_max = %.1f s | tol_tc = %.0e\n', ...
    1000*opts.t_clear_threshold, opts.tc_max, opts.tol_tc);

nom = studyB_pipeline(S0, opts);
fprintf('\nnominal:\n');
print_case_result('Nominal', nom, nom);

cases = define_cases(S0);
N = numel(cases);
fprintf('\n%d cases:\n', N);

results = repmat(nom, 1, N+1);
results(1) = nom;
labels       = strings(N+1,1);  labels(1)       = "Nominal";
families     = strings(N+1,1);  families(1)     = "B0";
descriptions = strings(N+1,1);  descriptions(1) = "baseline";
bands        = strings(N+1,1);  bands(1)        = "-";

for k = 1:N
    c = cases(k);
    try
        S_pert = c.apply(S0);
    catch ME
        warning('Case %s failed in apply(): %s', c.label, ME.message);
        results(k+1) = blank_result(nom);
        labels(k+1) = c.label; families(k+1) = c.family;
        descriptions(k+1) = c.description; bands(k+1) = c.band;
        continue;
    end
    try
        res = studyB_pipeline(S_pert, opts);
    catch ME
        warning('Case %s failed in pipeline: %s', c.label, ME.message);
        res = blank_result(nom);
    end
    results(k+1)      = res;
    labels(k+1)       = c.label;
    families(k+1)     = c.family;
    descriptions(k+1) = c.description;
    bands(k+1)        = c.band;
    print_case_result(c.label, res, nom);
end

gfm = S0.gfm_buses(:)';
n_buses = numel(gfm);

Msys_vec        = arrayfun(@(r) r.Msys,        results).';
CCT_sys_vec     = arrayfun(@(r) r.CCT_sys,     results).';
verdict_breaker = arrayfun(@(r) r.verdict_stable_to_breaker, results).';
verdict_eq      = arrayfun(@(r) r.verdict_post_eq_exists,    results).';

maxVthShiftPre   = nan(N+1,1);
maxVthShiftFault = nan(N+1,1);
maxVthShiftPost  = nan(N+1,1);
for k = 1:(N+1)
    dpre = zeros(1,n_buses); dfault = dpre; dpost = dpre;
    for bi = 1:n_buses
        b = gfm(bi);
        dpre(bi)   = abs(results(k).Vth_pre_map(b)   - nom.Vth_pre_map(b));
        dfault(bi) = abs(results(k).Vth_fault_map(b) - nom.Vth_fault_map(b));
        dpost(bi)  = abs(results(k).Vth_post_map(b)  - nom.Vth_post_map(b));
    end
    maxVthShiftPre(k)   = max(dpre);
    maxVthShiftFault(k) = max(dfault);
    maxVthShiftPost(k)  = max(dpost);
end

Msys_dev = (Msys_vec    - nom.Msys)    / max(1e-12, abs(nom.Msys));
CCT_dev  = (CCT_sys_vec - nom.CCT_sys) / max(1e-12, abs(nom.CCT_sys));

T = table(labels, families, descriptions, bands, ...
          Msys_vec, Msys_dev, CCT_sys_vec, CCT_dev, ...
          verdict_breaker, verdict_eq, ...
          maxVthShiftPre, maxVthShiftFault, maxVthShiftPost, ...
    'VariableNames', {'case','family','description','band', ...
                      'Msys_CL','Msys_relDev','CCT_sys_s','CCT_relDev', ...
                      'stable_to_breaker','post_eq_exists', ...
                      'maxVthShift_pre_pu','maxVthShift_fault_pu','maxVthShift_post_pu'});

writetable(T, 'studyB_robustness_table.csv');

fprintf('\nnominal: M_sys = %.4f | CCT_sys = %.4f s | above clearing reference = %d\n', ...
    nom.Msys, nom.CCT_sys, nom.verdict_stable_to_breaker);
fams = unique(families,'stable');
for ii = 1:numel(fams)
    f = fams(ii);
    if f == "B0", continue; end
    idx = strcmp(families, f);
    if ~any(idx), continue; end
    Mv = Msys_vec(idx); Cv = CCT_sys_vec(idx);
    vd = verdict_breaker(idx); ve = verdict_eq(idx);
    fprintf('%s (%d cases):\n', f, sum(idx));
    fprintf('  M_sys range : [%.4f, %.4f]  (max |dev| = %.2f%%)\n', ...
        min(Mv), max(Mv), 100*max(abs(Mv - nom.Msys))/abs(nom.Msys));
    fprintf('  CCT_sys range  : [%.4f, %.4f] s (max |dev| = %.2f%%)\n', ...
        min(Cv), max(Cv), 100*max(abs(Cv - nom.CCT_sys))/abs(nom.CCT_sys));
    nFlipBreaker = sum(vd ~= nom.verdict_stable_to_breaker);
    nFlipEq      = sum(ve ~= nom.verdict_post_eq_exists);
    fprintf('  verdict at the clearing reference changed : %d / %d\n', nFlipBreaker, sum(idx));
    fprintf('  post-fault equilibrium lost               : %d / %d\n', nFlipEq, sum(idx));
end

all_pert = (1:(N+1)) > 1;
Mv = Msys_vec(all_pert); Cv = CCT_sys_vec(all_pert);
vd = verdict_breaker(all_pert);
fprintf('\nall %d cases:\n', sum(all_pert));
fprintf('  M_sys : nominal=%.4f, range=[%.4f, %.4f] (max |dev| = %.2f%%)\n', ...
    nom.Msys, min(Mv), max(Mv), 100*max(abs(Mv-nom.Msys))/abs(nom.Msys));
fprintf('  CCT_sys  : nominal=%.4f s, range=[%.4f, %.4f] s (max |dev| = %.2f%%)\n', ...
    nom.CCT_sys, min(Cv), max(Cv), 100*max(abs(Cv-nom.CCT_sys))/abs(nom.CCT_sys));
fprintf('  verdict at the clearing reference: %d changes (nominal %d, CCT below it in %d cases)\n', ...
    sum(vd ~= nom.verdict_stable_to_breaker), nom.verdict_stable_to_breaker, sum(~vd));
fprintf('\nelapsed %.1f s\n', toc(t_total));
end

function cases = define_cases(S0)
% each case: label, family, description, band, and a function that returns
% the perturbed state
cases = struct([]);

% B1: quantities frozen at the operating point, no new power flow

cases = [cases struct( ...
    'label','B1a_S_slack_-5pct','family','B1', ...
    'description','Slack SG complex power x 0.95', ...
    'band','S_slack +/-5%', ...
    'apply',@(S) apply_S_gen_scale(S, S.sg_bus, 0.95))];
cases = [cases struct( ...
    'label','B1b_S_slack_+5pct','family','B1', ...
    'description','Slack SG complex power x 1.05', ...
    'band','S_slack +/-5%', ...
    'apply',@(S) apply_S_gen_scale(S, S.sg_bus, 1.05))];
cases = [cases struct( ...
    'label','B1c_S_GFL_-10pct','family','B1', ...
    'description','GFL complex power x 0.90', ...
    'band','S_GFL +/-10%', ...
    'apply',@(S) apply_S_gen_scale(S, S.gfl_buses(1), 0.90))];
cases = [cases struct( ...
    'label','B1d_S_GFL_+10pct','family','B1', ...
    'description','GFL complex power x 1.10', ...
    'band','S_GFL +/-10%', ...
    'apply',@(S) apply_S_gen_scale(S, S.gfl_buses(1), 1.10))];
cases = [cases struct( ...
    'label','B1e_Vop_mag_-2pct','family','B1', ...
    'description','|Vop| at all non-cluster buses x 0.98', ...
    'band','|V| +/-2%', ...
    'apply',@(S) apply_V_mag_scale_nonCluster(S, 0.98))];
cases = [cases struct( ...
    'label','B1f_Vop_mag_+2pct','family','B1', ...
    'description','|Vop| at all non-cluster buses x 1.02', ...
    'band','|V| +/-2%', ...
    'apply',@(S) apply_V_mag_scale_nonCluster(S, 1.02))];
cases = [cases struct( ...
    'label','B1g_Vop_ang_-1deg','family','B1', ...
    'description','ang(Vop) at all non-cluster buses - 1 deg', ...
    'band','ang(V) +/-1 deg', ...
    'apply',@(S) apply_V_angle_shift_nonCluster(S, -deg2rad(1)))];
cases = [cases struct( ...
    'label','B1h_Vop_ang_+1deg','family','B1', ...
    'description','ang(Vop) at all non-cluster buses + 1 deg', ...
    'band','ang(V) +/-1 deg', ...
    'apply',@(S) apply_V_angle_shift_nonCluster(S, +deg2rad(1)))];
cases = [cases struct( ...
    'label','B1i_joint_neg','family','B1', ...
    'description','Joint -: S_slack x0.95, S_GFL x0.90, |V| x0.98, ang(V) -1 deg', ...
    'band','joint worst-case', ...
    'apply',@(S) apply_joint(S, -1))];
cases = [cases struct( ...
    'label','B1j_joint_pos','family','B1', ...
    'description','Joint +: S_slack x1.05, S_GFL x1.10, |V| x1.02, ang(V) +1 deg', ...
    'band','joint worst-case', ...
    'apply',@(S) apply_joint(S, +1))];

% B2: one of eight lines out before the fault, near the clusters and
% far from them. Line 7-8 (it would island bus 8) and line 2-4 (tripped
% after the fault) are not among them.
B2_lines = [ ...
    1 2;
    2 3;
    2 5;
    3 4;
    4 5;
    6 11;
    10 11;
    13 14];
for ii = 1:size(B2_lines,1)
    a = B2_lines(ii,1); b = B2_lines(ii,2);
    lbl = sprintf('B2_pretrip_%d_%d', a, b);
    cases = [cases struct( ...
        'label',lbl,'family','B2', ...
        'description',sprintf('Pre-fault N-1: open line %d-%d, PF re-run', a, b), ...
        'band','N-1 line outage (PF re-run)', ...
        'apply',@(S) apply_pretrip_with_PF(S, [a b]))];
end

% B3: network parameters
cases = [cases struct( ...
    'label','B3a_Xdp_-25pct','family','B3', ...
    'description','Xdprime_sg x 0.75', ...
    'band','Xdprime +/-25%', ...
    'apply',@(S) apply_Xdprime_scale(S, 0.75))];
cases = [cases struct( ...
    'label','B3b_Xdp_+25pct','family','B3', ...
    'description','Xdprime_sg x 1.25', ...
    'band','Xdprime +/-25%', ...
    'apply',@(S) apply_Xdprime_scale(S, 1.25))];
cases = [cases struct( ...
    'label','B3c_branchX_-10pct','family','B3', ...
    'description','All branch reactances x 0.90 (PF re-run)', ...
    'band','branch X +/-10%', ...
    'apply',@(S) apply_branchX_scale_with_PF(S, 0.90))];
cases = [cases struct( ...
    'label','B3d_branchX_+10pct','family','B3', ...
    'description','All branch reactances x 1.10 (PF re-run)', ...
    'band','branch X +/-10%', ...
    'apply',@(S) apply_branchX_scale_with_PF(S, 1.10))];
end

function S = apply_S_gen_scale(S, gen_bus, scale)
    S.Sgen_bus(gen_bus) = S.Sgen_bus(gen_bus) * scale;
end

function S = apply_V_mag_scale_nonCluster(S, scale)
    nonC = setdiff(1:S.nb, S.gfm_buses);
    for bb = nonC
        S.Vop(bb) = abs(S.Vop(bb)) * scale * exp(1j*angle(S.Vop(bb)));
    end
end

function S = apply_V_angle_shift_nonCluster(S, dtheta)
    nonC = setdiff(1:S.nb, S.gfm_buses);
    for bb = nonC
        S.Vop(bb) = abs(S.Vop(bb)) * exp(1j*(angle(S.Vop(bb)) + dtheta));
    end
end

function S = apply_joint(S, sgn)
    if sgn < 0
        S = apply_S_gen_scale(S, S.sg_bus, 0.95);
        S = apply_S_gen_scale(S, S.gfl_buses(1), 0.90);
        S = apply_V_mag_scale_nonCluster(S, 0.98);
        S = apply_V_angle_shift_nonCluster(S, -deg2rad(1));
    else
        S = apply_S_gen_scale(S, S.sg_bus, 1.05);
        S = apply_S_gen_scale(S, S.gfl_buses(1), 1.10);
        S = apply_V_mag_scale_nonCluster(S, 1.02);
        S = apply_V_angle_shift_nonCluster(S, +deg2rad(1));
    end
end

function S = apply_pretrip_with_PF(S, pretrip_pair)
    define_constants;
    % S.mpc0 is case14 as loaded, before the dispatch of step1
    mpc = S.mpc0;
    mpc = apply_stressed_dispatch(mpc, S);
    br  = mpc.branch;
    a = pretrip_pair(1); b = pretrip_pair(2);
    idx = find((br(:,F_BUS)==a & br(:,T_BUS)==b) | (br(:,F_BUS)==b & br(:,T_BUS)==a), 1);
    if isempty(idx)
        error('apply_pretrip_with_PF: line %d-%d not found', a, b);
    end
    br(idx, BR_STATUS) = 0;
    mpc.branch = br;

    mpopt = mpoption('verbose',0,'out.all',0);
    [res, success] = runpf(mpc, mpopt);
    if success ~= 1
        error('apply_pretrip_with_PF: PF did not converge for pre-trip %d-%d', a, b);
    end

    S.mpc0     = mpc;
    S.Vop      = res.bus(:,VM) .* exp(1j*deg2rad(res.bus(:,VA)));
    S.Sgen_bus = zeros(S.nb,1);
    for g = 1:size(res.gen,1)
        bb = res.gen(g, GEN_BUS);
        S.Sgen_bus(bb) = S.Sgen_bus(bb) + (res.gen(g, PG) + 1j*res.gen(g, QG))/S.baseMVA;
    end
    S.Sload    = res.bus(:,PD)/S.baseMVA + 1j*res.bus(:,QD)/S.baseMVA;
    S.Yload    = conj(S.Sload) ./ max(1e-12, abs(S.Vop).^2);
    S.Yload_mat= spdiags(S.Yload, 0, S.nb, S.nb);
    [S.Ybus_pre, ~, ~] = makeYbus(S.baseMVA, res.bus, res.branch);
end

function S = apply_Xdprime_scale(S, scale)
    S.Xdprime_sg = S.Xdprime_sg * scale;
    S.Zs_sg = 1j*S.Xdprime_sg;
    S.Ys_sg = 1/S.Zs_sg;
end

function S = apply_branchX_scale_with_PF(S, scale)
    define_constants;
    mpc = S.mpc0;
    mpc = apply_stressed_dispatch(mpc, S);
    mpc.branch(:, BR_X) = mpc.branch(:, BR_X) * scale;
    mpopt = mpoption('verbose',0,'out.all',0);
    [res, success] = runpf(mpc, mpopt);
    if success ~= 1
        error('apply_branchX_scale_with_PF: PF did not converge for scale=%.3f', scale);
    end
    S.mpc0 = mpc;
    S.Vop  = res.bus(:,VM) .* exp(1j*deg2rad(res.bus(:,VA)));
    S.Sgen_bus = zeros(S.nb,1);
    for g = 1:size(res.gen,1)
        bb = res.gen(g, GEN_BUS);
        S.Sgen_bus(bb) = S.Sgen_bus(bb) + (res.gen(g, PG) + 1j*res.gen(g, QG))/S.baseMVA;
    end
    S.Sload    = res.bus(:,PD)/S.baseMVA + 1j*res.bus(:,QD)/S.baseMVA;
    S.Yload    = conj(S.Sload) ./ max(1e-12, abs(S.Vop).^2);
    S.Yload_mat= spdiags(S.Yload, 0, S.nb, S.nb);
    [S.Ybus_pre, ~, ~] = makeYbus(S.baseMVA, res.bus, res.branch);
end

function mpc = apply_stressed_dispatch(mpc, S)
% same dispatch as in step1
    define_constants;
    if ~isfield(S,'Pg_set_MW') || isempty(S.Pg_set_MW)
        return;
    end
    for k = 1:size(mpc.gen,1)
        bb = mpc.gen(k, GEN_BUS);
        if isKey(S.Pg_set_MW, bb)
            Pg_new = S.Pg_set_MW(bb);
            Pmin   = mpc.gen(k, PMIN);
            Pmax   = mpc.gen(k, PMAX);
            mpc.gen(k, PG) = min(max(Pg_new, Pmin), Pmax);
        end
    end
end

function print_case_result(label, res, nom)
    if isempty(res.gfm_buses)
        fprintf('  %-30s  [FAILED]\n', label);
        return;
    end
    fprintf('  %-30s  M_sys=%.4f (%+6.2f%%)  CCT_sys=%.4f s (%+6.2f%%)  brkr=%d  eq=%d\n', ...
        label, res.Msys, ...
        100*(res.Msys - nom.Msys)/max(1e-12,abs(nom.Msys)), ...
        res.CCT_sys, ...
        100*(res.CCT_sys - nom.CCT_sys)/max(1e-12,abs(nom.CCT_sys)), ...
        res.verdict_stable_to_breaker, res.verdict_post_eq_exists);
end

function r = blank_result(nom)
    r = nom;
    r.Msys = NaN; r.CCT_sys = NaN;
    r.MhCL = nan(size(nom.MhCL));
    r.CCT_h = nan(size(nom.CCT_h));
    r.verdict_stable_to_breaker = false;
    r.verdict_post_eq_exists = false;
end
