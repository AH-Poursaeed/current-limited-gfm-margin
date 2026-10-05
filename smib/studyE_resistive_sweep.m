function studyE_resistive_sweep()
% Does the d-priority > circular > q-priority ordering of Table III hold
% when the coupling has resistance? (Section IV.A.) R/X is swept from 0 to
% 2 for the four cases and the three limiter modes, in two ways:
%   A_constZ   |Z| kept at the R = 0 value, so only the X/R ratio changes
%   B_constX   X kept, R added on top
% The sweep uses a 20 s search window and a 1e-4 s tolerance. A point is
% not scored if a mode has no SEP/UEP pair, reaches the end of the window,
% or if the limiter never acts on the post-fault curve.
% Writes studyE_resistive_sweep.csv (rewritten after every R/X value).

t_total = tic;

% two periods of the angle, with the step of the 2000-point grid used
% for Table III. With resistance the UEP can move beyond pi.
dstep     = 2*pi/1999;
deltaGrid = (-pi + (0:3999)'*dstep);
faultModel = struct('kV',0.2,'kX',1.0);

optsRef  = struct('t_hi',2.0, 'tol',1e-3, 'f0',60,'delta_h',deltaGrid);
optsSweep = struct('t_hi',20.0,'tol',1e-4, 'f0',60,'delta_h',deltaGrid);

modes    = {'circular','d-priority','q-priority'};
numCases = 4;

rho_grid = [0 0.02 0.05 0.08 0.10 0.15 0.20 0.25 0.30 0.40 ...
            0.50 0.60 0.70 0.80 1.00 1.25 1.50 2.00];

params  = {'A_constZ','B_constX'};
csvFile = 'studyE_resistive_sweep.csv';

fprintf('\nPriority-limiter ordering against network resistance\n');
fprintf('%d parameterizations x %d R/X values x %d cases x %d modes = %d CCTs\n', ...
    numel(params), numel(rho_grid), numCases, numel(modes), ...
    numel(params)*numel(rho_grid)*numCases*numel(modes));
fprintf('check : t_hi = %.1f s, tol = %.0e (as in Table III)\n', optsRef.t_hi, optsRef.tol);
fprintf('sweep : t_hi = %.1f s, tol = %.0e\n', optsSweep.t_hi, optsSweep.tol);

% R = 0 has to give Table III exactly. The bisection starts from [0, 2] s
% and stops at 1e-3 s, so every CCT is a multiple of 1/1024 s; the table
% prints these rounded to five decimals.
target_CCT = [ 416 416 455 548
               578 649 660 594
               299 221 281 467 ]/1024;
target_M_case1_circ = 0.79404;

fprintf('\nR = 0 against Table III:\n');
ref_dev = zeros(numel(modes), numCases);
for h = 1:numCases
    p = ias_cluster_params_ccpac(h); p.Rnet = 0;
    for m = 1:numel(modes)
        c = studyE_compute_cluster_margin(p, deltaGrid, modes{m});
        o = studyE_compute_cct_smib(c, faultModel, optsRef, modes{m});
        ref_dev(m,h) = o.CCT - target_CCT(m,h);
        if m == 1 && h == 1, ref_M1 = c.M_CL; end
    end
end
fprintf('  max |CCT - Table III| over 12 entries : %.3e s\n', max(abs(ref_dev(:))));
fprintf('  entries reproduced exactly            : %d of 12\n', nnz(ref_dev == 0));
fprintf('  case 1 circular margin                : %.6f (%.5f)\n', ref_M1, target_M_case1_circ);
if max(abs(ref_dev(:))) ~= 0 || abs(ref_M1 - target_M_case1_circ) > 5e-6
    error('studyE_resistive_sweep:check', 'R = 0 does not reproduce Table III');
end

rows = {};
fprintf('\nsweep, flags: I no SEP/UEP pair, C end of search window, N limiter inactive\n');
for ip = 1:numel(params)
    parName = params{ip};
    fprintf('\n[%s]\n', parName);
    fprintf('  %-5s %-4s | %-9s %-9s %-9s | %-6s %-6s | %s\n', ...
        'rho','case','CCT_circ','CCT_d','CCT_q','negIP0','flags','verdict');
    fprintf('  %s\n', repmat('-', 1, 92));

    for ir = 1:numel(rho_grid)
        rho = rho_grid(ir);
        for h = 1:numCases
            p0 = ias_cluster_params_ccpac(h);
            X0 = p0.Xnet;
            if rho == 0
                Xn = X0; Rn = 0;
            else
                switch parName
                    case 'A_constZ', Xn = X0/sqrt(1+rho^2); Rn = rho*Xn;
                    case 'B_constX', Xn = X0;               Rn = rho*X0;
                end
            end

            CCTv = nan(1,numel(modes));  fea = false(1,numel(modes));
            cen  = false(1,numel(modes)); fneg = nan(1,numel(modes));
            act  = false(1,numel(modes));

            for m = 1:numel(modes)
                p = p0; p.Xnet = Xn; p.Rnet = Rn;
                c = studyE_compute_cluster_margin(p, deltaGrid, modes{m});
                o = studyE_compute_cct_smib(c, faultModel, optsSweep, modes{m});
                CCTv(m) = o.CCT;  fea(m) = o.feasible;
                cen(m)  = o.censored;  fneg(m) = o.fracNegIP0;
                % if the unlimited current never exceeds Imax the three
                % modes are the same curve and the point says nothing
                I0max  = max(c.I_unlim);
                act(m) = I0max > p.Imax;

                rows(end+1,:) = {parName, rho, h, X0, Xn, Rn, hypot(Rn,Xn), ...
                    modes{m}, o.CCT, double(o.feasible), double(o.censored), ...
                    double(act(m)), I0max/p.Imax, ...
                    c.M_CL, o.Vcrit, o.Aacc, o.fracNegIP0, ...
                    o.delta_SEP_post*180/pi, o.delta_UEP_post*180/pi, ...
                    c.Delta_lim*180/pi, c.Delta_ref*180/pi};
            end

            [verdict, flags] = local_verdict(CCTv, fea, cen, act);
            fprintf('  %-5.2f %-4d | %9.5f %9.5f %9.5f | %6.3f %-6s | %s\n', ...
                rho, h, CCTv(1), CCTv(2), CCTv(3), fneg(2), flags, verdict);
        end
        writetable(local_table(rows), csvFile);
    end
end

T = local_table(rows);
writetable(T, csvFile);
% last R/X at which d > q (or d > circular) and first R/X at which it
% no longer holds, per case
fprintf('\nsummary:\n');
nWell3 = 0; nInv3 = 0; nHold3 = 0; nInactive = 0;
maxNegHold = -inf; minNegInv = inf;
for ip = 1:numel(params)
    parName = params{ip};
    fprintf('\n[%s]\n', parName);
    fprintf('  %-4s | %-34s | %-34s\n', 'case', ...
        'd vs q  (last d>q / first q>=d)', 'd vs circ (last d>c / first c>=d)');
    fprintf('  %s\n', repmat('-', 1, 80));
    for h = 1:numCases
        lastDQ = NaN; firstQD = NaN; negAtQD = NaN;
        lastDC = NaN; firstCD = NaN; negAtCD = NaN;
        for ir = 1:numel(rho_grid)
            sel = strcmp(T.param,parName) & T.rho==rho_grid(ir) & T.case_id==h;
            gc = sel & strcmp(T.mode,'circular');
            gd = sel & strcmp(T.mode,'d-priority');
            gq = sel & strcmp(T.mode,'q-priority');
            live = T.limiter_active(gd)==1 && T.limiter_active(gq)==1;
            okc = T.feasible(gc)==1 && T.censored(gc)==0;
            okd = T.feasible(gd)==1 && T.censored(gd)==0;
            okq = T.feasible(gq)==1 && T.censored(gq)==0;
            cc = T.CCT_s(gc); cd = T.CCT_s(gd); cq = T.CCT_s(gq);

            if okd && okq && ~live, nInactive = nInactive + 1; end

            if okd && okq && live
                nWell3 = nWell3 + 1;
                if cd > cq
                    nHold3 = nHold3 + 1;
                    lastDQ = rho_grid(ir);
                    maxNegHold = max(maxNegHold, T.frac_negIP0(gd));
                else
                    nInv3 = nInv3 + 1;
                    minNegInv = min(minNegInv, T.frac_negIP0(gd));
                    if isnan(firstQD)
                        firstQD = rho_grid(ir); negAtQD = T.frac_negIP0(gd);
                    end
                end
            end
            if okd && okc && live
                if cd > cc
                    lastDC = rho_grid(ir);
                elseif isnan(firstCD)
                    firstCD = rho_grid(ir); negAtCD = T.frac_negIP0(gd);
                end
            end
        end
        fprintf('  %-4d | %-6.2f / %-6.2f (negIP0 %-5.3f)      | %-6.2f / %-6.2f (negIP0 %-5.3f)\n', ...
            h, lastDQ, firstQD, negAtQD, lastDC, firstCD, negAtCD);
    end
end

nAll = height(T)/3;
fprintf('\nd-priority against q-priority, whole grid:\n');
fprintf('  grid points                                : %d\n', nAll);
fprintf('  not scored, a mode without SEP/UEP pair    : %d\n', local_countflag(T,'feasible',0));
fprintf('  not scored, a mode stable up to %.0f s       : %d\n', optsSweep.t_hi, local_countflag(T,'censored',1));
fprintf('  not scored, limiter inactive               : %d\n', nInactive);
fprintf('  scored                                     : %d\n', nWell3);
fprintf('    d > q, as in Table III                   : %d\n', nHold3);
fprintf('    q >= d                                   : %d\n', nInv3);
fprintf('  largest frac_negIP0 among d > q rows       : %.3f\n', maxNegHold);
fprintf('  smallest frac_negIP0 among inversion rows  : %.3f\n', minNegInv);

fprintf('\nordering d > circular > q for R/X <= 0.25:\n');
selLow = T.rho <= 0.25;
nLow = 0; nLowHold = 0; nLowSkip = 0;
u = unique(strcat(T.param, '|', string(T.rho), '|', string(T.case_id)));
key = strcat(T.param, '|', string(T.rho), '|', string(T.case_id));
for i = 1:numel(u)
    s = key == u(i);
    if ~any(s & selLow), continue; end
    if any(T.feasible(s)==0) || any(T.censored(s)==1) || any(T.limiter_active(s)==0)
        nLowSkip = nLowSkip + 1; continue;
    end
    nLow = nLow + 1;
    cc = T.CCT_s(s & strcmp(T.mode,'circular'));
    cd = T.CCT_s(s & strcmp(T.mode,'d-priority'));
    cq = T.CCT_s(s & strcmp(T.mode,'q-priority'));
    if cd > cc && cc > cq, nLowHold = nLowHold + 1; end
end
fprintf('  grid points scored : %d (not scored: %d)\n', nLow, nLowSkip);
fprintf('  ordering holds in  : %d of %d\n', nLowHold, nLow);
fprintf('elapsed %.1f s, results in %s\n', toc(t_total), csvFile);
end

function T = local_table(rows)
T = cell2table(rows, 'VariableNames', {'param','rho','case_id','X0','Xnet','Rnet','Znet_mag', ...
    'mode','CCT_s','feasible','censored','limiter_active','I0max_over_Imax', ...
    'MhCL','Vcrit','A_acc_faulton','frac_negIP0', ...
    'delta_SEP_post_deg','delta_UEP_post_deg','Delta_lim_deg','Delta_ref_deg'});
end

function [verdict, flags] = local_verdict(CCTv, fea, cen, act)
flags = '';
if ~all(fea), flags = [flags 'I']; end
if any(cen),  flags = [flags 'C']; end
if ~all(act), flags = [flags 'N']; end
if isempty(flags), flags = '-'; end
if ~all(fea)
    verdict = 'no SEP/UEP pair in some mode';
elseif ~all(act)
    verdict = 'limiter inactive';
elseif any(cen)
    verdict = 'stable to the end of the window';
elseif CCTv(2) > CCTv(3) && CCTv(2) > CCTv(1)
    verdict = 'd > circ > q';
elseif CCTv(2) > CCTv(3)
    verdict = 'd > q, circ >= d';
else
    verdict = 'q >= d';
end
end

function n = local_countflag(T, col, val)
k = 0;
u = unique(strcat(T.param, '|', string(T.rho), '|', string(T.case_id)));
for i = 1:numel(u)
    sel = strcat(T.param, '|', string(T.rho), '|', string(T.case_id)) == u(i);
    if any(T.(col)(sel) == val), k = k + 1; end
end
n = k;
end
