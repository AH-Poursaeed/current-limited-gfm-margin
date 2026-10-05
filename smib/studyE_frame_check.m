function studyE_frame_check()
% CCT of cases 1-4 with the priority split taken on the Thevenin-voltage
% axis, as in (A.1), and on the converter's own angle (Table XIV). Nothing
% else differs between the two runs. Done at R/X = 0, which is what the
% table prints, and at R/X = 0.08. Search window 20 s, tolerance 1e-4 s.
% Writes studyE_frame_check.csv.

fprintf('\nPriority split on the Thevenin-voltage axis and on the converter axis\n');

dstep     = 2*pi/1999;
deltaGrid = (-pi + (0:3999)'*dstep);
optsCCT   = struct('t_hi',20.0,'tol',1e-4,'f0',60,'delta_h',deltaGrid);
faultModel = struct('kV',0.2,'kX',1.0);

modes  = {'circular','d-priority','q-priority'};
frames = {'thevenin','converter'};
rhos   = [0 0.08];

rows = {};
for ir = 1:numel(rhos)
    rho = rhos(ir);
    fprintf('\n-- R/X = %.2f --\n', rho);
    fprintf('  %-4s %-11s | %-9s %-9s %-9s | %s\n', ...
        'case','frame','CCT_circ','CCT_d','CCT_q','ordering');
    fprintf('  %s\n', repmat('-', 1, 74));
    for h = 1:4
        for jf = 1:numel(frames)
            p0 = ias_cluster_params_ccpac(h);
            if rho == 0
                p0.Rnet = 0;
            else
                % keep |Z| at its R = 0 value
                X0 = p0.Xnet;
                p0.Xnet = X0/sqrt(1+rho^2);
                p0.Rnet = rho*p0.Xnet;
            end
            p0.frame = frames{jf};

            CCTv = nan(1,3);
            cen  = false(1,3); fea = true(1,3);
            for m = 1:3
                c = studyE_compute_cluster_margin(p0, deltaGrid, modes{m});
                o = studyE_compute_cct_smib(c, faultModel, optsCCT, modes{m});
                CCTv(m) = o.CCT;
                cen(m) = o.censored; fea(m) = o.feasible;
                rows(end+1,:) = {rho, h, frames{jf}, modes{m}, o.CCT, ...
                    double(o.feasible), double(o.censored), c.M_CL, o.Vcrit};
            end

            if ~all(fea)
                ord = 'infeasible';
            elseif any(cen)
                ord = 'stable to the end of the window';
            else
                [~, i] = sort(CCTv, 'descend');
                ord = strjoin(modes(i), ' > ');
            end
            fprintf('  %-4d %-11s | %9.5f %9.5f %9.5f | %s\n', ...
                h, frames{jf}, CCTv(1), CCTv(2), CCTv(3), ord);
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'rho','case_id','frame','mode','CCT_s','feasible','censored','MhCL','Vcrit'});
writetable(T, 'studyE_frame_check.csv');

fprintf('\nd-priority against q-priority:\n');
for ir = 1:numel(rhos)
    rho = rhos(ir);
    for jf = 1:numel(frames)
        nHold = 0; nInv = 0; nSkip = 0;
        for h = 1:4
            s = T.rho==rho & T.case_id==h & strcmp(T.frame,frames{jf});
            if any(T.feasible(s)==0) || any(T.censored(s)==1), nSkip = nSkip+1; continue; end
            cd = T.CCT_s(s & strcmp(T.mode,'d-priority'));
            cq = T.CCT_s(s & strcmp(T.mode,'q-priority'));
            if cd > cq, nHold = nHold + 1; else, nInv = nInv + 1; end
        end
        fprintf('  R/X = %.2f, %-10s frame : d > q in %d of 4 cases, q >= d in %d (not scored %d)\n', ...
            rho, frames{jf}, nHold, nInv, nSkip);
    end
end
end
