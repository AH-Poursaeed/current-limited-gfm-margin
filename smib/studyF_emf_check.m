function studyF_emf_check()
% Pre-fault currents of the two converters for a spread m_E in the internal
% EMF magnitude, at three separations. An EMF difference drives a
% circulating current of about m_E * E / X_sep between the units; with
% Imax = 0.6 pu and E = 1 pu it reaches the limit near m_E = 0.6 * X_sep.
% This is the basis of the EMF remark under Table XIII. No dynamics.
% Writes studyF_emf_check.csv.

X_list = [0.05 0.20 0.70];
m_list = [0.001 0.002 0.005 0.01 0.02 0.03 0.05 0.075 0.10 0.125 0.15 0.20 0.25 0.30];
SUM_E  = 2.0;

fprintf('\nEMF spread and pre-fault circulating current\n');
rows = struct([]);
for x = X_list
    fprintf('\nX_sep = %.2f pu, limit expected near m_E = %.3f\n', ...
        x, 0.6*x/(SUM_E/2));
    fprintf('%8s %10s %10s %10s %10s %12s %s\n', ...
        'm_E','|I1| pu','|I2| pu','I_circ est','I_max','on limit','op. point');
    for m = m_list
        E1 = (SUM_E/2)*(1+m/2);  E2 = (SUM_E/2)*(1-m/2);
        P = studyC_make_2GFM('Xa',x/2,'Xb',x/2,'E1',E1,'E2',E2);
        [d1, d2] = studyC_initial_sep(P);
        sep_ok = ~isnan(d1);
        if sep_ok
            [~,~,I1,I2,~,ls] = studyC_per_converter_pe(d1, d2, P, 'pre');
            a1 = abs(I1); a2 = abs(I2); lim = any(ls);
        else
            a1 = NaN; a2 = NaN; lim = NaN;
        end
        icirc = m*(SUM_E/2)/x;
        if isnan(lim), limstr = '-'; else, limstr = char(string(logical(lim))); end
        fprintf('%8.4f %10.4f %10.4f %10.4f %10.2f %12s %s\n', ...
            m, a1, a2, icirc, P.Imax1, limstr, char(string(sep_ok)));
        r.Xsep = x; r.m_E = m; r.E1 = E1; r.E2 = E2;
        r.absI1 = a1; r.absI2 = a2; r.I_circ_est = icirc;
        r.Imax = P.Imax1;
        if isnan(lim), r.limited_prefault = NaN; else, r.limited_prefault = double(lim); end
        r.SEP_exists = double(sep_ok);
        if isempty(rows), rows = r; else, rows(end+1) = r; end
    end
end

T = struct2table(rows,'AsArray',true);
writetable(T,'studyF_emf_check.csv');

fprintf('\n');
for x = X_list
    s = T(T.Xsep == x,:);
    pred = 0.6*x/(SUM_E/2);
    hit = s.m_E(s.limited_prefault == 1 | s.SEP_exists == 0 | isnan(s.limited_prefault));
    if isempty(hit)
        fprintf('X_sep = %.2f : not on the limit up to m_E = %.3f (expected %.3f)\n', ...
            x, max(s.m_E), pred);
    else
        fprintf('X_sep = %.2f : on the limit from m_E = %.3f (expected %.3f)\n', ...
            x, min(hit), pred);
    end
end
end
