function studyF_envelope()
% Reads studyF_tolerance_sweep.csv and reports, per electrical separation
% and for error bounds of 1, 2, 5 and 10 %, the largest tested spread of
% each setting such that every tested level up to it has an error value
% within the bound. The next level tested is reported with it, because
% these are grid values and not thresholds. Table XIII is the 5 % bound at
% X_sep = 0.05 pu.
% Also printed: the sign of the error, the error from separation alone,
% and the difference from the Table XII data at the points the two grids
% share. Writes studyF_admissibility_envelope.csv.

CSV = 'studyF_tolerance_sweep.csv';
T = readtable(CSV, 'VariableNamingRule','preserve');
fprintf('\n%s: %d cases\n', CSV, height(T));

fam_all = string(T.family);
feas    = T.feasible == 1;

i0 = fam_all == "F1" & T.Xsep == 0;
if any(i0)
    floor_pct = 100*abs(T.dCCT_rel(i0));
else
    floor_pct = NaN;
end

fprintf('with an error value: %d, without: %d\n', sum(feas), sum(~feas));
fprintf('matched pair at the hub: %.4f %%\n', floor_pct);

n_over  = sum(feas & T.dCCT_rel < 0);
n_under = sum(feas & T.dCCT_rel > 0);
n_exact = sum(feas & T.dCCT_rel == 0);
fprintf('aggregate CCT longer in %d cases, shorter in %d, equal in %d\n', ...
    n_over, n_under, n_exact);
if n_under > 0
    fprintf('cases where the aggregate CCT is the shorter one:\n');
    idx = find(feas & T.dCCT_rel > 0);
    for q = 1:numel(idx)
        fprintf('    %-26s  m = %.4f  X_sep = %.3f  dCCT_rel = %+.4f %%\n', ...
            string(T.("case")(idx(q))), T.mismatch(idx(q)), T.Xsep(idx(q)), ...
            100*T.dCCT_rel(idx(q)));
    end
end
dd = 100*abs(T.dCCT_rel); dd(~feas) = -inf;
[mx, imx] = max(dd);
fprintf('largest |error|: %.3f %% in %s\n', mx, string(T.("case")(imx)));

fprintf('\nseparation alone, parameters matched (F1):\n');
T1 = sortrows(T(fam_all == "F1",:), 'Xsep');
fprintf('  %10s  %12s\n', 'X_sep (pu)', 'error (%)');
for q = 1:height(T1)
    fprintf('  %10.4f  %+12.4f\n', T1.Xsep(q), 100*T1.dCCT_rel(q));
end
bounds = [1 2 5 10];
fprintf('\n  largest tested X_sep within each bound [next level tested]:\n');
sep_env = struct([]);
for b = bounds
    [xok, eok, xnext, enext] = walk(T1.Xsep, T1.dCCT_rel, T1.feasible, b);
    fprintf('   <= %2d %% : X_sep = %.4f pu at %.3f %%   [%.4f pu at %.3f %%]\n', ...
        b, xok, eok, xnext, enext);
    s.quantity = 'X_sep  electrical separation'; s.family = 'F1';
    s.Xsep = NaN; s.bound_pct = b;
    s.max_admissible = xok; s.err_at_max = eok;
    s.next_tested = xnext; s.err_at_next = enext;
    if isempty(sep_env), sep_env = s; else, sep_env(end+1) = s; end
end

fams  = {'F2','F3','F4','F5','F7','F6','F8'};
names = {'H  inertia', 'D  damping', 'Imax  current limit', ...
         'Pref  pre-fault dispatch', 'E  internal EMF magnitude', ...
         'all four, aligned', 'all four, adverse'};

Xs = unique(T.Xsep(fam_all == "F2"));
env = struct([]);

for b = bounds
    fprintf('\nlargest admissible spread for |error| <= %d %% (error reached):\n', b);
    fprintf('  %-26s', 'setting');
    for x = Xs', fprintf('   X_sep=%.2f pu      ', x); end
    fprintf('\n');
    for f = 1:numel(fams)
        idx_f = fam_all == fams{f};
        if ~any(idx_f), continue; end
        fprintf('  %-26s', names{f});
        for x = Xs'
            sub = sortrows(T(idx_f & T.Xsep == x,:), 'mismatch');
            [mok, eok, mnext, enext] = walk(sub.mismatch, sub.dCCT_rel, sub.feasible, b);
            fprintf('   %6.3f (%5.2f %%)  ', mok, eok);
            e.quantity = names{f}; e.family = fams{f};
            e.Xsep = x; e.bound_pct = b;
            e.max_admissible = mok; e.err_at_max = eok;
            e.next_tested = mnext; e.err_at_next = enext;
            if isempty(env), env = e; else, env(end+1) = e; end
        end
        fprintf('\n');
    end
end

fprintf('\nlargest spread up to which every case has an error value (cases without):\n');
fprintf('  %-26s', 'setting');
for x = Xs', fprintf('   X_sep=%.2f', x); end
fprintf('\n');
for f = 1:numel(fams)
    idx_f = fam_all == fams{f};
    if ~any(idx_f), continue; end
    fprintf('  %-26s', names{f});
    for x = Xs'
        sub = sortrows(T(idx_f & T.Xsep == x,:), 'mismatch');
        lv = unique(sub.mismatch); mf = 0;
        for q = 1:numel(lv)
            if all(sub.feasible(sub.mismatch <= lv(q)) == 1), mf = lv(q); else, break; end
        end
        n_inf = sum(sub.feasible == 0);
        fprintf('   %6.3f (%d)   ', mf, n_inf);
    end
    fprintf('\n');
end

% the error is not always monotone in the spread, which is why the
% envelope asks for every level up to the reported one
fprintf('\nsettings where the error falls again as the spread grows:\n');
for f = 1:numel(fams)
    idx_f = fam_all == fams{f};
    if ~any(idx_f), continue; end
    for x = Xs'
        sub = sortrows(T(idx_f & T.Xsep == x,:), 'mismatch');
        sub = sub(sub.feasible == 1,:);
        if height(sub) < 3, continue; end
        e = 100*abs(sub.dCCT_rel);
        de = diff(e);
        if any(de < -1e-9)
            k = find(de < -1e-9);
            fprintf('  %-26s X_sep=%.2f : at m = %s (by up to %.3f %%)\n', ...
                names{f}, x, mat2str(round(sub.mismatch(k+1),4)), max(-de(k)));
        end
    end
end

% F6 and F8 mismatch H, D, Imax and Pref by the same amount and differ
% only in which unit gets the larger current limit
fprintf('\nall four settings together against one at a time:\n');
fprintf('  %7s %8s %14s %14s %14s %14s\n', ...
    'X_sep','m','aligned (%)','adverse (%)','worst single','sum of singles');
for x = unique(T.Xsep(fam_all == "F6"))'
    for m = unique(T.mismatch(fam_all == "F6"))'
        singles = [];
        for f = {'F2','F3','F4','F5'}
            j = fam_all == f{1} & T.Xsep == x & abs(T.mismatch - m) < 1e-9 & feas;
            if any(j), singles(end+1) = 100*abs(T.dCCT_rel(find(j,1))); end
        end
        s6 = cell_val(T, fam_all, "F6", x, m);
        s8 = cell_val(T, fam_all, "F8", x, m);
        if isempty(singles)
            fprintf('  %7.2f %8.3f %14s %14s %14s %14s\n', x, m, s6, s8, '-', '-');
        else
            fprintf('  %7.2f %8.3f %14s %14s %14.3f %14.3f\n', ...
                x, m, s6, s8, max(singles), sum(singles));
        end
    end
end

% the ratios r = 2, 4, 10 of studyC_run are on this grid as
% m = 2(r-1)/(r+1); its data also hold 1/r, the same case with the units
% swapped
fprintf('\nagainst the Table XII data at the shared grid points:\n');
Tc = readtable(studyF_reference_csv(), 'VariableNamingRule','preserve');
c_case = string(Tc.("case"));
map = {'F2','C2','H',5; 'F3','C3','D',10; 'F4','C4','I',1.2};
fprintf('  %-24s %12s %14s %14s %10s\n','case','r','this sweep (%)','Table XII (%)','diff');
maxgap = 0;
for q = 1:size(map,1)
    for x = [0.05 0.20 0.70]
        for r = [2 4 10]
            m = 2*(r-1)/(r+1);
            lblF = sprintf('%s_%s_m%06.4f_X%05.3f', map{q,1}, map{q,3}, m, x);
            jF = find(string(T.("case")) == lblF, 1);
            if isempty(jF), continue; end
            best = NaN;
            for rr = [r 1/r]
                lblC = sprintf('%s_%s_%05.3f_X_%05.3f', map{q,2}, map{q,3}, rr, x);
                jC = find(c_case == lblC, 1);
                if ~isempty(jC) && ~isnan(Tc.("dCCT_rel")(jC))
                    best = 100*abs(Tc.("dCCT_rel")(jC));
                end
            end
            if isnan(best) || T.feasible(jF) ~= 1, continue; end
            vF = 100*abs(T.dCCT_rel(jF));
            fprintf('  %-24s %12g %14.4f %14.4f %10.4f\n', lblF, r, vF, best, vF-best);
            maxgap = max(maxgap, abs(vF-best));
        end
    end
end
fprintf('\n  largest difference: %.4f percentage points\n', maxgap);

all_env = [sep_env(:); env(:)];
writetable(struct2table(all_env, 'AsArray', true), 'studyF_admissibility_envelope.csv');
fprintf('\nwrote studyF_admissibility_envelope.csv (%d rows)\n', numel(all_env));
end

function s = cell_val(T, fam_all, fam, x, m)
j = fam_all == fam & T.Xsep == x & abs(T.mismatch - m) < 1e-9;
if ~any(j)
    s = sprintf('%14s','-');
elseif T.feasible(find(j,1)) ~= 1
    s = sprintf('%14s','no value');
else
    s = sprintf('%14.3f', 100*abs(T.dCCT_rel(find(j,1))));
end
end

function [vok, eok, vnext, enext] = walk(v, d, feasible, bound_pct)
% largest level such that all cases up to it have an error value within the
% bound, and the first level that does not, with its error. vok = 0 with
% eok = 0 means that not even the smallest level is within the bound.
[v, ord] = sort(v); d = d(ord); feasible = feasible(ord);
lv = unique(v);
vok = 0; eok = 0; vnext = NaN; enext = NaN;
for q = 1:numel(lv)
    sel = v <= lv(q);
    if all(feasible(sel) == 1) && max(100*abs(d(sel))) <= bound_pct
        vok = lv(q); eok = max(100*abs(d(sel)));
    else
        vnext = lv(q);
        sel_q = v == lv(q);
        if all(feasible(sel_q) == 1)
            enext = max(100*abs(d(sel_q)));
        else
            enext = NaN;
        end
        break;
    end
end
end
