function [delta_SEP, delta_UEP, Delta_delta, ok] = ias_find_equilibria_fwd(delta, Pe, Pref)
% As ias_find_equilibria, except that the UEP is the first one above the
% SEP. With a resistive coupling the curve is no longer symmetric and the
% UEP closest to zero angle can lie below the SEP. ok is false when no
% such pair exists.

F = Pe - Pref;
N = numel(delta);

idx_bracket = find(F(1:N-1) .* F(2:N) <= 0);

roots_delta = [];
for k = 1:numel(idx_bracket)
    i  = idx_bracket(k);
    d1 = delta(i);   d2 = delta(i+1);
    F1 = F(i);       F2 = F(i+1);
    if abs(F1 - F2) < 1e-12
        d_root = 0.5 * (d1 + d2);
    else
        d_root = d1 - F1 * (d2 - d1) / (F2 - F1);
    end
    roots_delta(end+1,1) = d_root;
end

if isempty(roots_delta)
    delta_SEP = NaN; delta_UEP = NaN; Delta_delta = NaN; ok = false;
    return;
end

dPe_dDelta = gradient(Pe, delta(2) - delta(1));

SEP_candidates = [];
UEP_candidates = [];
for k = 1:numel(roots_delta)
    [~, idx_near] = min(abs(delta - roots_delta(k)));
    slope = dPe_dDelta(idx_near);
    if slope > 0
        SEP_candidates(end+1,1) = delta(idx_near);
    elseif slope < 0
        UEP_candidates(end+1,1) = delta(idx_near);
    end
end

if isempty(SEP_candidates) || isempty(UEP_candidates)
    delta_SEP = NaN; delta_UEP = NaN; Delta_delta = NaN; ok = false;
    return;
end

[~, idx_s] = min(abs(SEP_candidates));
delta_SEP  = SEP_candidates(idx_s);

fwd = UEP_candidates(UEP_candidates > delta_SEP);
if isempty(fwd)
    delta_UEP = NaN; Delta_delta = NaN; ok = false;
    return;
end
delta_UEP   = min(fwd);
Delta_delta = abs(delta_UEP - delta_SEP);
ok = true;
end
