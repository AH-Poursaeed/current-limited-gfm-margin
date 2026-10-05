function [delta_SEP, delta_UEP, Delta_delta] = ...
         ias_find_equilibria(delta, Pe, Pref)
% Equilibria of a power-angle curve: intersections of Pe(delta) with Pref,
% eqs. (25)-(26), classified by the sign of the slope, eqs. (28)-(29).
% Returns the SEP and the UEP closest to zero angle and their separation,
% eq. (30), or NaN if either is missing. Each equilibrium is given at the
% grid angle nearest to the interpolated root.

F = Pe - Pref;
N = numel(delta);

% sign changes of Pe - Pref, root by linear interpolation
idx_bracket = find(F(1:N-1) .* F(2:N) <= 0);

roots_delta = [];

for k = 1:numel(idx_bracket)
    i  = idx_bracket(k);
    d1 = delta(i);
    d2 = delta(i+1);
    F1 = F(i);
    F2 = F(i+1);

    if abs(F1 - F2) < 1e-12
        d_root = 0.5 * (d1 + d2);
    else
        d_root = d1 - F1 * (d2 - d1) / (F2 - F1);
    end

    roots_delta(end+1,1) = d_root;
end

if isempty(roots_delta)
    warning('ias_find_equilibria: no intersection with Pref found.');
    delta_SEP   = NaN;
    delta_UEP   = NaN;
    Delta_delta = NaN;
    return;
end

dPe_dDelta = gradient(Pe, delta(2) - delta(1));

SEP_candidates = [];
UEP_candidates = [];

for k = 1:numel(roots_delta)
    d_root = roots_delta(k);
    [~, idx_near] = min(abs(delta - d_root));
    slope = dPe_dDelta(idx_near);

    if slope > 0
        SEP_candidates(end+1,1) = delta(idx_near);
    elseif slope < 0
        UEP_candidates(end+1,1) = delta(idx_near);
    end
end

if isempty(SEP_candidates) || isempty(UEP_candidates)
    warning('ias_find_equilibria: could not find both SEP and UEP.');
    delta_SEP   = NaN;
    delta_UEP   = NaN;
    Delta_delta = NaN;
    return;
end

[~, idx_s] = min(abs(SEP_candidates));
[~, idx_u] = min(abs(UEP_candidates));

delta_SEP   = SEP_candidates(idx_s);
delta_UEP   = UEP_candidates(idx_u);
Delta_delta = abs(delta_UEP - delta_SEP);

end
