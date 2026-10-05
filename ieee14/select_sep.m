function [delta0, delta0_raw, reason] = select_sep(eq_pre, delta_guess, wrap2pi, circDist, wrap_dedup)
% Stable equilibrium closest to delta_guess on the circle. delta0 is
% wrapped to [0, 2*pi), delta0_raw is the root as found on the grid.
    stable = eq_pre.roots(eq_pre.K > 0);
    if isempty(stable)
        delta0 = NaN; delta0_raw = NaN; reason='no_stable_root'; return;
    end

    dg = wrap2pi(delta_guess);
    stab_wr = wrap2pi(stable);

    [stab_wr, stable] = dedup_wrapped(stab_wr, stable, wrap_dedup);

    d = arrayfun(@(x) circDist(x, dg), stab_wr);
    [~,i] = min(d);
    delta0 = stab_wr(i);
    delta0_raw = stable(i);
    reason='ok';
end
