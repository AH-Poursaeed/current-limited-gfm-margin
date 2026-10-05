function [sep_wr, uep_unwrapped, sep_wr_out, uep_wr_out, reason] = select_sep_uep(eq_post, delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup)
% Post-fault SEP closest to delta0, and the first UEP met from it in the
% direction dir of the fault-on swing. The UEP is returned unwrapped, so
% that dir*(uep - sep) is the SEP-UEP distance of eq. (30).
    stable   = eq_post.roots(eq_post.K > 0);
    unstable = eq_post.roots(eq_post.K < 0);

    if isempty(stable)
        sep_wr=NaN; uep_unwrapped=NaN; sep_wr_out=NaN; uep_wr_out=NaN; reason='no_stable_root_post'; return;
    end
    if isempty(unstable)
        sep_wr=NaN; uep_unwrapped=NaN; sep_wr_out=NaN; uep_wr_out=NaN; reason='no_unstable_root_post'; return;
    end

    d0 = wrap2pi(delta0);
    st_wr = wrap2pi(stable);
    un_wr = wrap2pi(unstable);

    [st_wr, stable] = dedup_wrapped(st_wr, stable, wrap_dedup);
    [un_wr, unstable] = dedup_wrapped(un_wr, unstable, wrap_dedup);

    d = arrayfun(@(x) circDist(x, d0), st_wr);
    [~,i] = min(d);
    sep_wr = st_wr(i);

    best_alpha = inf; best_uep_wr = NaN;
    for r = un_wr(:).'
        alpha = dir*(r - sep_wr);
        while alpha <= 0
            alpha = alpha + 2*pi;
        end
        if alpha > min_alpha && alpha < best_alpha
            best_alpha = alpha;
            best_uep_wr = r;
        end
    end

    if isnan(best_uep_wr) || isinf(best_alpha)
        uep_unwrapped=NaN; reason='no_unstable_ahead_in_dir';
        sep_wr_out=sep_wr; uep_wr_out=NaN; return;
    end

    uep_unwrapped = sep_wr + dir*best_alpha;
    sep_wr_out = sep_wr;
    uep_wr_out = wrap2pi(uep_unwrapped);
    reason='ok';
end
