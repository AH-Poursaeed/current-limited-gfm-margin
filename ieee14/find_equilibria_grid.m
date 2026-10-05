function eq = find_equilibria_grid(delta, Pe, Pm, fmin_tol, dedup_tol)
% Roots of Pe(delta) = Pm on the angle grid: sign changes, located by
% linear interpolation, and near-tangent points where |Pe - Pm| has a local
% minimum below fmin_tol. eq.roots holds the angles and eq.K the slope
% dPe/ddelta at each of them (K > 0 stable, K < 0 unstable), eqs. (25)-(29).
    y = Pe - Pm;

    idx_sc = find(y(1:end-1).*y(2:end) <= 0);
    roots_sc = nan(size(idx_sc));
    for k = 1:numel(idx_sc)
        i = idx_sc(k);
        y1 = y(i); y2 = y(i+1);
        d1 = delta(i); d2 = delta(i+1);
        if abs(y2 - y1) < 1e-15
            roots_sc(k) = 0.5*(d1+d2);
        else
            roots_sc(k) = d1 - y1*(d2-d1)/(y2-y1);
        end
    end

    absf = abs(y);
    idx_min = find(absf(2:end-1) <= absf(1:end-2) & absf(2:end-1) <= absf(3:end));
    idx_min = idx_min + 1;
    roots_tg = [];
    for k = 1:numel(idx_min)
        i = idx_min(k);
        if absf(i) < fmin_tol
            a = delta(max(1,i-2));
            c = delta(min(numel(delta), i+2));
            try
                xstar = fminbnd(@(x) abs(interp1(delta, y, x, 'linear', 'extrap')), a, c);
                if abs(interp1(delta, y, xstar, 'linear', 'extrap')) < fmin_tol
                    roots_tg(end+1,1) = xstar;
                end
            catch
            end
        end
    end

    roots = [roots_sc(:); roots_tg(:)];
    roots = roots(~isnan(roots));
    if isempty(roots)
        eq.roots = []; eq.K = []; return;
    end

    roots = sort(roots);
    roots_u = roots(1);
    for k = 2:numel(roots)
        if abs(roots(k)-roots_u(end)) > dedup_tol
            roots_u(end+1,1) = roots(k);
        end
    end
    roots = roots_u;

    % slope at the grid point nearest to each root
    dPe = gradient(Pe, delta);
    K = nan(size(roots));
    for k = 1:numel(roots)
        j = round(interp1(delta, 1:numel(delta), roots(k), 'nearest', 'extrap'));
        j = max(1, min(numel(delta), j));
        K(k) = dPe(j);
    end

    eq.roots = roots;
    eq.K = K;
end
