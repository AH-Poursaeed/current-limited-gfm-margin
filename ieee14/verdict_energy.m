function ok = verdict_energy(tc, t, X, sep, dir, alpha_u, alpha_grid, U, Ecrit, H, omega_s)
% True if the post-fault energy of the clearing state, eq. (34), does not
% exceed the barrier Ecrit of eq. (35), which is the test of eq. (36). The
% clearing state is interpolated on the fault-on trajectory (t, X). alpha
% is the angle from the SEP along the swing direction, on which the
% potential U is tabulated up to the UEP at alpha_u; a state beyond the
% UEP is unstable. With the speed deviation dw in pu and time in s the
% speed term of (34) is H*omega_s*dw^2.
    if tc <= 0
        deltac = X(1,1);
        dwc    = X(1,2);
    else
        deltac = interp1(t, X(:,1), tc, 'linear');
        dwc    = interp1(t, X(:,2), tc, 'linear');
    end

    alpha_c = dir*(deltac - sep);

    % behind the SEP no potential energy is counted
    if alpha_c < 0
        alpha_c = 0;
    end

    if alpha_c >= alpha_u
        ok = false; return;
    end

    Uc = interp1(alpha_grid, U, alpha_c, 'linear');
    if isnan(Uc)
        ok = false; return;
    end

    Eclear = (H*omega_s)*(dwc^2) + Uc;
    ok = (Eclear <= Ecrit);
end
