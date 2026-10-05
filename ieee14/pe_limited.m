function Pe = pe_limited(delta, Vth, Zth, Eabs, Zf, Imax)
% Electrical power at the GFM terminal for converter angle delta (relative
% to the Thevenin angle) with the circular current limit, eqs. (17)-(24).
% Zf is the converter reactance between the internal EMF and the terminal.

    theta_th = angle(Vth);
    delta_col = delta(:);

    E = Eabs .* exp(1j*(theta_th + delta_col));
    Iu = (E - Vth) ./ (Zth + Zf);

    I = Iu;
    Iabs_u = abs(Iu);
    lim = Iabs_u > Imax;
    I(lim) = Iu(lim) .* (Imax ./ Iabs_u(lim));

    Vpcc = Vth + Zth .* I;
    Pe_col = real(Vpcc .* conj(I));
    Pe = reshape(Pe_col, size(delta));
end
