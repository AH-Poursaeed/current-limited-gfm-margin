function [U_post, Vcrit] = ias_energy_barrier(delta_h, Pe_post, Pref, delta_SEP, delta_UEP)
% Post-fault potential energy U(delta), zero at the SEP, eq. (33), and its
% value at the UEP, the critical barrier of eq. (35).

integrand = Pe_post - Pref;
U_raw = cumtrapz(delta_h, integrand);

[~, iS] = min(abs(delta_h - delta_SEP));
U_post = U_raw - U_raw(iS);

[~, iU] = min(abs(delta_h - delta_UEP));
Vcrit = U_post(iU);
end
