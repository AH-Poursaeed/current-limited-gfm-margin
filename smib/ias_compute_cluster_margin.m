function cluster = ias_compute_cluster_margin(cluster, deltaVec)
% Margin of one cluster with the circular limiter, eqs. (17)-(31):
% M_CL = Delta_lim / Delta_ref, the SEP-UEP separation with the current
% limit over the separation without it. deltaVec is the angle grid in rad
% (default 2000 points on [-pi, pi]). The curves, the equilibria and the
% margin are added to the input struct.

if nargin < 2 || isempty(deltaVec)
    Npts   = 2000;
    delta_h = linspace(-pi, +pi, Npts).';
else
    delta_h = deltaVec(:);
end

Eh   = cluster.Eh_mag;
Vth  = cluster.Vth_mag;
thth = cluster.theta_th;
Xnet = cluster.Xnet;
Imax = cluster.Imax;

Ephas  = Eh  .* exp(1j * delta_h);
Vth_ph = Vth .* exp(1j * thth);

% unlimited current and power
I0      = (Ephas - Vth_ph) ./ (1j * Xnet);
I0_mag  = abs(I0);

Vbus0   = Vth_ph + 1j * Xnet .* I0;

Pe0     = real(Vbus0 .* conj(I0));

% circular limiter
Ih      = I0;
overLim = I0_mag > Imax + 1e-12;
if any(overLim)
    Ih(overLim) = Imax .* (I0(overLim) ./ I0_mag(overLim));
end
Ih_mag  = abs(Ih);

Vbus    = Vth_ph + 1j * Xnet .* Ih;

Pe      = real(Vbus .* conj(Ih));

Pref = cluster.Pref_h;

[delta_SEP_ref, delta_UEP_ref, Delta_ref] = ...
    ias_find_equilibria(delta_h, Pe0, Pref);

[delta_SEP, delta_UEP, Delta_lim] = ...
    ias_find_equilibria(delta_h, Pe, Pref);

M_CL = Delta_lim / Delta_ref;

cluster.delta_h       = delta_h;
cluster.Pe_unlim      = Pe0;
cluster.Pe_lim        = Pe;
cluster.I_unlim       = I0_mag;
cluster.I_lim         = Ih_mag;

cluster.delta_SEP_ref = delta_SEP_ref;
cluster.delta_UEP_ref = delta_UEP_ref;
cluster.Delta_ref     = Delta_ref;

cluster.delta_SEP     = delta_SEP;
cluster.delta_UEP     = delta_UEP;
cluster.Delta_lim     = Delta_lim;

cluster.M_CL          = M_CL;

end
