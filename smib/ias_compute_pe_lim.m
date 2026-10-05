function [Pe, I0_mag, Ih_mag] = ias_compute_pe_lim(cluster, delta_h, mode)
% Electrical power of a cluster with current limiting, eqs. (17)-(24).
%   mode   'circular' (default), 'd-priority' or 'q-priority'
% The priority rules split the unlimited current into the components in
% phase and in quadrature with the Thevenin voltage, eqs. (A.1)-(A.3).
% I0_mag and Ih_mag are the unlimited and limited current magnitudes.

if nargin < 3 || isempty(mode), mode = 'circular'; end

Eh   = cluster.Eh_mag;
Vth  = cluster.Vth_mag;
thth = cluster.theta_th;
Xnet = cluster.Xnet;
Imax = cluster.Imax;

Ephas  = Eh  .* exp(1j*delta_h);
Vth_ph = Vth .* exp(1j*thth);

I0     = (Ephas - Vth_ph) ./ (1j*Xnet);
I0_mag = abs(I0);

switch lower(mode)
    case 'circular'
        Ih      = I0;
        overLim = I0_mag > Imax + 1e-12;
        Ih(overLim) = Imax .* (I0(overLim) ./ I0_mag(overLim));
    case 'd-priority'
        % in-phase component first, what is left of the limit to the other
        I0_rot = I0 .* exp(-1j*thth);
        I_P    = real(I0_rot);
        I_Q    = imag(I0_rot);
        sgnP = sign(I_P);  sgnP(sgnP == 0) = 1;
        I_P_h = sgnP .* min(abs(I_P), Imax);
        headQ = sqrt(max(Imax.^2 - I_P_h.^2, 0));
        sgnQ = sign(I_Q);  sgnQ(sgnQ == 0) = 1;
        I_Q_h = sgnQ .* min(abs(I_Q), headQ);
        Ih = (I_P_h + 1j*I_Q_h) .* exp(1j*thth);
    case 'q-priority'
        % quadrature component first
        I0_rot = I0 .* exp(-1j*thth);
        I_P    = real(I0_rot);
        I_Q    = imag(I0_rot);
        sgnQ = sign(I_Q);  sgnQ(sgnQ == 0) = 1;
        I_Q_h = sgnQ .* min(abs(I_Q), Imax);
        headP = sqrt(max(Imax.^2 - I_Q_h.^2, 0));
        sgnP = sign(I_P);  sgnP(sgnP == 0) = 1;
        I_P_h = sgnP .* min(abs(I_P), headP);
        Ih = (I_P_h + 1j*I_Q_h) .* exp(1j*thth);
    otherwise
        error('ias_compute_pe_lim: unknown mode "%s"', mode);
end

Ih_mag = abs(Ih);

Vbus = Vth_ph + 1j*Xnet .* Ih;
Pe   = real(Vbus .* conj(Ih));
end
