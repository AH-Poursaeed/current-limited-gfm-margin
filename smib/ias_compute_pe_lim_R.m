function [Pe, I0_mag, Ih_mag, I_P_h, I_Q_h] = ias_compute_pe_lim_R(cluster, delta_h, mode)
% ias_compute_pe_lim with two additions:
%   cluster.Rnet    series resistance of the coupling (default 0)
%   cluster.frame   axis of the priority split: 'thevenin' (default), the
%                   Thevenin voltage angle as in (A.1), or 'converter',
%                   the converter angle delta
% With Rnet = 0 and the default frame the values are those of
% ias_compute_pe_lim. I_P_h and I_Q_h are the limited current components
% on the chosen axis.

if nargin < 3 || isempty(mode), mode = 'circular'; end

Eh   = cluster.Eh_mag;
Vth  = cluster.Vth_mag;
thth = cluster.theta_th;
Xnet = cluster.Xnet;
Imax = cluster.Imax;
if isfield(cluster, 'Rnet') && ~isempty(cluster.Rnet)
    Rnet = cluster.Rnet;
else
    Rnet = 0;
end

if isfield(cluster, 'frame') && ~isempty(cluster.frame)
    frame = lower(cluster.frame);
else
    frame = 'thevenin';
end

Ephas  = Eh  .* exp(1j*delta_h);
Vth_ph = Vth .* exp(1j*thth);

switch frame
    case 'thevenin',  ang = thth .* ones(size(delta_h));
    case 'converter', ang = delta_h;
    otherwise, error('ias_compute_pe_lim_R: unknown frame "%s"', frame);
end

% without resistance keep the expression of ias_compute_pe_lim, so the
% two functions agree to the last bit
if Rnet == 0
    I0 = (Ephas - Vth_ph) ./ (1j*Xnet);
else
    I0 = (Ephas - Vth_ph) ./ (Rnet + 1j*Xnet);
end
I0_mag = abs(I0);

switch lower(mode)
    case 'circular'
        Ih      = I0;
        overLim = I0_mag > Imax + 1e-12;
        Ih(overLim) = Imax .* (I0(overLim) ./ I0_mag(overLim));
    case 'd-priority'
        I0_rot = I0 .* exp(-1j*ang);
        I_P    = real(I0_rot);
        I_Q    = imag(I0_rot);
        sgnP = sign(I_P);  sgnP(sgnP == 0) = 1;
        I_P_h = sgnP .* min(abs(I_P), Imax);
        headQ = sqrt(max(Imax.^2 - I_P_h.^2, 0));
        sgnQ = sign(I_Q);  sgnQ(sgnQ == 0) = 1;
        I_Q_h = sgnQ .* min(abs(I_Q), headQ);
        Ih = (I_P_h + 1j*I_Q_h) .* exp(1j*ang);
    case 'q-priority'
        I0_rot = I0 .* exp(-1j*ang);
        I_P    = real(I0_rot);
        I_Q    = imag(I0_rot);
        sgnQ = sign(I_Q);  sgnQ(sgnQ == 0) = 1;
        I_Q_h = sgnQ .* min(abs(I_Q), Imax);
        headP = sqrt(max(Imax.^2 - I_Q_h.^2, 0));
        sgnP = sign(I_P);  sgnP(sgnP == 0) = 1;
        I_P_h = sgnP .* min(abs(I_P), headP);
        Ih = (I_P_h + 1j*I_Q_h) .* exp(1j*ang);
    otherwise
        error('ias_compute_pe_lim_R: unknown mode "%s"', mode);
end

Ih_mag = abs(Ih);

if Rnet == 0
    Vbus = Vth_ph + 1j*Xnet .* Ih;
else
    Vbus = Vth_ph + (Rnet + 1j*Xnet) .* Ih;
end
Pe = real(Vbus .* conj(Ih));

Ih_rot = Ih .* exp(-1j*ang);
I_P_h  = real(Ih_rot);
I_Q_h  = imag(Ih_rot);
end
