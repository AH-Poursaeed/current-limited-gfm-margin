function [Pe1, Pe2, I1, I2, V_HUB, lim_state] = studyC_per_converter_pe(delta1, delta2, P, stage)
% Network solution of the two-converter hub for given internal angles
% (rad, relative to the Thevenin angle) and stage 'pre', 'fault' or 'post'.
% Returns the active power and current of each unit, the hub voltage and a
% flag per unit telling whether it sits on its current limit. A unit over
% its limit is replaced by a current source of magnitude Imax along the
% direction of its unlimited current, and the hub voltage is solved again.

% keeps the admittances finite when a unit is placed at the hub
eps_X = 1e-6;
Xa = max(P.Xa, eps_X);
Xb = max(P.Xb, eps_X);

switch lower(stage)
    case 'pre',   Vth_mag = P.Vth_pre;   Xc = P.Xc_pre;
    case 'fault', Vth_mag = P.Vth_fault; Xc = P.Xc_fault;
    case 'post',  Vth_mag = P.Vth_post;  Xc = P.Xc_post;
    otherwise, error('stage must be pre/fault/post');
end

E1c = P.E1 * exp(1j*(P.theta_th + delta1));
E2c = P.E2 * exp(1j*(P.theta_th + delta2));
Vth = Vth_mag * exp(1j*P.theta_th);

Ya = 1/(1j*Xa); Yb = 1/(1j*Xb); Yc = 1/(1j*Xc);

% both units as voltage sources
V_HUB_A = (E1c*Ya + E2c*Yb + Vth*Yc) / (Ya + Yb + Yc);
I1A = (E1c - V_HUB_A) * Ya;
I2A = (E2c - V_HUB_A) * Yb;

over1 = abs(I1A) > P.Imax1 + 1e-12;
over2 = abs(I2A) > P.Imax2 + 1e-12;

if ~over1 && ~over2
    I1 = I1A; I2 = I2A; V_HUB = V_HUB_A;
    lim_state = [false false];
elseif over1 && ~over2
    I1c = P.Imax1 * (I1A / abs(I1A));
    V_HUB = (I1c + E2c*Yb + Vth*Yc) / (Yb + Yc);
    I2B = (E2c - V_HUB) * Yb;
    if abs(I2B) <= P.Imax2 + 1e-12
        I1 = I1c; I2 = I2B; lim_state = [true false];
    else
        I2c = P.Imax2 * (I2B / abs(I2B));
        V_HUB = Vth + 1j*Xc*(I1c + I2c);
        I1 = I1c; I2 = I2c; lim_state = [true true];
    end
elseif ~over1 && over2
    I2c = P.Imax2 * (I2A / abs(I2A));
    V_HUB = (E1c*Ya + I2c + Vth*Yc) / (Ya + Yc);
    I1B = (E1c - V_HUB) * Ya;
    if abs(I1B) <= P.Imax1 + 1e-12
        I1 = I1B; I2 = I2c; lim_state = [false true];
    else
        I1c = P.Imax1 * (I1B / abs(I1B));
        V_HUB = Vth + 1j*Xc*(I1c + I2c);
        I1 = I1c; I2 = I2c; lim_state = [true true];
    end
else
    I1c = P.Imax1 * (I1A / abs(I1A));
    I2c = P.Imax2 * (I2A / abs(I2A));
    V_HUB = Vth + 1j*Xc*(I1c + I2c);
    I1 = I1c; I2 = I2c; lim_state = [true true];
end

Pe1 = real(V_HUB * conj(I1));
Pe2 = real(V_HUB * conj(I2));
end
