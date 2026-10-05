function dx = fault_ode(~, x, Pm, H, D, omega_s, Vth, Zth, Eabs, Zf, Imax)
% Swing equation (14) on the fault-on curve. x = [delta; dw], delta in rad
% and dw the speed deviation in pu.
    delta = x(1);
    dw    = x(2);
    Pe = pe_limited(delta, Vth, Zth, Eabs, Zf, Imax);
    ddelta = omega_s * dw;
    ddw = (Pm - Pe - D*dw) / (2*H);
    dx = [ddelta; ddw];
end
