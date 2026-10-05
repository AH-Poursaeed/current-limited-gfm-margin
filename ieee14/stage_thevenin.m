function [Vth, Zth] = stage_thevenin(Ynet, Yload_mat, nb, ...
    sg_bus, Zs_sg, Ys_sg, ...
    gfl_buses, gfm_buses, ...
    Sgen_bus, Vop, Zf_gfm, Yf_gfm, ...
    tbus)
% Thevenin equivalent of the network seen from GFM terminal tbus for one
% stage admittance matrix, eqs. (A.4)-(A.5). The generator and the other
% GFM converters enter as Norton sources behind their reactances and the
% GFL converter as a current injection, all at the pre-fault operating
% point; loads are constant admittances.

    Yphys = Ynet + Yload_mat;
    Iphys = zeros(nb,1);

    Iop_sg = conj(Sgen_bus(sg_bus) ./ Vop(sg_bus));
    Eop_sg = Vop(sg_bus) + Zs_sg * Iop_sg;
    In_sg  = Eop_sg / Zs_sg;

    Yphys(sg_bus, sg_bus) = Yphys(sg_bus, sg_bus) + Ys_sg;
    Iphys(sg_bus) = Iphys(sg_bus) + In_sg;

    for bb = gfl_buses
        Iphys(bb) = Iphys(bb) + conj(Sgen_bus(bb) ./ Vop(bb));
    end

    for bb = gfm_buses
        if bb == tbus, continue; end
        Iop = conj(Sgen_bus(bb) ./ Vop(bb));
        Eop = Vop(bb) + Zf_gfm * Iop;
        In  = Eop / Zf_gfm;
        Yphys(bb,bb) = Yphys(bb,bb) + Yf_gfm;
        Iphys(bb)    = Iphys(bb)    + In;
    end

    % Kron reduction onto the terminal
    k = tbus;
    r = setdiff(1:nb, k);

    Yrr = Yphys(r,r);  Yrt = Yphys(r,k);
    Ytr = Yphys(k,r);  Ytt = Yphys(k,k);

    Yeq = Ytt - Ytr * (Yrr \ Yrt);
    Zth = 1 / Yeq;

    Ir  = Iphys(r);
    Vth = -(1 / Yeq) * (Ytr * (Yrr \ Ir));
end
