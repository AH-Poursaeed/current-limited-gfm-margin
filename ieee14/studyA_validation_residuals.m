function studyA_validation_residuals()
% Checks the stage Thevenin equivalents against the full network (Table
% VIII). For each stage and GFM terminal:
%   - terminal voltage for four test currents, full nodal solution against
%     Vth + Zth*I
%   - electrical power at seven converter angles: the limited converter
%     current, found with the equivalent, is injected into the full
%     network and the power at the terminal is compared with the reduced
%     expression
% Needs STEP2_state.mat. Writes studyA_validation_residuals.csv.

load('STEP2_state.mat','S');

Itests = [ ...
    0.20 + 0.10j
    0.20 - 0.10j
    0.35 + 0.00j
    0.10 + 0.25j];

delta_test = deg2rad([-120 -60 -20 0 20 60 120]).';

stages = {'pre','fault','post'};
buses  = S.gfm_buses(:);

nS = numel(stages);
nB = numel(buses);

maxAbsV  = zeros(nS, nB);
maxRelV  = zeros(nS, nB);
maxAbsPe = zeros(nS, nB);

for ss = 1:nS
    st = stages{ss};
    switch st
        case 'pre',   Ynet=S.Ybus_pre;   Vth_map=S.Vth_pre_map;   Zth_map=S.Zth_pre_map;
        case 'fault', Ynet=S.Ybus_fault; Vth_map=S.Vth_fault_map; Zth_map=S.Zth_fault_map;
        case 'post',  Ynet=S.Ybus_post;  Vth_map=S.Vth_post_map;  Zth_map=S.Zth_post_map;
    end

    for bi = 1:nB
        b = buses(bi);
        Vth = Vth_map(b);
        Zth = Zth_map(b);
        Eabs = S.Eabs_map(b);
        Imax = S.Imax_map(b);

        eA = zeros(numel(Itests),1);
        eR = zeros(numel(Itests),1);
        for k = 1:numel(Itests)
            It = Itests(k);
            Vfull = terminal_voltage_full(Ynet, S.Yload_mat, S.nb, ...
                S.sg_bus, S.Zs_sg, S.Ys_sg, S.gfl_buses, S.gfm_buses, ...
                S.Sgen_bus, S.Vop, S.Zf_gfm, S.Yf_gfm, b, It);
            Vpred = Vth + Zth*It;
            eA(k) = abs(Vfull - Vpred);
            eR(k) = eA(k) / max(1e-12, abs(Vfull));
        end
        maxAbsV(ss,bi) = max(eA);
        maxRelV(ss,bi) = max(eR);

        ePe = zeros(numel(delta_test),1);
        for k = 1:numel(delta_test)
            d = delta_test(k);
            Pe_red  = pe_reduced(d, Vth, Zth, Eabs, S.Zf_gfm, Imax);
            Pe_full = pe_full_network(d, Ynet, S.Yload_mat, S.nb, ...
                S.sg_bus, S.Zs_sg, S.Ys_sg, S.gfl_buses, S.gfm_buses, ...
                S.Sgen_bus, S.Vop, S.Zf_gfm, S.Yf_gfm, b, Vth, Zth, Eabs, Imax);
            ePe(k) = abs(Pe_full - Pe_red);
        end
        maxAbsPe(ss,bi) = max(ePe);
    end
end

fprintf('\nThevenin equivalents against the full network\n');
fprintf('voltage : max |V_full - (Vth + Zth*I)| over %d test currents\n', numel(Itests));
fprintf('power   : max |Pe_full - Pe_reduced| over %d test angles\n', numel(delta_test));
for ss = 1:nS
    fprintf('\n-- %s --\n', stages{ss});
    for bi = 1:nB
        fprintf('bus %d : maxAbsV=%.3e pu | maxRelV=%.3e | maxAbsPe=%.3e pu\n', ...
            buses(bi), maxAbsV(ss,bi), maxRelV(ss,bi), maxAbsPe(ss,bi));
    end
end
fprintf('\nlargest over all stages and buses: |dV| %.3e pu, relative %.3e, |dPe| %.3e pu\n', ...
    max(maxAbsV(:)), max(maxRelV(:)), max(maxAbsPe(:)));

rows = {};
for ss = 1:nS
    for bi = 1:nB
        rows(end+1,:) = {stages{ss}, buses(bi), maxAbsV(ss,bi), maxRelV(ss,bi), maxAbsPe(ss,bi)};
    end
end
T = cell2table(rows, 'VariableNames', {'stage','bus','maxAbsV_pu','maxRelV','maxAbsPe_pu'});
writetable(T, 'studyA_validation_residuals.csv');
end

% Full-network side of the two checks. The rest of the system is modelled
% as in stage_thevenin; the terminal carries the test current, or the
% current of the limited converter.

function Vt_full = terminal_voltage_full(Ynet, Yload_mat, nb, ...
    sg_bus, Zs_sg, Ys_sg, gfl_buses, gfm_buses, Sgen_bus, Vop, Zf_gfm, Yf_gfm, tbus, It_test)
    Yphys = Ynet + Yload_mat;
    Iphys = zeros(nb,1);
    Iop_sg = conj(Sgen_bus(sg_bus) ./ Vop(sg_bus));
    Eop_sg = Vop(sg_bus) + Zs_sg*Iop_sg;
    In_sg  = Eop_sg / Zs_sg;
    Yphys(sg_bus, sg_bus) = Yphys(sg_bus, sg_bus) + Ys_sg;
    Iphys(sg_bus) = Iphys(sg_bus) + In_sg;
    for bb = gfl_buses
        Iphys(bb) = Iphys(bb) + conj(Sgen_bus(bb) ./ Vop(bb));
    end
    for bb = gfm_buses
        if bb == tbus, continue; end
        Iop = conj(Sgen_bus(bb) ./ Vop(bb));
        Eop = Vop(bb) + Zf_gfm*Iop;
        In  = Eop / Zf_gfm;
        Yphys(bb,bb) = Yphys(bb,bb) + Yf_gfm;
        Iphys(bb)    = Iphys(bb)    + In;
    end
    Iinj = Iphys; Iinj(tbus) = It_test;
    V = Yphys \ Iinj;
    Vt_full = V(tbus);
end

function Pe = pe_reduced(delta, Vth, Zth, Eabs, Zf, Imax)
    theta_th = angle(Vth);
    E0 = Eabs * exp(1j*(theta_th + delta));
    Iu = (E0 - Vth) / (Zth + Zf);
    if abs(Iu) > Imax
        I = Iu * (Imax/abs(Iu));
    else
        I = Iu;
    end
    Vpcc = Vth + Zth*I;
    Pe = real(Vpcc * conj(I));
end

function Pe = pe_full_network(delta, Ynet, Yload_mat, nb, ...
    sg_bus, Zs_sg, Ys_sg, gfl_buses, gfm_buses, ...
    Sgen_bus, Vop, Zf_gfm, Yf_gfm, tbus, Vth_stage, Zth_stage, Eabs, Imax)
    [Yphys, Iphys] = build_rest_system(Ynet, Yload_mat, nb, ...
        sg_bus, Zs_sg, Ys_sg, gfl_buses, gfm_buses, Sgen_bus, Vop, Zf_gfm, Yf_gfm, tbus);
    theta_th = angle(Vth_stage);
    E0 = Eabs * exp(1j*(theta_th + delta));
    Iu = (E0 - Vth_stage) / (Zth_stage + Zf_gfm);
    if abs(Iu) > Imax
        Ik = Iu * (Imax/abs(Iu));
    else
        Ik = Iu;
    end
    Iinj = Iphys; Iinj(tbus) = Iinj(tbus) + Ik;
    V = Yphys \ Iinj;
    Vk = V(tbus);
    Pe = real(Vk * conj(Ik));
end

function [Yphys, Iphys] = build_rest_system(Ynet, Yload_mat, nb, ...
    sg_bus, Zs_sg, Ys_sg, gfl_buses, gfm_buses, Sgen_bus, Vop, Zf_gfm, Yf_gfm, tbus)
    Yphys = Ynet + Yload_mat;
    Iphys = zeros(nb,1);
    Iop_sg = conj(Sgen_bus(sg_bus) ./ Vop(sg_bus));
    Eop_sg = Vop(sg_bus) + Zs_sg*Iop_sg;
    In_sg  = Eop_sg / Zs_sg;
    Yphys(sg_bus, sg_bus) = Yphys(sg_bus, sg_bus) + Ys_sg;
    Iphys(sg_bus) = Iphys(sg_bus) + In_sg;
    for bb = gfl_buses
        Iphys(bb) = Iphys(bb) + conj(Sgen_bus(bb) ./ Vop(bb));
    end
    for bb = gfm_buses
        if bb == tbus, continue; end
        Iop = conj(Sgen_bus(bb) ./ Vop(bb));
        Eop = Vop(bb) + Zf_gfm*Iop;
        In  = Eop / Zf_gfm;
        Yphys(bb,bb) = Yphys(bb,bb) + Yf_gfm;
        Iphys(bb)    = Iphys(bb)    + In;
    end
end
