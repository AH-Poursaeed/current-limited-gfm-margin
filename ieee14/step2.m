% Step 2: admittance matrices of the pre-fault, fault-on and post-fault
% network and the Thevenin equivalent seen from each GFM terminal in each
% stage (Table VII). Writes STEP2_state.mat and thevenin_equivalents.csv.

clear; close all;
define_constants;

load('STEP1_state.mat','S');

wrapPi = @(x) atan2(sin(x), cos(x));

Ypre = S.Ybus_pre;

% fault-on: shunt at the faulted bus; post-fault: tripped line removed
Yfault = Ypre;
k = S.fault_bus;
Yfault(k,k) = Yfault(k,k) + 1/S.Zfault_pu;

mpc_post = S.mpc0;
br = mpc_post.branch;

for kk = 1:size(S.trip_pairs,1)
    a = S.trip_pairs(kk,1);
    b = S.trip_pairs(kk,2);
    idx = find((br(:,F_BUS)==a & br(:,T_BUS)==b) | (br(:,F_BUS)==b & br(:,T_BUS)==a), 1);
    assert(~isempty(idx), 'Requested trip line %d-%d not found in case.', a, b);
    br(idx, BR_STATUS) = 0;
end
mpc_post.branch = br;

[Ypost, ~, ~] = makeYbus(S.baseMVA, mpc_post.bus, mpc_post.branch);

S.Ybus_fault = Yfault;
S.Ybus_post  = Ypost;

S.Vth_pre_map   = containers.Map('KeyType','double','ValueType','any');
S.Zth_pre_map   = containers.Map('KeyType','double','ValueType','any');
S.Vth_fault_map = containers.Map('KeyType','double','ValueType','any');
S.Zth_fault_map = containers.Map('KeyType','double','ValueType','any');
S.Vth_post_map  = containers.Map('KeyType','double','ValueType','any');
S.Zth_post_map  = containers.Map('KeyType','double','ValueType','any');

for tbus = S.gfm_buses
    [Vth_pre,  Zth_pre]  = stage_thevenin(Ypre,   S.Yload_mat, S.nb, ...
        S.sg_bus, S.Zs_sg, S.Ys_sg, S.gfl_buses, S.gfm_buses, ...
        S.Sgen_bus, S.Vop, S.Zf_gfm, S.Yf_gfm, tbus);

    [Vth_f,    Zth_f]    = stage_thevenin(Yfault, S.Yload_mat, S.nb, ...
        S.sg_bus, S.Zs_sg, S.Ys_sg, S.gfl_buses, S.gfm_buses, ...
        S.Sgen_bus, S.Vop, S.Zf_gfm, S.Yf_gfm, tbus);

    [Vth_post, Zth_post] = stage_thevenin(Ypost,  S.Yload_mat, S.nb, ...
        S.sg_bus, S.Zs_sg, S.Ys_sg, S.gfl_buses, S.gfm_buses, ...
        S.Sgen_bus, S.Vop, S.Zf_gfm, S.Yf_gfm, tbus);

    S.Vth_pre_map(tbus)   = Vth_pre;    S.Zth_pre_map(tbus)   = Zth_pre;
    S.Vth_fault_map(tbus) = Vth_f;      S.Zth_fault_map(tbus) = Zth_f;
    S.Vth_post_map(tbus)  = Vth_post;   S.Zth_post_map(tbus)  = Zth_post;
end

% converter angle relative to the pre-fault Thevenin angle, used later to
% pick the operating SEP among the equilibria
S.delta_g_map = containers.Map('KeyType','double','ValueType','double');

for b = S.gfm_buses
    Iop = conj(S.Sgen_bus(b) ./ S.Vop(b));
    Eop = S.Vop(b) + S.Zf_gfm * Iop;
    Vth_pre = S.Vth_pre_map(b);
    S.delta_g_map(b) = wrapPi(angle(Eop) - angle(Vth_pre));
end

fprintf('\nStep 2: Thevenin equivalents per stage\n');
fprintf('fault at bus %d, Zf = %.3e%+.3ej pu, line %d-%d tripped\n', S.fault_bus, ...
    real(S.Zfault_pu), imag(S.Zfault_pu), S.trip_pairs(1,1), S.trip_pairs(1,2));

for b = S.gfm_buses
    Vpre = S.Vth_pre_map(b);   Zpre = S.Zth_pre_map(b);
    Vf   = S.Vth_fault_map(b); Zf   = S.Zth_fault_map(b);
    Vpo  = S.Vth_post_map(b);  Zpo  = S.Zth_post_map(b);

    fprintf('\nbus %d\n', b);
    fprintf('  pre   : Vth = %.6f pu, %8.3f deg | Zth = %.6f %+.6fj pu\n', abs(Vpre), rad2deg(angle(Vpre)), real(Zpre), imag(Zpre));
    fprintf('  fault : Vth = %.6f pu, %8.3f deg | Zth = %.6f %+.6fj pu\n', abs(Vf),   rad2deg(angle(Vf)),   real(Zf),   imag(Zf));
    fprintf('  post  : Vth = %.6f pu, %8.3f deg | Zth = %.6f %+.6fj pu\n', abs(Vpo),  rad2deg(angle(Vpo)),  real(Zpo),  imag(Zpo));
end

% Table VII as a csv file
stage = {'pre', 'fault', 'post'};
tab = {};
for b = S.gfm_buses
    V = [S.Vth_pre_map(b), S.Vth_fault_map(b), S.Vth_post_map(b)];
    Z = [S.Zth_pre_map(b), S.Zth_fault_map(b), S.Zth_post_map(b)];
    for s = 1:3
        tab(end+1,:) = {b, stage{s}, abs(V(s)), rad2deg(angle(V(s))), real(Z(s)), imag(Z(s))};
    end
end
writetable(cell2table(tab, 'VariableNames', ...
    {'bus', 'stage', 'Vth_pu', 'Vth_deg', 'Rth_pu', 'Xth_pu'}), 'thevenin_equivalents.csv');

save('STEP2_state.mat','S');
fprintf('\nsaved STEP2_state.mat and thevenin_equivalents.csv\n');
