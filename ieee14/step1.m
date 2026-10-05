% Step 1: scenario of the modified IEEE 14-bus system, pre-fault power flow
% and the operating-point quantities the later steps need.
% Requires MATPOWER. Writes STEP1_state.mat.

clear; close all;
define_constants;

S = struct();

% synchronous generator at bus 1, GFM converters at 2, 3 and 8, GFL at 6
S.mpc_name   = 'case14';
S.sg_bus     = 1;
S.gfm_buses  = [2 3 8];
S.gfl_buses  = 6;

% three-phase fault at bus 4, cleared by tripping line 2-4
S.fault_bus  = 4;
S.Zfault_pu  = 1e-3 + 1j*1e-3;
S.trip_pairs = [2 4];

% series reactance of the generator and of the GFM converters, pu
S.Xdprime_sg = 0.25;
S.Xf_gfm     = 0.12;

% swing parameters of the GFM converters
S.f0       = 60;
S.omega_s  = 2*pi*S.f0;
S.H_gfm    = 0.50;
S.D_gfm    = 0.05;

% current limit as a multiple of the rated current
S.kI       = 1.25;
S.Vrated   = 1.0;

S.Srated_MVA = containers.Map('KeyType','double','ValueType','double');
S.Srated_MVA(2) = 100;
S.Srated_MVA(3) = 100;
S.Srated_MVA(8) = 100;

% active power in MW. Buses 3, 6 and 8 are synchronous condensers in
% case14 and are given 10 MW each here.
S.Pg_set_MW = containers.Map('KeyType','double','ValueType','double');
S.Pg_set_MW(2) = 40;
S.Pg_set_MW(3) = 10;
S.Pg_set_MW(6) = 10;
S.Pg_set_MW(8) = 10;

mpc = loadcase(S.mpc_name);
S.mpc0 = mpc;

gen_buses = mpc.gen(:, GEN_BUS).';
assert(ismember(S.sg_bus, gen_buses), 'SG bus must be in mpc.gen(:,GEN_BUS).');
assert(all(ismember(S.gfm_buses, gen_buses)), 'Some GFM buses are not generator buses.');
assert(all(ismember(S.gfl_buses, gen_buses)), 'Some GFL buses are not generator buses.');
assert(isempty(intersect(S.gfm_buses, S.gfl_buses)), 'A bus cannot be both GFM and GFL.');

for k = 1:size(mpc.gen,1)
    b = mpc.gen(k, GEN_BUS);
    if isKey(S.Pg_set_MW, b)
        Pg_new = S.Pg_set_MW(b);
        Pmin = mpc.gen(k, PMIN);
        Pmax = mpc.gen(k, PMAX);
        mpc.gen(k, PG) = min(max(Pg_new, Pmin), Pmax);
    end
end

mpopt = mpoption('verbose', 0, 'out.all', 0);
[res, success] = runpf(mpc, mpopt);
assert(success == 1, 'Power flow did not converge.');

S.baseMVA = res.baseMVA;
S.nb      = size(res.bus,1);

S.Vop = res.bus(:, VM) .* exp(1j * deg2rad(res.bus(:, VA)));

fprintf('\nStep 1: pre-fault power flow\n');
fprintf('|V| range: [%.4f, %.4f] pu\n', min(abs(S.Vop)), max(abs(S.Vop)));

[Ybus, ~, ~] = makeYbus(S.baseMVA, res.bus, res.branch);
S.Ybus_pre = Ybus;

S.Sgen_bus = zeros(S.nb,1);
for g = 1:size(res.gen,1)
    b = res.gen(g, GEN_BUS);
    S.Sgen_bus(b) = S.Sgen_bus(b) + (res.gen(g, PG) + 1j*res.gen(g, QG)) / S.baseMVA;
end

Pd = res.bus(:, PD) / S.baseMVA;
Qd = res.bus(:, QD) / S.baseMVA;
S.Sload = Pd + 1j*Qd;

% loads as constant admittances
S.Yload = conj(S.Sload) ./ max(1e-12, abs(S.Vop).^2);
S.Yload_mat = spdiags(S.Yload, 0, S.nb, S.nb);

S.Zs_sg  = 1j*S.Xdprime_sg;     S.Ys_sg  = 1/S.Zs_sg;
S.Zf_gfm = 1j*S.Xf_gfm;         S.Yf_gfm = 1/S.Zf_gfm;

S.Eabs_map  = containers.Map('KeyType','double','ValueType','double');
S.Pm_map    = containers.Map('KeyType','double','ValueType','double');

% internal EMF and power reference of each GFM from the power flow
for b = S.gfm_buses
    Iop = conj(S.Sgen_bus(b) ./ S.Vop(b));
    Eop = S.Vop(b) + S.Zf_gfm * Iop;
    S.Eabs_map(b) = abs(Eop);
    S.Pm_map(b)   = real(S.Sgen_bus(b));
end

fprintf('\nGFM converters, pu on system base:\n');
for b = S.gfm_buses
    fprintf('  bus %d: Pm=%.6f pu | |Eop|=%.6f pu\n', b, S.Pm_map(b), S.Eabs_map(b));
end

S.Irated_map = containers.Map('KeyType','double','ValueType','double');
S.Imax_map   = containers.Map('KeyType','double','ValueType','double');
S.Sr_pu_map  = containers.Map('KeyType','double','ValueType','double');

for b = S.gfm_buses
    if isKey(S.Srated_MVA, b)
        Sr_pu = S.Srated_MVA(b) / S.baseMVA;
    else
        gidx = find(res.gen(:, GEN_BUS) == b, 1);
        Pmax = res.gen(gidx, PMAX) / S.baseMVA;
        Qmax = max(abs(res.gen(gidx, QMAX)), abs(res.gen(gidx, QMIN))) / S.baseMVA;
        Sr_pu = sqrt(Pmax^2 + Qmax^2);
    end
    Irated = Sr_pu / S.Vrated;
    Imax   = S.kI * Irated;

    S.Sr_pu_map(b)  = Sr_pu;
    S.Irated_map(b) = Irated;
    S.Imax_map(b)   = Imax;

    fprintf('  bus %d: Sr=%.6f pu | Irated=%.6f | Imax=%.6f (kI=%.2f)\n', ...
        b, Sr_pu, Irated, Imax, S.kI);
end

save('STEP1_state.mat','S');
fprintf('\nsaved STEP1_state.mat\n');
