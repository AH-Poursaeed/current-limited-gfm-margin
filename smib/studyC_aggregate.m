function [agg, CCT_agg, M_agg] = studyC_aggregate(P, optsCCT)
% Aggregated single-cluster model of the two-converter hub, eqs. (7)-(13):
% inertia, damping, power references and current limits are summed and the
% cluster EMF is taken from the pre-fault hub voltage, eq. (A.6). Returns
% the aggregate with its CCT and margin from the single-cluster routines.

if nargin < 2 || isempty(optsCCT)
    optsCCT = struct('t_hi',2.0,'tol',1e-3,'f0',P.f0);
end
if ~isfield(optsCCT,'f0'), optsCCT.f0 = P.f0; end

[delta1_op, delta2_op] = studyC_initial_sep(P);
if isnan(delta1_op) || isnan(delta2_op)
    error('studyC_aggregate:noOperatingPoint', ...
          'no pre-fault operating point found for the two converters');
end
[~, ~, ~, ~, V_HUB_op, ~] = studyC_per_converter_pe(delta1_op, delta2_op, P, 'pre');

E_h_lifted     = abs(V_HUB_op);
theta_h_lifted = angle(V_HUB_op) - P.theta_th;

H_h     = P.H1 + P.H2;
D_h     = P.D1 + P.D2;
Pref_h  = P.Pref1 + P.Pref2;
Qref_h  = P.Qref1 + P.Qref2;
Imax_h  = P.Imax1 + P.Imax2;

cluster.caseID   = 0;
cluster.Eh_mag   = E_h_lifted;
cluster.Vth_mag  = P.Vth;
cluster.theta_th = P.theta_th;
cluster.Xnet     = P.Xc;
cluster.Imax     = Imax_h;
cluster.Hh       = H_h;
cluster.Dh       = D_h;
cluster.Pref_h   = Pref_h;
cluster.Qref_h   = Qref_h;

% same angle grid as the single-converter cases
Npts = 2000;
deltaGrid = linspace(-pi, +pi, Npts).';
cluster = ias_compute_cluster_margin(cluster, deltaGrid);
M_agg = cluster.M_CL;

optsCCT.delta_h = deltaGrid;
faultModel.kV = P.kV;
faultModel.kX = P.kX;
out = ias_compute_cct_smib(cluster, faultModel, optsCCT);
if ~isfield(out, 'Vcrit')
    error('studyC_aggregate:noEquilibrium', ...
          'aggregated cluster has no post-fault SEP/UEP pair');
end
CCT_agg = out.CCT;

agg.H_h     = H_h;
agg.D_h     = D_h;
agg.Pref_h  = Pref_h;
agg.Qref_h  = Qref_h;
agg.Imax_h  = Imax_h;
agg.E_h     = E_h_lifted;
agg.theta_h = theta_h_lifted;
agg.V_HUB_op = V_HUB_op;
agg.delta1_op = delta1_op;
agg.delta2_op = delta2_op;
agg.cluster = cluster;
agg.smib_out = out;
agg.Vcrit = out.Vcrit;
agg.delta_SEP_pre  = out.delta_SEP_pre;
agg.delta_SEP_post = out.delta_SEP_post;
agg.delta_UEP_post = out.delta_UEP_post;
end
