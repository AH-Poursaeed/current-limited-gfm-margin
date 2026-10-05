function cluster = studyD_compute_cluster_margin(cluster, deltaVec, mode)
% ias_compute_cluster_margin for a chosen limiter mode: 'circular',
% 'd-priority' or 'q-priority'. The reference separation is always that of
% the unlimited curve.

if nargin < 2 || isempty(deltaVec)
    Npts = 2000; deltaVec = linspace(-pi, +pi, Npts).';
end
if nargin < 3 || isempty(mode), mode = 'circular'; end
delta_h = deltaVec(:);

[Pe_lim, I_unlim_mag, Ih_mag] = ias_compute_pe_lim(cluster, delta_h, mode);

% unlimited curve: a limit that never binds
cluster_ref       = cluster;
cluster_ref.Imax  = 1e6;
Pe_unlim          = ias_compute_pe_lim(cluster_ref, delta_h, 'circular');

[delta_SEP_ref, delta_UEP_ref, Delta_ref] = ias_find_equilibria(delta_h, Pe_unlim, cluster.Pref_h);
[delta_SEP,     delta_UEP,     Delta_lim] = ias_find_equilibria(delta_h, Pe_lim,   cluster.Pref_h);

if isnan(Delta_lim) || isnan(Delta_ref) || Delta_ref < 1e-12
    M_CL = NaN;
else
    M_CL = Delta_lim / Delta_ref;
end

cluster.delta_h        = delta_h;
cluster.Pe_unlim       = Pe_unlim;
cluster.Pe_lim         = Pe_lim;
cluster.I_unlim        = I_unlim_mag;
cluster.I_lim          = Ih_mag;
cluster.delta_SEP_ref  = delta_SEP_ref;
cluster.delta_UEP_ref  = delta_UEP_ref;
cluster.Delta_ref      = Delta_ref;
cluster.delta_SEP      = delta_SEP;
cluster.delta_UEP      = delta_UEP;
cluster.Delta_lim      = Delta_lim;
cluster.M_CL           = M_CL;
cluster.limiter_mode   = mode;
end
