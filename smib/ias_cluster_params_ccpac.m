function cluster = ias_cluster_params_ccpac(caseID)
% Parameters of the four single-converter cases of Table I, taken from
% reference [12] of the paper: a grid impedance of j0.4 pu in series with
% j0.3 pu (cases 1, 3, 4) or j0.1 pu (case 2). Vth = 1 pu, Pref = 0.5 pu,
% H = 5 s, D = 10 pu.

if nargin < 1
    caseID = 1;
end

Zg = 1j * 0.4;

switch caseID
    case 1
        Zv_mag = 0.3;
        Eh_mag = 1.0;
        Imax   = 1.2;
    case 2
        Zv_mag = 0.1;
        Eh_mag = 1.0;
        Imax   = 1.2;
    case 3
        Zv_mag = 0.3;
        Eh_mag = 1.5;
        Imax   = 1.2;
    case 4
        Zv_mag = 0.3;
        Eh_mag = 1.0;
        Imax   = 2.0;
    otherwise
        error('ias_cluster_params_ccpac: caseID must be 1..4.');
end

Zv = 1j * Zv_mag;
Zeq = Zg + Zv;

cluster.caseID   = caseID;

cluster.Eh_mag   = Eh_mag;
cluster.Vth_mag  = 1.0;
cluster.theta_th = 0.0;
cluster.Xnet     = imag(Zeq);

cluster.Imax     = Imax;

cluster.Hh       = 5.0;
cluster.Dh       = 10.0;

cluster.Pref_h   = 0.5;
cluster.Qref_h   = 0.0;

end
