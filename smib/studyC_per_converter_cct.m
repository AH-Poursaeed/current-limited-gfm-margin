function [CCT_per, info] = studyC_per_converter_cct(P, optsCCT)
% CCT of the two-converter hub by time-domain simulation. Both swing
% equations are integrated through the fault and for T_post seconds after
% clearing; the case is stable if neither converter moves more than
% instab_angle_rad (pi) away from its pre-fault angle. Bisection on the
% clearing time. No energy function and no aggregation enter here.
%   optsCCT.t_hi, .tol     search window and tolerance, s (2, 1e-3)
%   optsCCT.T_post         post-fault simulation time, s (10)
%   optsCCT.rtol, .atol    ode45 tolerances
% info.I_coh is the largest growth of the angle difference between the two
% converters during the fault at the CCT. CCT_per is NaN when there is no
% pre-fault operating point.

if nargin < 2 || isempty(optsCCT), optsCCT = struct(); end
if ~isfield(optsCCT,'t_hi'),    optsCCT.t_hi   = 2.0;  end
if ~isfield(optsCCT,'tol'),     optsCCT.tol    = 1e-3; end
if ~isfield(optsCCT,'T_post'),  optsCCT.T_post = 10.0; end
if ~isfield(optsCCT,'instab_angle_rad'), optsCCT.instab_angle_rad = pi; end
if ~isfield(optsCCT,'rtol'),    optsCCT.rtol   = 1e-7; end
if ~isfield(optsCCT,'atol'),    optsCCT.atol   = 1e-9; end

[delta1_op, delta2_op] = studyC_initial_sep(P);
if isnan(delta1_op) || isnan(delta2_op)
    CCT_per = NaN;
    info = struct('CCT_per',NaN,'N_ode45',0,'N_bisect',0, ...
                  'delta1_op',NaN,'delta2_op',NaN, ...
                  'I_coh',NaN,'peak_dw_diff',NaN, ...
                  'X_fault_traj',[], 'X_post_final',[], 't_post_final',[], ...
                  'max_dev_1',NaN,'max_dev_2',NaN, ...
                  'criterion','no pre-fault operating point');
    return;
end
sep_per = [delta1_op; delta2_op];
x0 = [delta1_op; 0; delta2_op; 0];

ws = P.ws;
odeF = @(t,x) per_conv_rhs(t, x, P, 'fault', ws);
odeP = @(t,x) per_conv_rhs(t, x, P, 'post',  ws);
odeopts = odeset('RelTol',optsCCT.rtol,'AbsTol',optsCCT.atol);

N_ode45 = 0;
N_bisect = 0;

    function [ok, info_local] = isStable(tc)
        n_local = 0;
        if tc <= 0
            x_clear = x0;
            t_fault_traj = 0;        X_fault_traj = x0(:).';
        else
            sol = ode45(odeF, [0 tc], x0, odeopts);
            n_local = n_local + 1;
            x_clear = sol.y(:,end);
            t_fault_traj = sol.x.'; X_fault_traj = sol.y.';
        end
        sol2 = ode45(odeP, [0 optsCCT.T_post], x_clear, odeopts);
        n_local = n_local + 1;
        X_post = sol2.y.';
        t_post = sol2.x.';
        delta1_t = X_post(:,1);
        delta2_t = X_post(:,3);
        max_dev_1 = max(abs(delta1_t - sep_per(1)));
        max_dev_2 = max(abs(delta2_t - sep_per(2)));
        ok = (max_dev_1 < optsCCT.instab_angle_rad) && ...
             (max_dev_2 < optsCCT.instab_angle_rad);
        info_local.x_clear      = x_clear;
        info_local.t_fault_traj = t_fault_traj;
        info_local.X_fault_traj = X_fault_traj;
        info_local.t_post       = t_post;
        info_local.X_post       = X_post;
        info_local.max_dev_1    = max_dev_1;
        info_local.max_dev_2    = max_dev_2;
        info_local.N_ode45_local = n_local;
    end

[ok_hi, info_hi] = isStable(optsCCT.t_hi);
N_ode45 = N_ode45 + info_hi.N_ode45_local;
if ok_hi
    CCT_per = optsCCT.t_hi;
else
    lo = 0.0; hi = optsCCT.t_hi;
    while (hi - lo) > optsCCT.tol
        mid = 0.5*(lo + hi);
        [ok_mid, info_mid] = isStable(mid);
        N_ode45 = N_ode45 + info_mid.N_ode45_local;
        N_bisect = N_bisect + 1;
        if ok_mid, lo = mid; else, hi = mid; end
    end
    CCT_per = lo;
end

% trajectory at the CCT, for I_coh and the speed difference
[~, info_lo] = isStable(CCT_per);
N_ode45 = N_ode45 + info_lo.N_ode45_local;

X_fault = info_lo.X_fault_traj;
if size(X_fault,1) < 2
    I_coh = abs(delta1_op - delta2_op);
    peak_dw_diff = 0;
else
    delta_diff   = X_fault(:,1) - X_fault(:,3);
    init_offset  = delta1_op - delta2_op;
    I_coh        = max(abs(delta_diff)) - abs(init_offset);
    peak_dw_diff = max(abs(X_fault(:,2) - X_fault(:,4)));
end

info.CCT_per       = CCT_per;
info.N_ode45       = N_ode45;
info.N_bisect      = N_bisect;
info.delta1_op     = delta1_op;
info.delta2_op     = delta2_op;
info.I_coh         = I_coh;
info.peak_dw_diff  = peak_dw_diff;
info.X_fault_traj  = X_fault;
info.X_post_final  = info_lo.X_post;
info.t_post_final  = info_lo.t_post;
info.max_dev_1     = info_lo.max_dev_1;
info.max_dev_2     = info_lo.max_dev_2;
info.criterion     = sprintf('max|delta_k - delta_k_SEP| < %.3f rad over T_post = %.1f s', ...
                             optsCCT.instab_angle_rad, optsCCT.T_post);
end

function dx = per_conv_rhs(~, x, P, stage, ws)
    d1 = x(1); dw1 = x(2);
    d2 = x(3); dw2 = x(4);
    [Pe1, Pe2] = studyC_per_converter_pe(d1, d2, P, stage);
    dx = [dw1;
          (ws/(2*P.H1))*(P.Pref1 - Pe1) - (P.D1/(2*P.H1))*dw1;
          dw2;
          (ws/(2*P.H2))*(P.Pref2 - Pe2) - (P.D2/(2*P.H2))*dw2];
end
