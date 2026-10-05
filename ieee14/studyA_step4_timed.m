function studyA_step4_timed()
% step4 with timers around its parts: equilibria, energy barrier, fault-on
% integration and bisection, per cluster (Table XI, reduced-order row).
% Wall-clock times depend on the machine and on whether MATLAB has already
% run the code once. Needs STEP3_state.mat. Writes studyA_step4_timing.csv.

load('STEP3_state.mat','S');

tc_max = 10.0;
rtol   = 1e-7;
atol   = 1e-9;
tol_tc = 1e-6;
max_it = 70;

Nscan       = 8001;
delta_grid  = linspace(-2*pi, 2*pi, Nscan).';
fmin_tol    = 1e-4;
root_dedup  = 1e-6;
wrap_dedup  = 1e-6;
min_alpha   = 1e-7;
eq_tol      = 1e-6;

wrap2pi  = @(x) mod(x, 2*pi);
circDist = @(a,b) abs(atan2(sin(a-b), cos(a-b)));

nb_g = numel(S.gfm_buses);
T_total_h   = nan(nb_g,1);
T_eq_h      = nan(nb_g,1);
T_barrier_h = nan(nb_g,1);
T_ode45_h   = nan(nb_g,1);
T_bisect_h  = nan(nb_g,1);
N_ode45_h   = zeros(nb_g,1);
N_bisect_h  = zeros(nb_g,1);
CCT_h_t     = nan(nb_g,1);

tic_total = tic;

for ii = 1:nb_g
    b = S.gfm_buses(ii);
    Pm   = S.Pm_map(b);
    Eabs = S.Eabs_map(b);
    Vth_pre   = S.Vth_pre_map(b);   Zth_pre   = S.Zth_pre_map(b);
    Vth_fault = S.Vth_fault_map(b); Zth_fault = S.Zth_fault_map(b);
    Vth_post  = S.Vth_post_map(b);  Zth_post  = S.Zth_post_map(b);
    Imax = S.Imax_map(b);

    t_bus = tic;

    t0 = tic;
    Pe_pre  = pe_limited(delta_grid, Vth_pre,  Zth_pre,  Eabs, S.Zf_gfm, Imax);
    Pe_post = pe_limited(delta_grid, Vth_post, Zth_post, Eabs, S.Zf_gfm, Imax);
    eq_pre  = find_equilibria_grid(delta_grid, Pe_pre,  Pm, fmin_tol, root_dedup);
    eq_post = find_equilibria_grid(delta_grid, Pe_post, Pm, fmin_tol, root_dedup);
    [delta0, ~, ~] = select_sep(eq_pre, S.delta_g_map(b), wrap2pi, circDist, wrap_dedup);
    if isnan(delta0)
        T_total_h(ii) = toc(t_bus);
        continue;
    end

    Pe_fault_0 = pe_limited(delta0, Vth_fault, Zth_fault, Eabs, S.Zf_gfm, Imax);
    acc0 = Pm - Pe_fault_0; dir = sign(acc0); if dir==0, dir=+1; end

    Pe_post_max = max(Pe_post); Pe_post_min = min(Pe_post);
    if ~(Pm <= Pe_post_max + eq_tol && Pm >= Pe_post_min - eq_tol)
        CCT_h_t(ii) = 0.0;
        T_eq_h(ii) = toc(t0);
        T_total_h(ii) = toc(t_bus);
        continue;
    end

    [sep, uep, ~, ~, ~] = select_sep_uep(eq_post, delta0, dir, wrap2pi, circDist, min_alpha, wrap_dedup);
    T_eq_h(ii) = toc(t0);
    if isnan(sep) || isnan(uep)
        T_total_h(ii) = toc(t_bus);
        continue;
    end

    t0 = tic;
    alpha_u = dir*(uep - sep);
    while alpha_u <= 0, alpha_u = alpha_u + 2*pi; end
    alpha_grid = linspace(0, alpha_u, 10001).';
    delta_path = sep + dir*alpha_grid;
    Pe_post_path = pe_limited(delta_path, Vth_post, Zth_post, Eabs, S.Zf_gfm, Imax);
    g = (Pe_post_path - Pm) * dir;
    U = cumtrapz(alpha_grid, g); U = U - U(1);
    Ecrit = U(end);
    T_barrier_h(ii) = toc(t0);

    t0 = tic;
    odefun = @(tt,xx) fault_ode(tt, xx, Pm, S.H_gfm, S.D_gfm, S.omega_s, Vth_fault, Zth_fault, Eabs, S.Zf_gfm, Imax);
    opts = odeset('RelTol',rtol,'AbsTol',atol);
    [tt, X] = ode45(odefun, [0 tc_max], [delta0; 0], opts);
    N_ode45_h(ii) = 1;
    T_ode45_h(ii) = toc(t0);

    t0 = tic;
    verdict_ok = @(tc) verdict_energy(tc, tt, X, sep, dir, alpha_u, alpha_grid, U, Ecrit, S.H_gfm, S.omega_s);
    n_iter = 0;
    ok_max = verdict_ok(tc_max);
    if ok_max
        CCT = tc_max;
    else
        lo = 0.0; hi = tc_max;
        if ~verdict_ok(lo)
            CCT = 0.0;
        else
            for it = 1:max_it
                n_iter = n_iter + 1;
                if (hi - lo) <= tol_tc, break; end
                mid = 0.5*(lo+hi);
                if verdict_ok(mid), lo = mid; else, hi = mid; end
            end
            CCT = lo;
        end
    end
    T_bisect_h(ii) = toc(t0);
    N_bisect_h(ii) = n_iter;
    CCT_h_t(ii) = CCT;
    T_total_h(ii) = toc(t_bus);
end

T_total_all = toc(tic_total);

fprintf('\nTiming of the CCT computation, s\n');
fprintf('MATLAB %s\n', version);
fprintf('tc_max=%.3f s | tol_tc=%.1e s | RelTol=%.1e | AbsTol=%.1e\n', tc_max, tol_tc, rtol, atol);
for ii = 1:nb_g
    fprintf('bus %d: tot=%.4f | eq=%.4f | bar=%.4f | ode45=%.4f | bis=%.4f | Node45=%d | Nbis=%d | CCT=%.6f\n', ...
        S.gfm_buses(ii), T_total_h(ii), T_eq_h(ii), T_barrier_h(ii), ...
        T_ode45_h(ii), T_bisect_h(ii), N_ode45_h(ii), N_bisect_h(ii), CCT_h_t(ii));
end
fprintf('all %d clusters: %.4f s\n', nb_g, T_total_all);

bus_vec = S.gfm_buses(:);
T = table(bus_vec, T_total_h, T_eq_h, T_barrier_h, T_ode45_h, T_bisect_h, ...
    N_ode45_h, N_bisect_h, CCT_h_t, ...
    'VariableNames', {'bus','t_total_s','t_eq_s','t_barrier_s','t_ode45_s', ...
                      't_bisect_s','n_ode45','n_bisect','CCT_s'});
writetable(T, 'studyA_step4_timing.csv');
end
