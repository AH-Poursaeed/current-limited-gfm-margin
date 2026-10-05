% Fig. 8: post-fault potential energy of the cluster that sets the system
% CCT (bus 2) with its barrier, and the energy of the clearing state
% against the clearing time. The second curve meets the barrier at the
% CCT. The time marked in the figure is interpolated between two solver
% steps, so it agrees with the bisection result of step4 (Table X) to four
% decimals only.
% Needs STEP2_state.mat, STEP4_state.mat and ieee14_results.mat.
% Writes fig8_energy_barrier.png.

clear;

load('STEP2_state.mat',   'S');
load('ieee14_results.mat', 'rows');
load('STEP4_state.mat',   'tc_max');

set(groot,'defaultTextInterpreter',             'latex');
set(groot,'defaultLegendInterpreter',           'latex');
set(groot,'defaultAxesTickLabelInterpreter',    'latex');
set(groot,'defaultAxesFontName',                'Times New Roman');
set(groot,'defaultTextFontName',                'Times New Roman');
set(groot,'defaultAxesFontSize',                11);
set(groot,'defaultTextFontSize',                11);
set(groot,'defaultLineLineWidth',               1.5);
set(groot,'defaultAxesLineWidth',               0.75);
set(groot,'defaultAxesBox',                     'on');
set(groot,'defaultAxesXMinorTick',              'on');
set(groot,'defaultAxesYMinorTick',              'on');

C_Uh      = [0.13  0.47  0.71];
C_Vcrit   = [0.84  0.15  0.16];
C_Vcl     = [0.17  0.63  0.17];
C_tCCT    = [0.55  0.05  0.55];
C_safe    = [0.13  0.47  0.71];
C_unsafe  = [0.84  0.15  0.16];
C_grey    = [0.45  0.45  0.45];

% cluster with the smallest CCT
[~, idxLim] = min(rows.CCT_s);
b   = rows.bus(idxLim);
r   = rows(rows.bus == b, :);

P_h_ref = r.Pm_pu;
dir     = r.dir;
delta_h_SEP = r.sep_post;
alpha_u = r.alpha_u;

Vth_fault = S.Vth_fault_map(b);   Zth_fault = S.Zth_fault_map(b);
Vth_post  = S.Vth_post_map(b);    Zth_post  = S.Zth_post_map(b);
E_h       = S.Eabs_map(b);
I_h_max   = S.Imax_map(b);

H_h     = S.H_gfm;
D_h     = S.D_gfm;
omega_s = S.omega_s;

% potential energy along the SEP-UEP path, eq. (33), barrier at the UEP
Nalpha    = 14000;
alpha_h   = linspace(0, alpha_u, Nalpha).';
delta_path = delta_h_SEP + dir * alpha_h;

Pe_h_post = pe_limited(delta_path, Vth_post, Zth_post, E_h, S.Zf_gfm, I_h_max);

U_h = cumtrapz(alpha_h, (Pe_h_post - P_h_ref) * dir);
U_h = U_h - U_h(1);

V_h_crit = U_h(end);

% fault-on trajectory, and the post-fault energy of each of its states
% taken as the clearing state, eq. (34)
odefun = @(t,x) fault_ode(t, x, P_h_ref, H_h, D_h, omega_s, ...
                            Vth_fault, Zth_fault, E_h, S.Zf_gfm, I_h_max);
opts   = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t_ode, X_ode] = ode45(odefun, [0 tc_max], [r.delta0; 0], opts);

V_h_cl = nan(size(t_ode));

for k = 1:numel(t_ode)
    delta_h_cl = X_ode(k,1);
    omega_h_cl = X_ode(k,2);

    alpha_cl = dir * (delta_h_cl - delta_h_SEP);
    if alpha_cl < 0,        continue; end
    if alpha_cl >= alpha_u, V_h_cl(k) = NaN;            continue; end

    U_h_cl = interp1(alpha_h, U_h, alpha_cl, 'linear');

    V_h_cl(k) = H_h * omega_s * omega_h_cl^2 + U_h_cl;
end

% first crossing of the barrier
ok_idx   = isfinite(V_h_cl);
t_ok     = t_ode(ok_idx);
Vcl_ok   = V_h_cl(ok_idx);
dV       = Vcl_ok - V_h_crit;
cross_idx = find(dV >= 0, 1, 'first');

if isempty(cross_idx)
    t_CCT_h   = tc_max;
    CCT_stable = true;
elseif cross_idx == 1
    t_CCT_h   = 0;
    CCT_stable = false;
else
    t1 = t_ok(cross_idx-1);  t2 = t_ok(cross_idx);
    V1 = Vcl_ok(cross_idx-1); V2 = Vcl_ok(cross_idx);
    t_CCT_h   = t1 + (t2-t1)*(V_h_crit - V1)/(V2 - V1);
    CCT_stable = false;
end

fig = figure('Color','w','Units','inches','Position',[1 1 12.0 4.8]);

tl = tiledlayout(fig, 1, 2, ...
     'TileSpacing','compact', ...
     'Padding','compact');

ax1 = nexttile(tl,1);
hold(ax1,'on');

well_mask = U_h < V_h_crit;
aw = alpha_h(well_mask);  Uw = U_h(well_mask);
if ~isempty(aw)
    fill(ax1, [aw(1); aw; aw(end); aw(1)], ...
              [0;     Uw; 0;        0    ], ...
         C_safe,'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
end

over_mask = U_h >= V_h_crit;
ao = alpha_h(over_mask);  Uo = U_h(over_mask);
if ~isempty(ao)
    yhi_ax1  = V_h_crit * 1.20;
    fill(ax1, [ao(1); ao; ao(end); ao(1)], ...
        [V_h_crit; min(Uo, yhi_ax1); V_h_crit; V_h_crit], ...
        C_unsafe,'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
end

fill(ax1, [0; alpha_u; alpha_u; 0], ...
     V_h_crit + 0.003*abs(V_h_crit)*[-1 -1 1 1].', ...
     C_Vcrit,'FaceAlpha',0.22,'EdgeColor','none','HandleVisibility','off');

plot(ax1, delta_path, U_h, '-', 'Color',C_Uh, 'LineWidth',2.0, ...
     'DisplayName','$U_h(\delta_h)$');

hl_Vc1 = yline(ax1, V_h_crit, '--', 'Color',C_Vcrit, 'LineWidth',1.4, ...
               'Label',sprintf('$V_h^{\\mathrm{crit}} = %.4f$ pu', V_h_crit), ...
               'LabelHorizontalAlignment','right', ...
               'LabelVerticalAlignment','bottom', ...
               'Interpreter','latex','FontName','Times New Roman','FontSize',10);
hl_Vc1.Annotation.LegendInformation.IconDisplayStyle = 'off';

hl_z = yline(ax1, 0, ':', 'Color',C_grey,'LineWidth',0.8);
hl_z.Annotation.LegendInformation.IconDisplayStyle = 'off';

hl_SEP = xline(ax1, delta_h_SEP, ':', 'Color',C_grey, 'LineWidth',0.9, ...
               'Label','$\delta_h^{\mathrm{SEP,post}}$', ...
               'LabelOrientation','horizontal', ...
               'LabelVerticalAlignment','bottom', ...
               'Interpreter','latex','FontName','Times New Roman','FontSize',10);
hl_SEP.Annotation.LegendInformation.IconDisplayStyle = 'off';

hl_UEP = xline(ax1, delta_h_SEP + dir*alpha_u, ':', 'Color',C_Vcrit, 'LineWidth',0.9, ...
               'Label','$\delta_h^{\mathrm{UEP,post}}$', ...
               'LabelOrientation','horizontal', ...
               'LabelVerticalAlignment','bottom', ...
               'Interpreter','latex','FontName','Times New Roman','FontSize',10);
hl_UEP.Annotation.LegendInformation.IconDisplayStyle = 'off';

text(ax1, delta_h_SEP + dir*0.45*alpha_u, 0.07*V_h_crit, ...
     '$\leftarrow\Delta\delta_h$', ...
     'Interpreter','latex','FontName','Times New Roman','FontSize',10, ...
     'Color',C_grey);

text(ax1, alpha_u*0.28, V_h_crit*0.32, ...
     'Potential well', ...
     'Interpreter','latex','FontName','Times New Roman','FontSize',9.5, ...
     'Color',C_Uh*0.75,'HorizontalAlignment','center');

text(ax1, alpha_u*0.85, V_h_crit*1.12, ...
     'Barrier region', ...
     'Interpreter','latex','FontName','Times New Roman','FontSize',9.5, ...
     'Color',C_Vcrit*0.85,'HorizontalAlignment','center');

grid(ax1,'on');
ylim(ax1, [min(U_h) - 0.02*V_h_crit,  yhi_ax1]);
ax1.GridAlpha      = 0.16;
ax1.MinorGridAlpha = 0.07;
ax1.GridLineStyle  = ':';
ax1.Layer          = 'top';
ax1.TickDir        = 'out';

delta_h_UEP = delta_h_SEP + dir*alpha_u;
set(ax1,'XTick',[delta_h_SEP, delta_h_SEP + dir*(alpha_u/2), delta_h_UEP], ...
    'XTickLabel',{'$\delta_h^{\mathrm{SEP,post}}$', ...
    '$\delta_h^{\mathrm{SEP,post}}+\Delta\delta_h/2$', ...
    '$\delta_h^{\mathrm{UEP,post}}$'}, ...
    'TickLabelInterpreter','latex', ...
    'FontName','Times New Roman','FontSize',11,'LineWidth',0.75);

xlabel(ax1, '$\delta_h$ (rad)', ...
       'Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel(ax1, '$U_h(\delta_h)$ (pu)', ...
       'Interpreter','latex','FontName','Times New Roman','FontSize',12);
title(ax1, ...
    sprintf('Post-fault potential energy $U_h(\\delta_h)$, bus %d', b), ...
    'Interpreter','latex','FontName','Times New Roman', ...
    'FontSize',12,'FontWeight','normal');

lg1 = legend(ax1,'Location','northwest');
set(lg1,'Box','off','Interpreter','latex', ...
        'FontName','Times New Roman','FontSize',10.5);

ax2 = nexttile(tl,2);
hold(ax2,'on');

if ~CCT_stable
    smask = t_ok <= t_CCT_h & Vcl_ok < V_h_crit;
    if any(smask)
        ts = t_ok(smask);  Vs = Vcl_ok(smask);
        fill(ax2, [ts(1); ts; ts(end); ts(1)], ...
                  [0;     Vs; 0;        0    ], ...
             C_safe,'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
    end

    umask = Vcl_ok >= V_h_crit;
    if any(umask)
        tu = t_ok(umask);  Vu = Vcl_ok(umask);
        yhi_ax2  = max(Vcl_ok(isfinite(Vcl_ok))) * 1.10;
        fill(ax2, [tu(1); tu; tu(end); tu(1)], ...
            [V_h_crit; min(Vu, yhi_ax2); V_h_crit; V_h_crit], ...
            C_unsafe,'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
    end
end

fill(ax2, [t_ok(1); t_ok(end); t_ok(end); t_ok(1)], ...
     V_h_crit + 0.003*abs(V_h_crit)*[-1 -1 1 1].', ...
     C_Vcrit,'FaceAlpha',0.22,'EdgeColor','none','HandleVisibility','off');

hl_Vc2 = yline(ax2, V_h_crit, '--', 'Color',C_Vcrit, 'LineWidth',1.4, ...
               'Label',sprintf('$V_h^{\\mathrm{crit}} = %.4f$ pu', V_h_crit), ...
               'LabelHorizontalAlignment','right', ...
               'LabelVerticalAlignment','bottom', ...
               'Interpreter','latex','FontName','Times New Roman','FontSize',10);
hl_Vc2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(ax2, t_ok, Vcl_ok, '-', 'Color',C_Vcl, 'LineWidth',2.0, ...
     'DisplayName','$V_{h,\mathrm{cl}}(t_{\mathrm{cl}})$');

if ~CCT_stable
    hl_CCT = xline(ax2, t_CCT_h, ':', 'Color',C_tCCT, 'LineWidth',1.6, ...
        'Label',sprintf('$t^{\\star}_{\\mathrm{CCT},h} = %.5f~\\mathrm{s}$', t_CCT_h), ...
        'LabelOrientation','horizontal', ...
        'LabelVerticalAlignment','top', ...
        'Interpreter','latex','FontName','Times New Roman','FontSize',10);
    hl_CCT.Annotation.LegendInformation.IconDisplayStyle = 'off';

    plot(ax2, t_CCT_h, V_h_crit, 'o', ...
         'MarkerSize',7,'MarkerFaceColor',C_tCCT, ...
         'MarkerEdgeColor',C_tCCT*0.7,'LineWidth',1.0, ...
         'HandleVisibility','off');
end

if ~CCT_stable && t_CCT_h > t_ok(1) + 0.05*(t_ok(end)-t_ok(1))
    text(ax2, t_ok(1) + 0.35*t_CCT_h, min(Vcl_ok)*0.5 + V_h_crit*0.2, ...
         '$V_{h,\mathrm{cl}} < V_h^{\mathrm{crit}}$  (stable)', ...
         'Interpreter','latex','FontName','Times New Roman','FontSize',9.5, ...
         'Color',[0.08 0.45 0.12],'HorizontalAlignment','center');

    text(ax2, t_CCT_h + 0.30*(t_ok(end)-t_CCT_h), V_h_crit*1.10, ...
         '$V_{h,\mathrm{cl}} \geq V_h^{\mathrm{crit}}$  (unstable)', ...
         'Interpreter','latex','FontName','Times New Roman','FontSize',9.5, ...
         'Color',C_Vcrit*0.85,'HorizontalAlignment','center');
end

grid(ax2,'on');
ylim(ax2, [0,  yhi_ax2]);
ax2.GridAlpha      = 0.16;
ax2.MinorGridAlpha = 0.07;
ax2.GridLineStyle  = ':';
ax2.Layer          = 'top';
ax2.TickDir        = 'out';

set(ax2,'FontName','Times New Roman','FontSize',11, ...
        'TickLabelInterpreter','latex','LineWidth',0.75);

xlabel(ax2, 'Clearing time $t_{\mathrm{cl}}$ (s)', ...
    'Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel(ax2, '$V_{h,\mathrm{cl}}(t_{\mathrm{cl}})$ (pu)', ...
    'Interpreter','latex','FontName','Times New Roman','FontSize',12);
title(ax2, ...
    'Clearing-time verdict: $V_{h,\mathrm{cl}}(t_{\mathrm{cl}})$ vs.\ $V_h^{\mathrm{crit}}$', ...
    'Interpreter','latex','FontName','Times New Roman', ...
    'FontSize',12,'FontWeight','normal');

lg2 = legend(ax2,'Location','northwest');
set(lg2,'Box','off','Interpreter','latex', ...
        'FontName','Times New Roman','FontSize',10.5);

drawnow;
exportgraphics(fig, 'fig8_energy_barrier.png', 'Resolution', 300);

% take the plot defaults set above off again
props = {'TextInterpreter','LegendInterpreter','AxesTickLabelInterpreter', ...
         'AxesFontName','TextFontName','AxesFontSize','TextFontSize', ...
         'LineLineWidth','AxesLineWidth','AxesBox','AxesXMinorTick','AxesYMinorTick'};
for k = 1:numel(props)
    set(groot, ['default' props{k}], 'remove');
end
