% Single-converter cases 1-4 of Section IV.A: SEP-UEP separations and
% margins (Table I, Figs. 3-5), then energy barrier and CCT (Table II).
% Writes ias_ccpac_margins.csv, ias_ccpac_cct.csv and three png files.

clear; close all;

fsAxes   = 12;
fsLabel  = 14;
fsTitle  = 14;
fsLegend = 12;

set(groot, ...
    'defaultAxesFontSize',           fsAxes, ...
    'defaultTextFontSize',           fsAxes, ...
    'defaultAxesLabelFontSizeMultiplier', fsLabel/fsAxes, ...
    'defaultAxesTitleFontSizeMultiplier', fsTitle/fsAxes, ...
    'defaultLegendFontSize',         fsLegend, ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultTextInterpreter',        'latex', ...
    'defaultLegendInterpreter',      'latex');

%% margins
Npts   = 2000;
deltaGrid = linspace(-pi, +pi, Npts).';

numCases = 4;

param_h     = ias_cluster_params_ccpac(1);
clusters(1) = ias_compute_cluster_margin(param_h, deltaGrid);
fprintf('Case 1:  M_CL = %.3f\n', clusters(1).M_CL);

for h = 2:numCases
    param_h     = ias_cluster_params_ccpac(h);
    clusters(h) = ias_compute_cluster_margin(param_h, deltaGrid);
    fprintf('Case %d:  M_CL = %.3f\n', clusters(h).caseID, clusters(h).M_CL);
end

sysRes = ias_system_margin(clusters);
fprintf('\nsmallest margin: %.3f (case %d)\n\n', ...
        sysRes.M_sys, clusters(sysRes.idx_crit).caseID);

deg = @(x) x * 180/pi;

%% Fig. 3: power-angle curves with and without the current limit
fig3 = figure('Color','w','Position',[100 100 1100 700]);

for k = 1:numCases
    c = clusters(k);

    subplot(2,2,k); hold on; box on;

    h_unlim = plot(deg(c.delta_h), c.Pe_unlim, 'k--', 'LineWidth', 1.2);
    h_lim   = plot(deg(c.delta_h), c.Pe_lim,   'b-',  'LineWidth', 1.8);
    h_pref  = yline(c.Pref_h, 'r--', 'LineWidth', 1.2);

    h_sep = plot(deg(c.delta_SEP), c.Pref_h, 'go', ...
        'MarkerSize', 6, 'LineWidth', 1.2, 'MarkerFaceColor', 'g');
    h_uep = plot(deg(c.delta_UEP), c.Pref_h, 'ro', ...
        'MarkerSize', 6, 'LineWidth', 1.2, 'MarkerFaceColor', 'r');

    Delta_ref_deg = deg(c.Delta_ref);
    Delta_lim_deg = deg(c.Delta_lim);

    yl = get(gca,'YLim');

    x1_lim = deg(c.delta_SEP);
    x2_lim = deg(c.delta_UEP);
    y_lim  = c.Pref_h;

    % SEP-UEP distance on the limited curve, and below it on the
    % reference curve
    darrow(x1_lim, x2_lim, y_lim, [0 0.5 0], 1.2);

    xc_lim = (x1_lim + x2_lim)/2;
    yc_lim = y_lim + 0.06*diff(yl);
    text(xc_lim, yc_lim, '$\Delta\delta_h$', ...
        'Interpreter','latex', 'HorizontalAlignment','center', ...
        'Color', [0 0.5 0], 'FontSize', 9);

    x1_ref = deg(c.delta_SEP_ref);
    x2_ref = deg(c.delta_UEP_ref);
    y_ref  = c.Pref_h - 0.08*diff(yl);

    darrow(x1_ref, x2_ref, y_ref, [0.3 0.8 0.3], 1.0);

    xc_ref = (x1_ref + x2_ref)/2;
    yc_ref = y_ref - 0.06*diff(yl);
    text(xc_ref, yc_ref, '$\Delta\delta_h^{\mathrm{ref}}$', ...
        'Interpreter','latex', 'HorizontalAlignment','center', ...
        'Color', [0.3 0.8 0.3], 'FontSize', 9);

    x_txt = 0.55;
    y_txt1 = 0.24;
    y_txt2 = 0.14;

    text(x_txt, y_txt1, ...
        sprintf('$\\Delta\\delta_h^{\\mathrm{ref}} = %.1f^{\\circ}$', Delta_ref_deg), ...
        'Units','normalized', 'Interpreter','latex', ...
        'Color', [0.3 0.8 0.3], 'FontSize', 9);

    text(x_txt, y_txt2, ...
        sprintf('$\\Delta\\delta_h = %.1f^{\\circ}$', Delta_lim_deg), ...
        'Units','normalized', 'Interpreter','latex', ...
        'Color', [0 0.5 0], 'FontSize', 9);

    xlabel('$\delta_h~[\mathrm{deg}]$','Interpreter','latex');
    if k == 1 || k == 3
        ylabel('$P_{e,h}(\delta_h)~[\mathrm{pu}]$','Interpreter','latex');
    end

    title(sprintf('Case %d: $M_h^{\\mathrm{CL}} = %.2f$', ...
          c.caseID, c.M_CL), 'Interpreter','latex');

    if k == 1
        lg = legend([h_unlim, h_lim, h_pref, h_sep, h_uep], ...
               {'$P_{e,h}^0(\delta_h)$', '$P_{e,h}(\delta_h)$', ...
                '$P_h^{\mathrm{ref}}$', 'SEP', 'UEP'}, ...
               'Location','NorthWest');
        set(lg,'Interpreter','latex');
    end
end

%% Fig. 4: current magnitude, shaded where the limiter is inactive
fig4 = figure('Color','w','Position',[120 120 1100 700]);

for k = 1:numCases
    c = clusters(k);

    subplot(2,2,k); hold on; box on;

    h_unlim = plot(deg(c.delta_h), c.I_unlim, 'k--', 'LineWidth', 1.4);
    h_lim   = plot(deg(c.delta_h), c.I_lim,   'b-',  'LineWidth', 1.8);
    h_imax  = yline(c.Imax, 'r--', 'LineWidth', 1.4);

    idx_unlim = c.I_unlim <= c.Imax + 1e-6;
    if any(idx_unlim)
        xOmega = deg(c.delta_h(idx_unlim));
        yOmega = c.I_unlim(idx_unlim);
        xx = [xOmega(1); xOmega; xOmega(end)];
        yy = [0;        yOmega; 0];
        h_omega = patch(xx, yy, [0.9 0.9 1.0], ...
                        'EdgeColor','none', 'FaceAlpha',0.4);
        uistack(h_omega,'bottom');
    else
        h_omega = patch(nan,nan,[0.9 0.9 1.0], ...
                        'EdgeColor','none','FaceAlpha',0.4);
    end

    xlabel('$\delta_h~[\mathrm{deg}]$');
    if k == 1 || k == 3
        ylabel('$|I_h(\delta_h)|~[\mathrm{pu}]$');
    end
    title(sprintf('Case %d: cluster current magnitude', c.caseID));

    if k == 1
        legend([h_unlim, h_lim, h_imax, h_omega], ...
               {'$|I_h^0(\delta_h)|$', '$|I_h(\delta_h)|$', ...
                '$I_h^{\max}$', '$\Omega_h^{\mathrm{unlim}}$'}, ...
               'Location','SouthWest');
    end
end

Xnet = arrayfun(@(c)c.Xnet, clusters).';
Eh   = arrayfun(@(c)c.Eh_mag, clusters).';
Imax = arrayfun(@(c)c.Imax,  clusters).';
Mcl  = sysRes.M_vec;

%% Fig. 5: margin against Xnet, |E| and Imax
fig5 = figure('Color','w','Position',[140 140 900 500]);

subplot(3,1,1); hold on; box on;
idx_X = [1 2];
plot(Xnet(idx_X), Mcl(idx_X), 'o-', 'LineWidth', 1.6, 'MarkerSize',8);
xlabel('$X_h^{\mathrm{net}}~[\mathrm{pu}]$');
ylabel('$M_h^{\mathrm{CL}}$');
title('Effect of cluster reactance $X_h^{\mathrm{net}}$ (Cases 1 \& 2)');
grid on;

subplot(3,1,2); hold on; box on;
idx_E = [1 3];
plot(Eh(idx_E), Mcl(idx_E), 's-', 'LineWidth', 1.6, 'MarkerSize',8);
xlabel('$|E_h|~[\mathrm{pu}]$');
ylabel('$M_h^{\mathrm{CL}}$');
title('Effect of EMF magnitude $|E_h|$ (Cases 1 \& 3)');
grid on;

subplot(3,1,3); hold on; box on;
idx_I = [1 4];
plot(Imax(idx_I), Mcl(idx_I), 'd-', 'LineWidth', 1.6, 'MarkerSize',8);
xlabel('$I_h^{\max}~[\mathrm{pu}]$');
ylabel('$M_h^{\mathrm{CL}}$');
title('Effect of current limit $I_h^{\max}$ (Cases 1 \& 4)');
grid on;

%% Table I
caseVec      = (1:numCases).';
XnetVec      = Xnet;
EhVec        = Eh;
ImaxVec      = Imax;
DeltaRef_deg = deg([clusters.Delta_ref].');
DeltaLim_deg = deg([clusters.Delta_lim].');

T = table(caseVec, XnetVec, EhVec, ImaxVec, ...
          DeltaRef_deg, DeltaLim_deg, Mcl, ...
    'VariableNames', {'Case', 'Xnet_pu', 'Eh_pu', 'Imax_pu', ...
                      'DeltaRef_deg', 'DeltaLim_deg', 'M_CL'});

disp(T);
writetable(T, 'ias_ccpac_margins.csv');

%% Table II: energy barrier and CCT, fault-on voltage at 0.2 pu
faultModel.kV = 0.2;
faultModel.kX = 1.0;

optsCCT.t_hi = 2.0;
optsCCT.tol  = 1e-3;
optsCCT.f0   = 60;
optsCCT.delta_h = deltaGrid;

CCT = NaN(numCases,1);
Vcrit = NaN(numCases,1);

for k = 1:numCases
    outCCT = ias_compute_cct_smib(clusters(k), faultModel, optsCCT);
    CCT(k)   = outCCT.CCT;
    Vcrit(k) = outCCT.Vcrit;
end

T2 = table((1:numCases).', Vcrit, CCT, [clusters.M_CL].', ...
    'VariableNames', {'Case','Vcrit','CCT_sec','M_CL'});
disp(T2);
writetable(T2, 'ias_ccpac_cct.csv');

exportgraphics(fig3, 'fig3_power_angle.png', 'Resolution', 300);
exportgraphics(fig4, 'fig4_current.png', 'Resolution', 300);
exportgraphics(fig5, 'fig5_margin.png', 'Resolution', 300);

% take the plot defaults set above off again
props = {'AxesFontSize','TextFontSize','AxesLabelFontSizeMultiplier', ...
         'AxesTitleFontSizeMultiplier','LegendFontSize', ...
         'AxesTickLabelInterpreter','TextInterpreter','LegendInterpreter'};
for k = 1:numel(props)
    set(groot, ['default' props{k}], 'remove');
end

function darrow(x1, x2, y, col, lw)
% horizontal double arrow from x1 to x2 at height y, in data units
    plot([x1 x2], [y y], '-', 'Color', col, 'LineWidth', lw);
    plot(x1, y, '<', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 5);
    plot(x2, y, '>', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 5);
end
