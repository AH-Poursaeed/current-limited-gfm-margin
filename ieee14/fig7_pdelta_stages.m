% Fig. 7: power-angle curves of the bus 2 cluster in the three stages with
% the current limit, and the post-fault curve without it, on one angle
% axis referred to the pre-fault Thevenin angle, with two zoomed insets.
% Needs STEP2_state.mat and ieee14_results.mat (steps 1 to 5).
% Writes fig7_pdelta_stages.png.

clear;

load('STEP2_state.mat',   'S');
load('ieee14_results.mat', 'rows');

set(groot,'defaultTextInterpreter',          'latex');
set(groot,'defaultLegendInterpreter',        'latex');
set(groot,'defaultAxesTickLabelInterpreter', 'latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
set(groot,'defaultAxesFontName',   'Times New Roman');
set(groot,'defaultTextFontName',   'Times New Roman');
set(groot,'defaultAxesFontSize',   11);
set(groot,'defaultTextFontSize',   11);
set(groot,'defaultLineLineWidth',  1.5);
set(groot,'defaultAxesLineWidth',  0.8);
set(groot,'defaultAxesBox',        'on');
set(groot,'defaultAxesXMinorTick', 'on');
set(groot,'defaultAxesYMinorTick', 'on');

C.pre   = [0.13 0.47 0.71];
C.fault = [0.84 0.15 0.16];
C.postL = [0.17 0.63 0.17];
C.postR = [0.58 0.40 0.74];
C.Pm    = [0.50 0.50 0.50];

LWIDTH_THIN  = 1.0;
LWIDTH_THICK = 1.8;

% each stage has its own Thevenin angle; the curves are shifted onto the
% pre-fault one
b = 2;  Ng = 4000;  Iref = 1e6;

r = rows(rows.bus==b,:);
assert(height(r)==1,'Row for bus %d not found.',b);

Eabs = S.Eabs_map(b);   Imax = S.Imax_map(b);   Pm = S.Pm_map(b);
Vpre = S.Vth_pre_map(b);    Zpre = S.Zth_pre_map(b);
Vf   = S.Vth_fault_map(b);  Zf   = S.Zth_fault_map(b);
Vpo  = S.Vth_post_map(b);   Zpo  = S.Zth_post_map(b);

theta_ref = angle(Vpre);
delta_h   = linspace(0, 2*pi, Ng).';

mapDelta = @(V) delta_h - (angle(V) - theta_ref);
Pe_pre   = pe_limited(mapDelta(Vpre), Vpre, Zpre, Eabs, S.Zf_gfm, Imax);
Pe_fault = pe_limited(mapDelta(Vf),   Vf,   Zf,   Eabs, S.Zf_gfm, Imax);
Pe_postL = pe_limited(mapDelta(Vpo),  Vpo,  Zpo,  Eabs, S.Zf_gfm, Imax);
Pe_postR = pe_limited(mapDelta(Vpo),  Vpo,  Zpo,  Eabs, S.Zf_gfm, Iref);

shift_post        = angle(Vpo) - theta_ref;
dSEP              = mod(r.sep_post               + shift_post, 2*pi);
dUEP_lim          = mod(r.uep_post_wrapped       + shift_post, 2*pi);
dUEP_ref          = mod(mod(r.sep_post + r.DeltaRef_rad*r.dir, 2*pi) + shift_post, 2*pi);

fig1 = figure('Color','w','Units','inches','Position',[1 1 8.6 5.2]);
ax1  = axes(fig1,'Position',[0.10 0.13 0.86 0.78]);
hold(ax1,'on');  grid(ax1,'on');

fill(ax1, [delta_h; flipud(delta_h)], [Pe_postR; flipud(Pe_postL)], ...
     C.postR,'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');

p1 = plot(ax1, delta_h, Pe_pre,   '-',  'Color',C.pre,   'LineWidth',LWIDTH_THICK, ...
          'DisplayName','$P_{e,h}^{\mathrm{pre}}(\delta_h)$');
p2 = plot(ax1, delta_h, Pe_fault, '-',  'Color',C.fault, 'LineWidth',LWIDTH_THICK, ...
          'DisplayName','$P_{e,h}^{\mathrm{fault}}(\delta_h)$');
p3 = plot(ax1, delta_h, Pe_postL, '-',  'Color',C.postL, 'LineWidth',LWIDTH_THICK, ...
          'DisplayName','$P_{e,h}^{\mathrm{post}}(\delta_h)$  (limited)');
p4 = plot(ax1, delta_h, Pe_postR, '--', 'Color',C.postR, 'LineWidth',LWIDTH_THIN,  ...
          'DisplayName','$P_{e,h}^{\mathrm{post,ref}}(\delta_h)$  (no limit)');

hl_Pm = yline(ax1, Pm, ':', 'Color',C.Pm, 'LineWidth',LWIDTH_THIN, ...
              'Label','$P_{m,h}$','LabelHorizontalAlignment','right', ...
              'Interpreter','latex','FontName','Times New Roman','FontSize',10);
hl_Pm.Annotation.LegendInformation.IconDisplayStyle = 'off';

hl_SEP  = xline(ax1, dSEP,     '-',  'Color',[0.2 0.2 0.2],'LineWidth',1.2, ...
                'Label','$\delta_h^{\mathrm{SEP,post}}$', ...
                'LabelOrientation','horizontal','LabelVerticalAlignment','top', ...
                'Interpreter','latex','FontName','Times New Roman','FontSize',9.5);
hl_UEPl = xline(ax1, dUEP_lim, '-',  'Color',C.postL,'LineWidth',1.2, ...
                'Label','$\delta_h^{\mathrm{UEP,post,lim}}$', ...
                'LabelOrientation','horizontal','LabelVerticalAlignment','top', ...
                'LabelHorizontalAlignment','left', ...
                'Interpreter','latex','FontName','Times New Roman','FontSize',9.5);
hl_UEPr = xline(ax1, dUEP_ref, '--', 'Color',C.postR,'LineWidth',1.0, ...
                'Label','$\delta_h^{\mathrm{UEP,post,ref}}$', ...
                'LabelOrientation','horizontal','LabelVerticalAlignment','top', ...
                'Interpreter','latex','FontName','Times New Roman','FontSize',9.5);
hl_SEP.Annotation.LegendInformation.IconDisplayStyle  = 'off';
hl_UEPl.Annotation.LegendInformation.IconDisplayStyle = 'off';
hl_UEPr.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlabel(ax1,'$\delta_h$ (rad),\ referenced to $\angle V_{th,h}^{\mathrm{pre}}$', ...
       'Interpreter','latex','FontName','Times New Roman','FontSize',11);
ylabel(ax1,'$P_{e,h}^{(s)}(\delta_h)$ (pu)', ...
       'Interpreter','latex','FontName','Times New Roman','FontSize',11);
title(ax1,sprintf('Stage-dependent $P$--$\\delta_h$ for cluster $h$ (bus %d): limited vs.\\ reference',b), ...
      'Interpreter','latex','FontName','Times New Roman','FontSize',12,'FontWeight','normal');

lg1 = legend(ax1,[p1 p2 p3 p4],'Location','southeast');
set(lg1,'Box','on','Interpreter','latex','FontName','Times New Roman','FontSize',10);
polish(ax1);
set(ax1, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi], ...
         'XTickLabel', {'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'}, ...
         'TickLabelInterpreter','latex');
% fixed axis limits; the insets below are placed for them
xlim(ax1, [0 7]);
ylim(ax1, [-6 6]);

% zoomed regions: x-range, title and position of the inset
zoomRegions = struct( ...
    'xLim',   {[0.00, 0.30],    [5.90, 2*pi]}, ...
    'label',  {'$\delta_h \in [0,\,0.3]$ rad',  ...
               '$\delta_h \in [5.9,\,2\pi]$ rad'}, ...
    'pos',    {[0.185 0.200 0.22 0.27],   ...
               [0.600 0.520 0.22 0.27]});

boxColors = {[0.2 0.2 0.6], [0.6 0.2 0.2]};

for iz = 1:2
    xZ = zoomRegions(iz).xLim;

    mask  = delta_h >= xZ(1) & delta_h <= xZ(2);
    allPe = [Pe_pre(mask); Pe_fault(mask); Pe_postL(mask); Pe_postR(mask)];
    ylo   = min(allPe) - 0.04*(max(allPe) - min(allPe));
    yhi   = max(allPe) + 0.04*(max(allPe) - min(allPe));

    fill(ax1, [xZ(1) xZ(2) xZ(2) xZ(1)], [ylo ylo yhi yhi], ...
         boxColors{iz}, 'FaceAlpha',0.06,'EdgeColor',boxColors{iz}, ...
         'LineWidth',1.1,'LineStyle','--','HandleVisibility','off');

    axI = axes(fig1,'Position', zoomRegions(iz).pos);
    hold(axI,'on');  grid(axI,'on');

    fill(axI,[delta_h; flipud(delta_h)],[Pe_postR; flipud(Pe_postL)], ...
         C.postR,'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');

    plot(axI, delta_h, Pe_pre,   '-',  'Color',C.pre,   'LineWidth',1.3);
    plot(axI, delta_h, Pe_fault, '-',  'Color',C.fault, 'LineWidth',1.3);
    plot(axI, delta_h, Pe_postL, '-',  'Color',C.postL, 'LineWidth',1.3);
    plot(axI, delta_h, Pe_postR, '--', 'Color',C.postR, 'LineWidth',0.9);

    yline(axI, Pm,       ':',  'Color',C.Pm,          'LineWidth',0.9, ...
          'HandleVisibility','off');
    xline(axI, dSEP,     '-',  'Color',[0.2 0.2 0.2], 'LineWidth',0.9, ...
          'HandleVisibility','off');
    xline(axI, dUEP_lim, '-',  'Color',C.postL,       'LineWidth',0.9, ...
          'HandleVisibility','off');
    xline(axI, dUEP_ref, '--', 'Color',C.postR,       'LineWidth',0.8, ...
          'HandleVisibility','off');

    xlim(axI, xZ);
    ylim(axI, [ylo yhi]);

    set(axI,'LineWidth',1.3,'XColor',boxColors{iz},'YColor',boxColors{iz});

    set(axI,'TickLabelInterpreter','latex','FontName','Times New Roman','FontSize',8.5, ...
            'TickDir','out','XMinorTick','on','YMinorTick','on', ...
            'Color',[0.98 0.98 0.98]);

    xlabel(axI, '$\delta_h$ (rad)', ...
           'Interpreter','latex','FontName','Times New Roman','FontSize',8.5);
    ylabel(axI, '$P_{e,h}$ (pu)', ...
           'Interpreter','latex','FontName','Times New Roman','FontSize',8.5);
    title(axI, zoomRegions(iz).label, ...
          'Interpreter','latex','FontName','Times New Roman','FontSize',9, ...
          'Color',boxColors{iz},'FontWeight','normal');

    % lines from the marked region to the inset
    axIPos = axI.Position;

    dataToFig = @(ax,xd,yd) [ ...
        ax.Position(1) + (xd - ax.XLim(1))/diff(ax.XLim)*ax.Position(3), ...
        ax.Position(2) + (yd - ax.YLim(1))/diff(ax.YLim)*ax.Position(4)];

    if iz == 1
        cMain = dataToFig(ax1, xZ(2), yhi);
        cInst = [axIPos(1), axIPos(2)+axIPos(4)];
        cMain2 = dataToFig(ax1, xZ(2), ylo);
        cInst2 = [axIPos(1), axIPos(2)];
    else
        cMain = dataToFig(ax1, xZ(1), yhi);
        cInst = [axIPos(1)+axIPos(3), axIPos(2)+axIPos(4)];
        cMain2 = dataToFig(ax1, xZ(1), ylo);
        cInst2 = [axIPos(1)+axIPos(3), axIPos(2)];
    end

    annotation(fig1,'line',[cMain(1) cInst(1)], [cMain(2) cInst(2)], ...
               'Color',boxColors{iz},'LineStyle','--','LineWidth',0.8);
    annotation(fig1,'line',[cMain2(1) cInst2(1)],[cMain2(2) cInst2(2)], ...
               'Color',boxColors{iz},'LineStyle','--','LineWidth',0.8);
end

drawnow;
exportgraphics(fig1, 'fig7_pdelta_stages.png', 'Resolution', 300);

% take the plot defaults set above off again
props = {'TextInterpreter','LegendInterpreter','AxesTickLabelInterpreter', ...
         'ColorbarTickLabelInterpreter','AxesFontName','TextFontName', ...
         'AxesFontSize','TextFontSize','LineLineWidth','AxesLineWidth', ...
         'AxesBox','AxesXMinorTick','AxesYMinorTick'};
for k = 1:numel(props)
    set(groot, ['default' props{k}], 'remove');
end

function polish(ax, fs)
    if nargin < 2, fs = 11; end
    set(ax, 'FontName','Times New Roman', 'FontSize',fs, ...
            'TickLabelInterpreter','latex', ...
            'LineWidth',0.8, 'TickDir','out', ...
            'XMinorTick','on','YMinorTick','on');
end
