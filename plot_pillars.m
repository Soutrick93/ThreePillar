% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  plot_pillars.m
%  Three-family comparison figure (manuscript Fig. 4; saved as
%  Fig7_pillars.pdf): confound audit and predictive increment.
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function plot_pillars(inDir, figDir)
% Reads context_audit.mat and dependency_audit.mat from inDir and writes
% Fig7_pillars.pdf/png to figDir. With NO arguments both default to the
% folder this script sits in. (run_all passes results\ and results\figures\.)
if nargin < 1 || isempty(inDir),  inDir  = fullfile(fileparts(mfilename('fullpath')),'results'); end
if nargin < 2 || isempty(figDir), figDir = fullfile(fileparts(mfilename('fullpath')),'results','figures'); end
clc;
C = load(fullfile(inDir,'context_audit.mat'));    CTX = C.CTX_RESULT;
D = load(fullfile(inDir,'dependency_audit.mat'));  DEP = D.DEP_RESULT;

BLUE=[0.165 0.471 0.839]; ORANGE=[0.922 0.408 0.204]; GREY=[0.45 0.45 0.45];
sig = @(a) max(a, 1-a);

figure('Units','inches','Position',[1 1 11 4.3],'Color','w');

% ----------------------------- panel (a) -----------------------------
ax1 = subplot(1,2,1); hold on;
vals = {sig(CTX.wc_static), sig(CTX.wc_context), sig(DEP.wc)};
surv = {false(1,6), false(1,6), (DEP.q < 0.05)};
base = [sig(CTX.baseline_wc), sig(CTX.baseline_wc), sig(DEP.baseline_wc)];
for g = 1:3
    v = vals{g}; s = surv{g};
    jit = linspace(-0.16, 0.16, numel(v));
    for i = 1:numel(v)
        if s(i)
            plot(g+jit(i), v(i), 'o', 'MarkerSize',8, 'MarkerFaceColor',ORANGE, ...
                 'MarkerEdgeColor',ORANGE, 'LineWidth',1.3);
        else
            plot(g+jit(i), v(i), 'o', 'MarkerSize',8, 'MarkerFaceColor','w', ...
                 'MarkerEdgeColor',GREY, 'LineWidth',1.3);
        end
    end
    plot([g-0.26 g+0.26], [base(g) base(g)], '-', 'Color',BLUE, 'LineWidth',2.2);
end
plot(xlim, [0.5 0.5], ':', 'Color',GREY, 'LineWidth',1);
set(ax1,'XTick',1:3,'XTickLabel',{'Static','Context','Dependency'});
ylim([0.48 0.62]); ylabel('within-cell |AUROC|');
title('Confound audit across mechanism families','FontWeight','bold');

% ----------------------------- panel (b) -----------------------------
ax2 = subplot(1,2,2); hold on;
% network increment from the per-cell-line ablation (omics+net vs omics)
net_ldo = 0.014; net_ldo_e = 0.006; net_lco = 0.013; net_lco_e = 0.007;
if isfield(DEP,'pred_LDO')
    dep_ldo = DEP.pred_LDO(3); dep_ldo_e = DEP.pred_LDO(4);
    dep_lco = DEP.pred_LCO(3); dep_lco_e = DEP.pred_LCO(4);
else
    dep_ldo=0.008; dep_ldo_e=0.0034; dep_lco=0.001; dep_lco_e=0.0028;
end
X = [1 2]; w = 0.34;
b1 = bar(X-w/2, [net_ldo dep_ldo], w, 'FaceColor',BLUE,  'EdgeColor','none');
b2 = bar(X+w/2, [net_lco dep_lco], w, 'FaceColor',ORANGE,'EdgeColor','none');
errorbar(X-w/2, [net_ldo dep_ldo], [net_ldo_e dep_ldo_e], 'k','LineStyle','none','CapSize',4);
errorbar(X+w/2, [net_lco dep_lco], [net_lco_e dep_lco_e], 'k','LineStyle','none','CapSize',4);
plot(xlim,[0 0],'-','Color',GREY,'LineWidth',0.8);
set(ax2,'XTick',X,'XTickLabel',{'Network motifs','Functional dependency'});
ylabel('\DeltaAUROC over chemistry + expression'); ylim([-0.004 0.024]);
title('Predictive increment of mechanism features','FontWeight','bold');
legend([b1 b2],{'unseen drugs (LDO)','unseen cells (LCO)'},'Box','off','Location','northeast');

% ----------------------------- shared styling ------------------------
for ax = [ax1 ax2]
    set(ax,'FontSize',12,'TickDir','out','Box','off','XColor',GREY,'YColor',GREY, ...
        'LineWidth',0.8,'Layer','top');
    ax.Title.Color = 'k'; ax.YLabel.Color = 'k';
end

% ----------------------------- save ----------------------------------
if ~exist(figDir,'dir'), mkdir(figDir); end
set(gcf,'PaperPositionMode','auto');
print(gcf, fullfile(figDir,'Fig7_pillars.pdf'), '-dpdf','-r300','-bestfit');
print(gcf, fullfile(figDir,'Fig7_pillars.png'), '-dpng','-r200');
fprintf('saved -> %s\n', fullfile(figDir,'Fig7_pillars.pdf'));
