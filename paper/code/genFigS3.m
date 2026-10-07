function genFigS3
% ------------------------------------------------------------------------------------------------
% Supplementary Figure 3. False-positive control under anisotropic and nonstationary noise.
%
% Layout (panels exported separately for Illustrator assembly):
%   row 1  exemplar anisotropic fields,  rho = 1, 2, 3, 5, 10   (rho=1 is isotropic)
%   row 2  A-P mixing weight + exemplar nonstationary fields, q = 1, 2, 4, 6, 8  (q=1 is stationary)
%   row 3  FWER vs rho  and  FWER vs q, SMOOTH vs naive shuffle, y in [0 1], 95% binomial CI
% Exemplar fields are illustrative single draws at a fixed spectral slope (so the panels differ only
% in anisotropy/nonstationarity) -- NOT the realizations behind the FWER curves; the caption states
% this. The anisotropy row uses beta = 8 on the 1500-mode anisotropic bases, verified to reproduce the
% simulated fields' own spatial autocorrelation at every rho (within ~0.02 electrode-level Moran's I).
% The nonstationarity row uses beta = 29 on the 5000-mode isotropic basis -- the simulations' actual
% condition -- confirmed to match their recorded autocorrelation to ~0.01 across q (at beta = 8 / 1500
% modes this row would be measurably too rough: 15.6mm vs. the simulations' 18.8mm correlation
% half-decay).
%
% FWER panels are reported at one-tailed p <= 0.05, matching the main FWER figure's convention.
%
% Nugget condition: a = 1.0 (no nugget, the adversarial pure-1/f case).
% Data: paper/data/anisononstat_fig_data_v3.mat (produced by a separate simulation pipeline, not distributed).
% Run WITH a display. Panels -> paper/figs/S3/*
% ------------------------------------------------------------------------------------------------
set(0,'DefaultFigureVisible','on');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths  = smooth_setup('quiet');
D      = load(fullfile(paths.data,'anisononstat_fig_data_v3.mat'));
outdir = fullfile(paths.figs,'S3'); if ~exist(outdir,'dir'), mkdir(outdir); end

gi = get(groot,{'defaultTextInterpreter'});
set(groot,'defaultTextInterpreter','none','defaultAxesFontName','Arial','defaultTextFontName','Arial');
restore = onCleanup(@() set(groot,'defaultTextInterpreter',gi{1})); %#ok<NASGU>

cortex.pos = double(D.lh_pos); cortex.tri = double(D.lh_tri); cortex.curv = double(D.lh_curv(:));
try rdbu = flipud(slanCM('RdBu')); catch, rdbu = jet(256); end
FCL  = [-2.5 2.5];              % exemplar field colour limits (z-scored fields)
VW   = [-90 0];                 % lh lateral
COL  = struct('smooth',[0.082 0.396 0.753], 'naive',[0.776 0.157 0.157]);
thr  = double(D.thresh);

% ================= exemplar fields =================
for k = 1:numel(D.aniso_field_rho)
    overlay(cortex, double(D.aniso_fields(:,k)), FCL, rdbu, VW);
    savepanel(outdir, sprintf('field_aniso_rho%02d', D.aniso_field_rho(k)));
end
for k = 1:numel(D.nonstat_field_q)
    overlay(cortex, double(D.nonstat_fields(:,k)), FCL, rdbu, VW);
    savepanel(outdir, sprintf('field_nonstat_q%02d', D.nonstat_field_q(k)));
end
% A-P mixing weight (posterior short -> anterior long), its own scale + colorbar
overlay(cortex, double(D.ap_weight(:)), [0 1], parula(256), VW);
savepanel(outdir, 'field_ap_weight');

cbar(outdir,'cbar_field',    rdbu,        FCL,   'field amplitude (SD)');
cbar(outdir,'cbar_ap_weight',parula(256), [0 1], 'A-P weight w');

% ================= FWER panels =================
fwerpanel(outdir, 'fwer_aniso',   double(D.aniso_rho), ...
          double(D.aniso_smooth), double(D.aniso_smooth_lo), double(D.aniso_smooth_hi), ...
          double(D.aniso_naive),  double(D.aniso_naive_lo),  double(D.aniso_naive_hi), ...
          thr, COL, 'anisotropy ratio rho  (L_par / L_perp)');
fwerpanel(outdir, 'fwer_nonstat', double(D.nonstat_q), ...
          double(D.nonstat_smooth), double(D.nonstat_smooth_lo), double(D.nonstat_smooth_hi), ...
          double(D.nonstat_naive),  double(D.nonstat_naive_lo),  double(D.nonstat_naive_hi), ...
          thr, COL, 'nonstationarity ratio q  (L_long / L_short)');

fprintf('wrote panels to %s\n', outdir);
fprintf('  threshold: %s | nugget: %s | nsim = %d\n', char(D.thresh_label), char(D.nugget), D.nsim);
fprintf('  exemplar slopes: aniso beta=%g (%g modes) | nonstat beta=%g (%g modes)\n', ...
        D.aniso_beta, D.aniso_nmodes, D.nonstat_beta, D.nonstat_nmodes);
end

% ==================================== helpers ====================================
function fwerpanel(outdir,name,x,s,slo,shi,nv,nlo,nhi,thr,COL,xl)
f = figure('color','w','position',[100 100 470 430]); hold on; box off
xi = 1:numel(x);                                        % categorical spacing (grids are uneven)
errorbar(xi, nv, nv-nlo, nhi-nv, 'o-','LineWidth',2,'CapSize',4, ...
         'Color',COL.naive,'MarkerFaceColor',COL.naive,'MarkerEdgeColor','none','DisplayName','naive shuffle');
errorbar(xi, s,  s-slo,  shi-s,  'o-','LineWidth',2,'CapSize',4, ...
         'Color',COL.smooth,'MarkerFaceColor',COL.smooth,'MarkerEdgeColor','none','DisplayName','SMOOTH');
yline(thr,':','Color',[.45 .45 .45],'LineWidth',1.3,'HandleVisibility','off');
text(xi(end)+0.3, thr, sprintf('alpha = %.2f',thr), 'FontSize',9,'Color',[.45 .45 .45], ...
     'HorizontalAlignment','right','VerticalAlignment','bottom');   % right-aligned: at ylim [0 1] the
                                                                    % curves sit low and a left-hand
                                                                    % label collides with SMOOTH
set(gca,'XTick',xi,'XTickLabel',compose('%g',x),'FontSize',12);
xlim([xi(1)-0.4 xi(end)+0.4]); ylim([0 1]);
xlabel(xl); ylabel('FWER'); legend('Location','northeast','Box','off');
exportgraphics(f, fullfile(outdir,[name '.pdf']), 'ContentType','vector'); close(f);
end

function overlay(cortex, field, cl, cmap, vw)
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(cortex, 'vertexcolor', field(:), 'edgecolor','none');
colormap(cmap); clim(cl); view(vw); axis off;
end

function savepanel(outdir,name)
% Transparent-background raster PNG (house method; see genFigS6).
fw = fullfile(outdir,[name '_w.png']); fkp = fullfile(outdir,[name '_k.png']);
set(gcf,'InvertHardcopy','off','PaperPositionMode','auto'); set(findall(gcf,'type','axes'),'Color','none');
set(gcf,'Color','w'); print(gcf, fw, '-dpng','-r300');
set(gcf,'Color','k'); print(gcf, fkp,'-dpng','-r300'); close(gcf);
W = double(imread(fw)); K = double(imread(fkp));
alpha = min(max(1 - mean(W-K,3)/255, 0), 1);
C = min(max(K./max(alpha,1e-3), 0), 255);
[rr,cc] = find(alpha > 0.02);
if ~isempty(rr), r = min(rr):max(rr); c = min(cc):max(cc); C = C(r,c,:); alpha = alpha(r,c); end
imwrite(uint8(C), fullfile(outdir,[name '.png']), 'Alpha', alpha); delete(fw); delete(fkp);
end

function cbar(outdir,name,cmap,cl,label)
f = figure('color','w','position',[100 100 150 430]);
ax = axes('Parent',f,'Position',[0.05 0.06 0.02 0.88]); axis(ax,'off'); colormap(ax,cmap); set(ax,'CLim',cl);
cb = colorbar(ax,'Position',[0.32 0.06 0.22 0.88]); cb.Limits = cl; set(cb,'FontSize',13); ylabel(cb,label,'FontSize',13);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
