function genFigS4
% ------------------------------------------------------------------------------------------------
% Supplementary Figure 4. Simulation of ground-truth effects for sensitivity analysis. Illustrates the
% detection-sensitivity simulation end to end: geodesic Gaussian bump + per-subject 1/f-GRF surface
% noise + field=bump+noise + electrode berries + group SMOOTH t-map/cluster + coverage + the
% sensitivity map, plus a panel showing the 1/f-GRF fit to the empirical spatial autocorrelation.
%
% The noise panels and the rho(d) fit panel use the 1/f model (a, beta). Data: paper/data/methods_fig_data.mat
% (produced by the same separate simulation pipeline as the 1/f sims, not distributed). Plotting is genFig-style:
% ft_plot_mesh(...,'vertexcolor','curv') backdrop + overlay + colormap + clim + view -- NO camlight.
% Run WITH a display. Panels -> paper/figs/S4/*.pdf (named by view for easy assembly).
% ------------------------------------------------------------------------------------------------
set(0,'DefaultFigureVisible','on');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths  = smooth_setup('quiet');
D      = load(fullfile(paths.data,'methods_fig_data.mat'));
outdir = fullfile(paths.figs,'S4'); if ~exist(outdir,'dir'), mkdir(outdir); end

cortex.pos = double(D.pial_pos); cortex.tri = double(D.pial_tri); cortex.curv = double(D.curv(:));
seedp = double(D.seed_pos); crop = double(D.geo_from_seed(:)) < double(D.crop_radius);
A0 = D.params.A0;
bwr = diverge_bwr(256); mag = magma(256); grn = seqmap([0.97 0.98 0.95],[0.0 0.30 0.13],256);

VIEWS = {'lateral',[-90 0]};
for vv = 1:size(VIEWS,1)
    vn = VIEWS{vv,1}; vw = VIEWS{vv,2};

    % ---- METHOD (bump cropped to 1.5*FWHM around the seed; noise/field over the full surface) ----
    allv = true(size(cortex.pos,1),1);
    overlay(cortex,double(D.bump_surf),crop,[-A0 A0],bwr,vw,seedp); title(sprintf('Gaussian bump (%s)',vn)); savepanel(outdir,['01_bump_' vn]);
    for s = 1:size(D.noise_surf_show,1)
        overlay(cortex,double(D.noise_surf_show(s,:))',allv,[-A0 A0],bwr,vw,[]);   title(sprintf('1/f GRF noise - subj %d (%s)',D.show_subj(s),vn)); savepanel(outdir,sprintf('02_noise_s%d_%s',s,vn));
        overlay(cortex,double(D.field_surf_show(s,:))',allv,[-A0 A0],bwr,vw,seedp); title(sprintf('field=bump+noise - subj %d (%s)',D.show_subj(s),vn)); savepanel(outdir,sprintf('03_field_s%d_%s',s,vn));
    end

    % ---- BERRIES ----
    for s = 1:numel(D.elec_show)
        berry(cortex,double(D.elec_show{s}),double(D.elec_show_dv{s}),[-A0 A0],bwr,vw); title(sprintf('berry subj %d (%s)',D.show_subj(s),vn)); savepanel(outdir,sprintf('04_berry_s%d_%s',s,vn));
    end
    berry(cortex,double(D.elec_all),double(D.elec_all_dv),[-A0 A0],bwr,vw); title(sprintf('all subjects - electrode DV (%s)',vn)); savepanel(outdir,['04b_berry_all_' vn]);
    berry(cortex,double(D.elec_all),0.5*ones(size(D.elec_all,1),1),[0 1],[0.85 0.55 0.15],vw); title(sprintf('all electrodes (%s)',vn)); savepanel(outdir,['05_coverage_berry_' vn]);

    % ---- DETECT ----
    tm = double(D.tmap);
    overlay(cortex,tm,isfinite(tm),[-A0 A0],bwr,vw,seedp); title(sprintf('group t-map (%s)',vn)); savepanel(outdir,['06_tmap_' vn]);
    overlay(cortex,tm,logical(D.cluster_mask),[-4 4],bwr,vw,seedp); title(sprintf('significant cluster - masked t (%s)',vn)); savepanel(outdir,['07_cluster_' vn]);

    % ---- COVERAGE (# subjects and # electrodes) ----
    cs = double(D.coverage);      cs(cs < D.params.minnbsub) = NaN;
    overlay(cortex,cs,isfinite(cs),[3 20],grn,vw,[]); title(sprintf('coverage: # subjects (%s)',vn)); savepanel(outdir,['08_coverage_subj_' vn]);
    ce = double(D.coverage_elec); ce(ce < 1) = NaN;
    overlay(cortex,ce,isfinite(ce),[0 50],grn,vw,[]); title(sprintf('coverage: # electrodes (%s)',vn)); savepanel(outdir,['09_coverage_elec_' vn]);

    % ---- SENSITIVITY (LH) ----
    overlay(cortex,double(D.sens_map),isfinite(double(D.sens_map)),[0 1],mag,vw,[]); title(sprintf('sensitivity LH (%s)',vn)); savepanel(outdir,['10_sensitivity_lh_' vn]);
end

% ---- SENSITIVITY (RH): results map on the other hemisphere (rest of the figure is LH) ----
rhc.pos = double(D.rh_pial_pos); rhc.tri = double(D.rh_pial_tri); rhc.curv = double(D.rh_curv(:));
overlay(rhc,double(D.rh_sens_map),isfinite(double(D.rh_sens_map)),[0 1],mag,[90 0],[]); title('sensitivity RH (lateral)'); savepanel(outdir,'10_sensitivity_rh_lateral');

% ---- standalone colorbars (per plot type; per-subject bump/noise/field/berry share one) ----
cbar(outdir,'cbar_persubject',    bwr, [-A0 A0], 'field  (noise SD)');
cbar(outdir,'cbar_tmap',          bwr, [-A0 A0], 'group t');
cbar(outdir,'cbar_cluster',       bwr, [-4 4],   'cluster t (masked)');
cbar(outdir,'cbar_coverage_subj', grn, [3 20],   '# subjects');
cbar(outdir,'cbar_coverage_elec', grn, [0 50],   '# electrodes');
cbar(outdir,'cbar_sensitivity',   mag, [0 1],    'sensitivity');

% ---- 2D curves ----
curve2d(outdir,'11_sens_vs_subjects',   double(D.cov_curve_x),double(D.cov_curve_y),double(D.cov_curve_ci),'coverage (# subjects)','detection sensitivity','sensitivity vs # subjects');
hist2d(outdir,'11b_sens_vs_subjects_hist',  double(D.cov_pts),double(D.sens_pts),double(D.cov_max),round(double(D.cov_max)), double(D.cov_curve_x),double(D.cov_curve_y),'coverage (# subjects)','detection sensitivity','sensitivity vs # subjects (all vertices)');
curve2d(outdir,'12_sens_vs_electrodes', double(D.ne_curve_x), double(D.ne_curve_y), double(D.ne_curve_ci), 'coverage (# electrodes)','detection sensitivity','sensitivity vs # electrodes');
hist2d(outdir,'12b_sens_vs_electrodes_hist', double(D.ne_pts), double(D.sens_pts), double(D.ne_max), 45, double(D.ne_curve_x), double(D.ne_curve_y), 'coverage (# electrodes)','detection sensitivity','sensitivity vs # electrodes (all vertices)');
bump1d(outdir, A0, double(D.params.FWHM));

% ---- 1/f-GRF autocorrelation fit (empirical rho(d) vs the fitted 1/f model) ----
ok = logical(D.rho_ok); f = figure('color','w','position',[100 100 480 380]);
plot(double(D.rho_cen(ok)),double(D.rho_emp(ok)),'o','MarkerSize',7,'MarkerFaceColor',[0.20 0.40 0.65],'MarkerEdgeColor','none'); hold on
plot(double(D.rho_dd),double(D.rho_model),'-','LineWidth',2,'Color',[0.85 0.35 0.13]); grid on; box off
xlabel('electrode separation d (mm)'); ylabel('correlation \rho(d)'); xlim([0 60]);
legend({'empirical (socialcog)',sprintf('1/f GRF model  a=%.2f, \\beta=%.1f',D.sa_fit.a,D.sa_fit.beta)},'Location','northeast');
title('1/f GRF spatial autocorrelation matches the data'); set(gca,'FontSize',12);
exportgraphics(f,fullfile(outdir,'13_grf_spectrum.pdf'),'ContentType','vector'); close(f);

fprintf('wrote panels to %s\n', outdir);
end

% ==================================== helpers (genFig-style: no camlight) ====================================
function overlay(cortex, field, showmask, cl, cmap, vw, markpos)
field = field(:);
figure('color','w','position',[100 100 560 500]); axes('Position',[0.02 0.06 0.96 0.9]);
ft_plot_mesh(cortex, 'vertexcolor','curv', 'edgecolor','none'); hold on;
keep = showmask(:) & isfinite(field);
if any(keep)
    sub = submesh(cortex, keep);
    ft_plot_mesh(sub, 'vertexcolor', field(keep), 'edgecolor','none');
end
colormap(cmap); clim(cl); view(vw);
if ~isempty(markpos)
    d = markpos - mean(cortex.pos,1); d = d/norm(d);        % outward from cortex centroid
    mp = markpos + 12*d;                                    % float star ~12mm in front so it isn't occluded
    % plot3(mp(1),mp(2),mp(3),'p','MarkerSize',15,'MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor','k','LineWidth',0.7);
end
end

function berry(cortex, elec, dv, cl, cmap, vw)
dv = dv(:);
figure('color','w','position',[100 100 560 500]); axes('Position',[0.02 0.06 0.96 0.9]);
ft_plot_mesh(cortex, 'vertexcolor','curv', 'edgecolor','none'); hold on;
ft_plot_cloud(elec, dv, 'cloudtype','surf','radius',2.5,'scalerad','no','colormap',cmap,'clim',cl);
view(vw);
end

function sub = submesh(cortex, keep)
vmap = zeros(size(cortex.pos,1),1); vmap(keep) = 1:nnz(keep);
fk = all(keep(cortex.tri),2);
sub.pos = cortex.pos(keep,:); sub.tri = vmap(cortex.tri(fk,:));
end

function savepanel(outdir,name)
exportgraphics(gcf, fullfile(outdir,[name '.pdf']), 'ContentType','image','Resolution',300); close(gcf);
end

function cbar(outdir,name,cmap,cl,label)
f = figure('color','w','position',[100 100 150 430]);
ax = axes('Parent',f,'Position',[0.05 0.06 0.02 0.88]); axis(ax,'off'); colormap(ax,cmap); set(ax,'CLim',cl);
cb = colorbar(ax,'Position',[0.32 0.06 0.22 0.88]); cb.Limits = cl; set(cb,'FontSize',14); ylabel(cb,label,'FontSize',14);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end

function curve2d(outdir,name,x,y,e,xl,yl,ttl)
x=x(:); y=y(:); e=e(:); ok=isfinite(x)&isfinite(y);
f = figure('color','w','position',[100 100 460 380]);
errorbar(x(ok),y(ok),e(ok),'-o','LineWidth',2,'Color',[0.20 0.40 0.65], ...
         'MarkerFaceColor',[0.20 0.40 0.65],'MarkerEdgeColor','none','CapSize',4); grid on; box off
xlabel(xl); ylabel(yl); ylim([0 1]); set(gca,'FontSize',12); title(ttl);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end

function bump1d(outdir,A0,FWHM)
sig = FWHM/(2*sqrt(2*log(2)));
x = linspace(-40,40,600); y = A0*exp(-x.^2/(2*sig^2)); hm = A0/2;
f = figure('color','w','position',[100 100 500 380]);
plot(x,y,'-','LineWidth',2.2,'Color',[0.75 0.15 0.15]); hold on; box off
plot([-FWHM/2 FWHM/2],[hm hm],'k-','LineWidth',1.4);
plot([-FWHM/2 -FWHM/2],[0 hm],'k:','LineWidth',1); plot([FWHM/2 FWHM/2],[0 hm],'k:','LineWidth',1);
plot([-FWHM/2 FWHM/2],[hm hm],'k.','MarkerSize',12);
text(0,hm+0.18,sprintf('FWHM = %g mm',FWHM),'HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
text(1.5,A0-0.05,sprintf('peak = A_0 = %g SD',A0),'FontSize',11,'VerticalAlignment','top');
xlabel('geodesic distance from center (mm)'); ylabel('bump amplitude (noise SD)');
ylim([0 A0*1.12]); xlim([-40 40]); set(gca,'FontSize',12); title('planted Gaussian bump (1D profile)');
exportgraphics(f,fullfile(outdir,'14_bump_fwhm.pdf'),'ContentType','vector'); close(f);
end

function hist2d(outdir,name,x,y,xmax,nx,curx,cury,xl,yl,ttl)
x=x(:); y=y(:); k=isfinite(x)&isfinite(y)&(x<=xmax);
f = figure('color','w','position',[100 100 500 400]);
h = histogram2(x(k),y(k),'XBinEdges',linspace(0,xmax,nx+1),'YBinEdges',linspace(0,1,41), ...
           'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none'); hold on
colormap(magma(256)); set(gca,'ColorScale','log'); clim([1 max(h.Values(:))]); cb=colorbar; ylabel(cb,'# vertices');
plot(curx,cury,'-','LineWidth',2,'Color',[1 1 1]);
plot(curx,cury,'o','MarkerSize',4,'MarkerFaceColor','w','MarkerEdgeColor','k');
xlabel(xl); ylabel(yl); xlim([0 xmax]); ylim([0 1]); set(gca,'FontSize',12); title(ttl);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end

function cm = diverge_bwr(n)
b=[0.23 0.30 0.75]; w=[1 1 1]; r=[0.75 0.15 0.15]; x=linspace(0,1,n)';
cm=[interp1([0 .5 1],[b(1) w(1) r(1)],x) interp1([0 .5 1],[b(2) w(2) r(2)],x) interp1([0 .5 1],[b(3) w(3) r(3)],x)];
end
function cm = seqmap(a,b,n)
x=linspace(0,1,n)'; cm=[interp1([0 1],[a(1) b(1)],x) interp1([0 1],[a(2) b(2)],x) interp1([0 1],[a(3) b(3)],x)];
end
