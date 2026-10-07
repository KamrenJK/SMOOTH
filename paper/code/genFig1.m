function genFig1
% ------------------------------------------------------------------------------------------------
% Figure 1. SMOOTH enables continuous group-level inference from sparsely sampled intracranial data.
% Real-data workflow schematic (socialcog rTPJ): every panel needed for the numbered inset (steps
% 1-7) plus the all-electrode-by-subject/DK/Schaefer scatter panels. Teal mono-diverging colour
% scheme, transparent PNG brain maps (house convention, see genFig2_seedcorr.m's savepanel).
% See genFig1_inputs.m for the companion toy-DV-icon panel.
%
% Data: paper/data/fig1_data.mat (export_fig1.py -- depends on a separate simulation pipeline that is
% not distributed; only the cached .mat is needed to regenerate this figure).
% Run WITH a display. Panels -> paper/figs/1/*
% ------------------------------------------------------------------------------------------------
set(0,'DefaultFigureVisible','on');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
D = load(fullfile(paths.data,'fig1_data.mat'));
outdir = fullfile(paths.figs,'1'); if ~exist(outdir,'dir'), mkdir(outdir); end
set(groot,'defaultTextInterpreter','none','defaultAxesFontName','Arial','defaultTextFontName','Arial');

c.pos = double(D.pial_pos); c.tri = double(D.pial_tri); c.curv = double(D.pial_curv(:));
VW = [90 0];                                                     % rh lateral (rTPJ)
TEAL=[14 124 134]/255; VIOLET=[106 79 163]/255; ROSE=[195 60 124]/255;
tealm = mono_div(TEAL,256); violetm = mono_div(VIOLET,256); rosem = mono_div(ROSE,256);
violetb = brighten(violetm,-0.7);                                % eigenmode berries only: darker/more
EIGCL = [-0.75 0.75];                                            % saturated ramp + narrowed clim, so
                                                                 % contacts read against the cortex
                                                                 % (white centre is a fixed point of
                                                                 % brighten, so 0 stays white)
BVL = [-4 4]; SVL = [-1 1]; TCL = [-3 3];                        % coverage / smoothed / t-map clims
ids = double(D.sub_ids);
for i = 1:double(D.n_sub)
    berry(c, double(D.sub_elec{i}), double(D.sub_dv{i}), BVL, tealm, VW); savepanel(outdir, sprintf('cov_s%d', ids(i)));
    overlay(c, double(D.sub_emp(:,i)), SVL, tealm, VW); savepanel(outdir, sprintf('smooth_s%d', ids(i)));
    for k = 1:size(D.sub_surr,3)                                 % 3 sign-flip surrogate maps per subject (rose)
        overlay(c, squeeze(double(D.sub_surr(i,:,k)))', SVL, rosem, VW); savepanel(outdir, sprintf('surr_s%d_%d', ids(i), k));
    end
end
overlay(c, double(D.tmap),  TCL, tealm, VW); savepanel(outdir,'tmap');
overlay(c, double(D.sig_t), TCL, tealm, VW); savepanel(outdir,'cluster');
cbar(outdir,'cbar_tmap', tealm, TCL, 't-value');
nulldist(double(D.posdist), double(D.tpj_mass), double(D.tpj_p), ROSE, TEAL, fullfile(outdir,'nulldist.pdf'));

% ---- graph + eigenmodes + spectra for the first two exemplars ----
midx = double(D.mode_idx);
for ii = [1 2]
    id = ids(ii);
    graphviz(c, double(D.sub_elec{ii}), double(D.graph_edges{ii}), VW); savepanel(outdir, sprintf('graph_s%d', id));
    E = double(D.eig{ii});
    for k = 1:numel(midx)
        berry(c, double(D.sub_elec{ii}), E(:,k), EIGCL, violetb, VW); savepanel(outdir, sprintf('eig_s%d_m%d', id, midx(k)));
    end
    mc = double(D.mc{ii}); mcs = double(D.mc_surr{ii}); n = numel(mc); M = max(abs(mc))*1.06;
    coefstem1(1:n, mc,  TEAL, 'mode index (k)', 'coefficient', M, fullfile(outdir, sprintf('spec_emp_s%d.pdf', id)));
    coefstem1(1:n, mcs, ROSE, 'mode index (k)', 'coefficient', M, fullfile(outdir, sprintf('spec_surr_s%d.pdf', id)));
end

% ---- all electrodes: by subject | by DK | by Schaefer-300 ----
ae = snap(c, double(D.all_elec)); sub = double(D.all_sub); Ns = max(sub)+1;
hh = mod((0:Ns-1)'*0.618, 1); subcol = hsv2rgb([hh, 0.62*ones(Ns,1), 0.88*ones(Ns,1)]);   % golden-angle distinct hues
elecscatter(c, ae, subcol(sub+1,:),        VW); savepanel(outdir,'allelec_subject');
elecscatter(c, ae, double(D.all_dk_rgb),   VW); savepanel(outdir,'allelec_dk');
elecscatter(c, ae, double(D.all_schaefer_rgb), VW); savepanel(outdir,'allelec_schaefer');
fprintf('wrote Figure 1 panels to %s (exemplars %s, rTPJ p=%.3f)\n', outdir, mat2str(ids(:)'), double(D.tpj_p));
end
% ==================================== helpers ====================================
function overlay(c, field, cl, cmap, vw)
field = field(:); figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(c,'vertexcolor','curv','edgecolor','none'); hold on;
k = isfinite(field);
if any(k), vm=zeros(size(c.pos,1),1); vm(k)=1:nnz(k); fk=all(k(c.tri),2);
  s.pos=c.pos(k,:); s.tri=vm(c.tri(fk,:)); ft_plot_mesh(s,'vertexcolor',field(k),'edgecolor','none'); end
colormap(cmap); clim(cl); view(vw); axis off;
end
function berry(c, elec, dv, cl, cmap, vw)
idx = knnsearch(c.pos, elec); e = c.pos(idx,:);
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(c,'vertexcolor','curv','edgecolor','none'); hold on;
ft_plot_cloud(e, dv(:), 'cloudtype','surf','radius',2.5,'scalerad','no','colormap',cmap,'clim',cl);
view(vw); axis off;
end
function savepanel(outdir,name)
fw=fullfile(outdir,[name '_w.png']); fk=fullfile(outdir,[name '_k.png']);
set(gcf,'InvertHardcopy','off','PaperPositionMode','auto'); set(findall(gcf,'type','axes'),'Color','none');
set(gcf,'Color','w'); print(gcf,fw,'-dpng','-r300'); set(gcf,'Color','k'); print(gcf,fk,'-dpng','-r300'); close(gcf);
W=double(imread(fw)); K=double(imread(fk)); alpha=min(max(1-mean(W-K,3)/255,0),1); C=min(max(K./max(alpha,1e-3),0),255);
[rr,cc]=find(alpha>0.02); if ~isempty(rr), r=min(rr):max(rr); c2=min(cc):max(cc); C=C(r,c2,:); alpha=alpha(r,c2); end
imwrite(uint8(C),fullfile(outdir,[name '.png']),'Alpha',alpha); delete(fw); delete(fk);
end
function nulldist(sm, obs, p, crose, cteal, fn)
sm=sm(:);
f = figure('color','w','position',[100 100 560 370]); hold on; box off
edges = linspace(0, max([sm;obs])*1.05, 40);
histogram(sm, edges, 'FaceColor',crose, 'EdgeColor','none', 'FaceAlpha',0.8);
xline(prctile(sm,95),'--','95th percentile','Color',crose,'LabelVerticalAlignment','bottom');
xline(obs,'-',sprintf('observed cluster  p=%.3f',p),'Color',cteal,'LineWidth',2,'LabelOrientation','horizontal');
set(gca,'FontSize',11,'FontName','Arial','TickDir','out'); xlabel('max cluster mass'); ylabel('count');
exportgraphics(f, fn, 'ContentType','vector'); close(f);
end
function graphviz(c, elec, edges, vw)
idx = knnsearch(c.pos, elec); e = c.pos(idx,:); margin = 6;
if mean(c.pos(:,1)) < 0, delta = min(c.pos(:,1)) - max(e(:,1)) - margin; else, delta = max(c.pos(:,1)) - min(e(:,1)) + margin; end
ev = e; ev(:,1) = e(:,1) + delta;                                          % rigid medial-lateral shift, in front
edges = edges(edges(:,3) >= 0.1, :); ecol = [0.4 0.4 0.45];
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(c,'vertexcolor','curv','edgecolor','none'); hold on;
w = edges(:,3); lw = 0.25 + 2.0*(w - min(w))/(max(w) - min(w) + 1e-9);
for m = 1:size(edges,1), ij = edges(m,1:2); plot3(ev(ij,1), ev(ij,2), ev(ij,3), '-', 'Color',ecol, 'LineWidth',lw(m)); end
scatter3(ev(:,1), ev(:,2), ev(:,3), 18, ecol, 'filled', 'MarkerEdgeColor','none'); view(vw); axis off;
end
function coefstem1(x, y, col, xl, yl, M, fn)
f = figure('color','w','Units','centimeters','Position',[3 3 3.7 3.6]);
ax = axes('Parent',f,'Units','centimeters','Position',[0.8 0.85 2.6 2.354]); hold(ax,'on'); box(ax,'off');
stem(ax, x, y, 'filled','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',1.2,'LineWidth',0.25);
yline(ax,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.25);
set(ax,'FontSize',8,'FontName','Arial','TickDir','out','XTick',[],'YTick',[]); xlim(ax,[x(1)-0.5 x(end)+0.5]); ylim(ax,[-M M]);
xlabel(ax,xl,'FontSize',9); ylabel(ax,yl,'FontSize',9); exportgraphics(f,fn,'ContentType','vector'); close(f);
end
function e = snap(c, elec), idx = knnsearch(c.pos, elec); e = c.pos(idx,:); end
function elecscatter(c, e, col, vw)
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(c,'vertexcolor','curv','edgecolor','none'); hold on;
scatter3(e(:,1),e(:,2),e(:,3), 14, col, 'filled', 'MarkerEdgeColor','none'); view(vw); axis off;
end
function cbar(outdir,name,cmap,cl,label)
f=figure('color','w','position',[100 100 150 430]);
ax=axes('Parent',f,'Position',[0.05 0.06 0.02 0.88]); axis(ax,'off'); colormap(ax,cmap); set(ax,'CLim',cl);
cb=colorbar(ax,'Position',[0.32 0.06 0.22 0.88]); cb.Limits=cl; cb.Ticks=[cl(1) mean(cl) cl(2)]; set(cb,'FontSize',14); ylabel(cb,label,'FontSize',14);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector'); close(f);
end
function cm = mono_div(hi, n)
lo=[0.29 0.29 0.29]; mid=[1 1 1]; x=linspace(0,1,n)';
cm=[interp1([0 .5 1],[lo(1) mid(1) hi(1)],x) interp1([0 .5 1],[lo(2) mid(2) hi(2)],x) interp1([0 .5 1],[lo(3) mid(3) hi(3)],x)];
end
