function genFig2_seedcorr_exemplars(chans)
% ------------------------------------------------------------------------------------------------
% Locator panel for the Figure 2A trace inset: where on the grid the four exemplar contacts sit.
% Two variants, both on P19's native pial, left lateral, same berry conventions as
% genFig2_seedcorr.m so they drop straight into the same figure.
%
%   *_context.png   every contact small and grey, the four exemplars full size and coloured by
%                   their correlation to the seed (viridis, clim [0 1]) -- the same colour each one
%                   carries in the trace panel, so trace and location can be matched by eye alone
%   *_onmap.png     the full seed-correlation berry map, with the four exemplars drawn enlarged so
%                   they can be picked out of the gradient they belong to
%
% Contacts and distances to the seed:
%   TG36-TG37   seed        r = 1.000
%   TG28-TG36    7.60 mm    r = 0.517
%   TG35-TG36   10.46 mm    r = 0.300
%   TG1-TG2     48.73 mm    r = 0.031
%
% No text is baked into either panel -- labels and leader lines belong in Illustrator.
% Data: paper/data/corrsource.mat, paper/data/natcortex_P19.mat (see export_natcortex.py).
% Run WITH a display. Output: paper/figs/2_seedcorr/exemplars_P19_*.png (transparent)
% ------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(chans)
    chans = {'TG36-TG37','TG28-TG36','TG35-TG36','TG1-TG2'};
end
SUB  = 'P19'; CSUB = 18;
CLIM = [0 1];
RBG  = 1.5;      % background contacts
REX  = 3.2;      % exemplars -- deliberately larger than the 2.5 used in the berry map, so they read
RMAP = 2.5;      % as picked out rather than just differently coloured
GREY = [0.62 0.62 0.62];

set(0,'DefaultFigureVisible','on');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths  = smooth_setup('quiet');
outdir = fullfile(paths.figs,'2_seedcorr'); if ~exist(outdir,'dir'), mkdir(outdir); end

S = load(fullfile(paths.data,'corrsource.mat')); s = S.corrsource{CSUB};
N = load(fullfile(paths.data,sprintf('natcortex_%s.mat',SUB)));
cortex.pos = double(N.natcortex.pos); cortex.tri = double(N.natcortex.tri);
cortex.curv = double(N.natcortex.curv(:));

elec = double(s.elec.nativechanpos); lab = string(s.label(:));
R    = double(s.HFBcorr);
[tf, ci] = ismember(string(chans(:)), lab);
if ~all(tf), error('genFig2_seedcorr_exemplars:chan','missing: %s', strjoin(chans(~tf),', ')); end
r = R(ci(1),:).';

isL   = elec(ci(1),1) < 0;
keepv = (1:size(cortex.pos,1))' > double(N.n_lh); if isL, keepv = ~keepv; end
cortex = submesh(cortex, keepv);
keepe = (elec(:,1) < 0) == isL;
elec = elec(keepe,:); r = r(keepe); lab = lab(keepe);
[~, ci] = ismember(string(chans(:)), lab);
vw = [90*(1-2*isL) 0];

idx = knnsearch(cortex.pos, elec); e = cortex.pos(idx,:);
cmap = viridis; nc = size(cmap,1);
cix  = min(max(round((r(ci) - CLIM(1))/(CLIM(2)-CLIM(1)) * (nc-1)) + 1, 1), nc);
excol = cmap(cix,:);

d = sqrt(sum((elec - elec(ci(1),:)).^2, 2));
fprintf('\n%s exemplars | %s, %d contacts\n', SUB, hemistr(isL), numel(lab));
fprintf('  %-12s %12s %10s   colour\n','channel','d to seed','r to seed');
for k = 1:numel(chans)
    fprintf('  %-12s %9.2f mm %10.3f   [%.2f %.2f %.2f]\n', chans{k}, d(ci(k)), r(ci(k)), excol(k,:));
end

% --- variant 1: grey context, exemplars coloured ---
newfig(cortex, vw);
other = true(size(r)); other(ci) = false;
ft_plot_cloud(e(other,:), zeros(nnz(other),1), 'cloudtype','surf', 'radius',RBG, ...
              'scalerad','no', 'colormap',repmat(GREY,256,1), 'clim',[-1 1]);
for k = 1:numel(chans)     % one call each: 'surf' bakes FaceColor, so per-contact colours survive
    ft_plot_cloud(e(ci(k),:), 0, 'cloudtype','surf', 'radius',REX, ...
                  'scalerad','no', 'colormap',repmat(excol(k,:),256,1), 'clim',[-1 1]);
end
finish(vw); savepanel(outdir, sprintf('exemplars_%s_context', SUB));

% --- variant 2: on the full seed-correlation map ---
newfig(cortex, vw);
ft_plot_cloud(e(other,:), r(other), 'cloudtype','surf', 'radius',RMAP, ...
              'scalerad','no', 'colormap',cmap, 'clim',CLIM);
for k = 1:numel(chans)
    ft_plot_cloud(e(ci(k),:), 0, 'cloudtype','surf', 'radius',REX, ...
                  'scalerad','no', 'colormap',repmat(excol(k,:),256,1), 'clim',[-1 1]);
end
finish(vw); savepanel(outdir, sprintf('exemplars_%s_onmap', SUB));

fprintf('wrote exemplar panels to %s\n', outdir);
end

% ==================================== helpers ====================================
function newfig(c, ~)
figure('color','w','position',[100 100 520 470]); axes('Position',[0.02 0.04 0.96 0.9]);
ft_plot_mesh(c,'vertexcolor','curv','edgecolor','none'); hold on;
end

function finish(vw)
view(vw); axis off vis3d; camzoom(0.90);   % spheres sit proud of the surface; zoom out so the
end                                        % cloud is not clipped at the frame

function s = hemistr(isL), if isL, s = 'lh'; else, s = 'rh'; end, end

function sub = submesh(c, keep)
vmap = zeros(size(c.pos,1),1); vmap(keep) = 1:nnz(keep);
fk = all(keep(c.tri),2);
sub.pos = c.pos(keep,:); sub.tri = vmap(c.tri(fk,:)); sub.curv = c.curv(keep);
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
if ~isempty(rr), rw = min(rr):max(rr); cw = min(cc):max(cc); C = C(rw,cw,:); alpha = alpha(rw,cw); end
imwrite(uint8(C), fullfile(outdir,[name '.png']), 'Alpha', alpha); delete(fw); delete(fkp);
end
