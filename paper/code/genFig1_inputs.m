function genFig1_inputs
% ------------------------------------------------------------------------------------------------
% Figure 1 inset: toy illustrations of candidate DV inputs to SMOOTH ("many dependent variables"
% panel) -- scatter + linear fit | temporal-generalization matrix (8x8, autocorrelated) | PAC
% phase-amplitude histogram (sinusoidal) | beta bursts ("burst rate"). Clean vector icons, teal.
% Companion to genFig1.m; illustrative only, not derived from any dataset.
% Panels -> paper/figs/1_inputs/*
% ------------------------------------------------------------------------------------------------
set(0,'DefaultFigureVisible','on');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
set(groot,'defaultTextInterpreter','none','defaultAxesFontName','Arial','defaultTextFontName','Arial');
od = fullfile(paths.figs,'1_inputs'); if ~exist(od,'dir'), mkdir(od); end
TEAL=[14 124 134]/255; GREY=[0.55 0.55 0.58];

% 1) scatter + linear fit (10 dots)
rng(1); x=linspace(0.05,0.95,10)'; y=1.05*x + 0.11*randn(10,1) + 0.12;
f=fig(3.6,3.2); ax=axescm(f,[0.75 0.7 2.6 2.3]);
p=polyfit(x,y,1); plot(ax,[0 1],polyval(p,[0 1]),'-','Color',GREY,'LineWidth',1.2); hold(ax,'on');
scatter(ax,x,y,28,TEAL,'filled');
box(ax,'off'); set(ax,'XTick',[],'YTick',[],'TickDir','out'); xlim(ax,[-0.03 1.03]); ylim(ax,[min(y)-0.12 max(y)+0.12]);
xlabel(ax,'x','FontSize',11); ylabel(ax,'y','FontSize',11); save_vec(f,od,'in_scatter');

% 2) temporal-generalization matrix (8x8, autocorrelated diagonal band)
[X,Y]=meshgrid(1:8,1:8); rng(2);
M = 0.5 + 0.5*exp(-((X-Y).^2)/(2*1.7^2)) + 0.05*randn(8);         % broad diagonal band + a little noise
f=fig(3.0,3.0); ax=axescm(f,[0.25 0.25 2.5 2.5]);
imagesc(ax,M); axis(ax,'square'); colormap(ax, seq(TEAL,256));
set(ax,'XTick',[],'YTick',[],'Box','on','XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3]); save_vec(f,od,'in_tgm');

% 3) PAC phase-amplitude histogram (sinusoidal modulation)
nb=18; phc=(0.5:nb-0.5)/nb*2*pi; amp=1 + 0.55*cos(phc - pi/3);
f=fig(3.6,3.0); ax=axescm(f,[0.75 0.7 2.6 2.1]);
bar(ax,phc,amp,1,'FaceColor',TEAL,'EdgeColor','w','LineWidth',0.3); hold(ax,'on');
pp=linspace(0,2*pi,200); plot(ax,pp, 1+0.55*cos(pp-pi/3),'-','Color',GREY,'LineWidth',1.2);
box(ax,'off'); set(ax,'XTick',[],'YTick',[],'TickDir','out'); xlim(ax,[0 2*pi]); ylim(ax,[0 1.7]);
xlabel(ax,'phase','FontSize',11); ylabel(ax,'amplitude','FontSize',11); save_vec(f,od,'in_pac');

% 4) beta bursts ("burst rate")
t=linspace(0,1,900); sig=zeros(size(t));
for cc=[0.22 0.55 0.82], sig = sig + exp(-((t-cc)/0.028).^2).*sin(2*pi*22*(t-cc)); end
f=fig(4.4,2.0); ax=axescm(f,[0.12 0.12 4.15 1.75]);
plot(ax,t,sig,'-','Color',TEAL,'LineWidth',0.9); axis(ax,'off'); ylim(ax,[-1.15 1.15]); save_vec(f,od,'in_burst');

fprintf('wrote input icons to %s\n', od);
end
function f=fig(w,h), f=figure('color','w','Units','centimeters','Position',[3 3 w h]); end
function ax=axescm(f,pos), ax=axes('Parent',f,'Units','centimeters','Position',pos); end
function save_vec(f,od,name), exportgraphics(f,fullfile(od,[name '.pdf']),'ContentType','vector'); close(f); end
function cm=seq(hi,n), w=[1 1 1]; x=linspace(0,1,n)'; cm=[interp1([0 1],[w(1) hi(1)],x) interp1([0 1],[w(2) hi(2)],x) interp1([0 1],[w(3) hi(3)],x)]; end
