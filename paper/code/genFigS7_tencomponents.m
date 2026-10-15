function genFigS7_tencomponents(K)
% ------------------------------------------------------------------------------------------------
% Supplementary Figure 7, ten-component comparison: how much of each dataset's observed group t-map
% the first K (default 10) principal components of each null ensemble capture, SMOOTH vs RSF, for all
% three datasets. Captured variance is |B_K'y|^2 / |y|^2 for the centred t-map y and the ensemble's
% first K PCs B_K (an orthogonal projection, see genFig5 pc_reconstruction); r is the correlation
% between y and its reconstruction, as plotted in the reconstruction curves.
%
% Reads paper/data/genFig5_results.mat (fields recon_var, recon_r, recon_k); runs no analysis.
% genFig5 calls this at the end of a full three-dataset run. Writes paper/figs/S7/tencomponents.pdf.
% ------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(K), K = 10; end
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..'));
paths = smooth_setup('quiet');
RES   = load(fullfile(paths.data, 'genFig5_results.mat'));
names = {'primary', 'berezutskaya', 'cogitate'};
labs  = {'Communication', 'Audiovisual language', 'Visual perception'};
V = nan(3, 2); r = nan(3, 2);                      % columns: SMOOTH, RSF (pc_reconstruction order)
for i = 1:3
    R = RES.(names{i}); k = find(R.recon_k == K, 1);
    V(i,:) = 100 * R.recon_var(k,:); r(i,:) = R.recon_r(k,:);
end
col = [0.20 0.35 0.65; 0.85 0.35 0.15];           % SMOOTH, RSF (as in the reconstruction curves)
f = figure('color','w','units','centimeters','position',[0 0 14 9]);
b = bar(V, 'grouped', 'BarWidth', 0.8, 'EdgeColor', 'none');
for e = 1:2, b(e).FaceColor = col(e,:); end
hold on
for e = 1:2
    x = b(e).XEndPoints;
    for i = 1:3
        text(x(i), V(i,e) + 2, sprintf('%.0f%%', V(i,e)), 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontName', 'Arial', 'FontSize', 9, 'Interpreter', 'none');
    end
end
set(gca, 'XTickLabel', labs, 'FontName', 'Arial', 'FontSize', 10, 'TickLabelInterpreter', 'none', 'Box', 'off');
ylim([0 100]); ylabel('Observed t-map variance captured (%)', 'Interpreter', 'none');
title(sprintf('First %d principal components of each null ensemble', K), 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'none');
legend({'SMOOTH', 'RSF'}, 'Location', 'northwest', 'Box', 'off', 'Interpreter', 'none');
grid on; set(gca, 'YGrid', 'on', 'XGrid', 'off');
outdir = fullfile(paths.figs, 'S7'); if ~isfolder(outdir), mkdir(outdir); end
exportgraphics(f, fullfile(outdir, 'tencomponents.pdf'), 'ContentType', 'vector'); close(f);
fprintf('\nfirst %d PCs, variance of the observed t-map captured (r):\n', K);
for i = 1:3
    fprintf('  %-22s SMOOTH %5.1f%% (r = %.3f) | RSF %5.1f%% (r = %.3f)\n', labs{i}, V(i,1), r(i,1), V(i,2), r(i,2));
end
fprintf('wrote %s\n', fullfile(outdir, 'tencomponents.pdf'));
end
