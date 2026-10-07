function x = variogram(v,vi,vf,D)

% step = (vf-vi)/20;
% 
% % compute variogram given distmance mat and dv
% x = nan(length(vi:step:vf),1);
% count = 1;
% for d = vi:step:vf
%     % graph
%     A = D<d & D>=(d-step);
%     A = zeros(size(A)) + triu(A);  % make upper triangular
%     A(logical(eye(size(A)))) = 0;  % remove recurrent connections
% 
%     % difference matrix
%     v_col = v*ones(1,length(v));   % compute difference between all obs
%     V = v_col - v_col';
% 
%     % compute variance
%     V(~logical(A)) = 0;
%     S = sparse(V);
%     S = S.^2;
%     semivar = full(nansum(S(:)))/(2*nansum(A(:)));
%     x(count) = semivar;
%     count = count + 1;
% end

nBins = 20;
step = (vf - vi) / nBins;
x = nan(nBins + 1, 1);

% Loop over distance bins
for i = 1:(nBins + 1)
    d = vi + (i - 1) * step;

    % Create mask for distances in the current bin
    A = (D < d) & (D >= (d - step));
    A = triu(A, 1);

    % Compute squared differences
    V = v(:) - v(:)';   % full pairwise difference matrix
    V(~A) = 0;          % mask by adjacency matrix
    V = sparse(V).^2;

    % Compute semivariance
    semivar = full(nansum(V(:))) / (2 * nansum(A(:)));
    x(i) = semivar;
end