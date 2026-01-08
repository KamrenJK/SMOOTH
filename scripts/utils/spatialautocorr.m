function I = spatialautocorr(A,D)
% compute spatial autocorrelation (Moran's I)
% A: n x 1 spatial map
% D: n x n pairwise distance matrix

% clean out NaNs (should be done before calling function)
% idx      = isnan(A);
% A(idx)   = [];
% D(idx,:) = [];
% D(:,idx) = [];

D(D==0) = inf;
Dinv    = 1./D;
Dinv    = Dinv./max(Dinv,[],1);  % weight by inverse distance
A = A - mean(A,'all');  % center on zero
I = (numel(A)/sum(Dinv,'all'))*(A(:)'*Dinv*A(:))/(A(:)'*A(:));