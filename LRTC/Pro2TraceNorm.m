function [X, n, Sigma2] = Pro2TraceNorm(Z, tau, TraceNormToRankParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [S, V, D, Sigma2] = MySVDtau(Z, tau);
% V = max(diag(V) - tau, 0);
% n = sum(V > 0);
% X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

%TraceNormToRankParam : - Set to zero (default) for trace norm minimization (soft thresholding of singular values).
%                       - Set to inf for actual rank minimization (hard thresholding of singular values).
%                       - Other positive values should give intermediate results (keeping some high frequencies in the final matrix).

%% new
[m, n] = size(Z);

if(nargin<3)
    TraceNormToRankParam=0;
end

if 2*m < n
    AAT = Z*Z';
    [S, Sigma2, D] = svd(AAT);
    Sigma2 = diag(Sigma2);
    V = sqrt(Sigma2);
    %tol = max(size(Z)) * eps(max(V));
    %n = sum(V > max(tol, tau));
    n = sum(V > tau);
    
    if(isinf(TraceNormToRankParam))
        sub = zeros(1,n);
    else
        sub = ([0 [1:n-1]/(n-1)].^TraceNormToRankParam)*tau;
    end
    
    mid = max(V(1:n)-sub', 0) ./ V(1:n) ;
    X = S(:, 1:n) * diag(mid) * ( S(:, 1:n)' * Z );
    return;
end
if m > 2*n
    ATA = Z'*Z;
    [S, Sigma2, D] = svd(ATA);
    Sigma2 = diag(Sigma2);
    V = sqrt(Sigma2);
    %tol = max(size(Z)) * eps(max(V));
    %n = sum(V > max(tol, tau));
    n = sum(V > tau);
    
    if(isinf(TraceNormToRankParam))
        sub = zeros(1,n);
    else
        sub = ([0 [1:n-1]/(n-1)].^TraceNormToRankParam)*tau;
    end
    
    mid = max(V(1:n)-sub', 0) ./ V(1:n) ;
    X = (Z *  D(:, 1:n)) * bsxfun(@times, D(:, 1:n)', mid);
    %X = Z * ( D(:, 1:n) * diag(mid) * D(:, 1:n)' );
    return;
    
    %[X, n, Sigma2] = Pro2TraceNorm(Z', tau, TraceNormToRankParam);
    %X = X';
    %return;
end

%%{
[S,V,D] = svd(Z,'econ');
Sigma2 = diag(V).^2;
n = sum(diag(V) > tau);

if(isinf(TraceNormToRankParam))
    sub = zeros(1,n);
else
    sub=([0 [1:n-1]/(n-1)].^TraceNormToRankParam)*tau;
end
V(1:n,1:n) = diag(diag(V(1:n,1:n))-sub');
%V(2:n,2:n) = max(V(2:n,2:n)-subtractval, 0);
X = S(:, 1:n) * V(1:n,1:n) * D(:, 1:n)';
%}

