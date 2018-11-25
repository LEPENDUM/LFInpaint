%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low rank tensor/matrix completion of the paper : Light Field inpainting via low rank matrix completion.
% The code was modified from the tensor completion code of the HaLRTC method in the paper: 
% "Tensor completion for estimating missing values in visual data"
% The code may still be compatible for tensors but not tested. (minimizing the averaged ranks of
% matrices unfolded along the different tensors dimensions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% - T : input matrix to complete/approximate.
% - Omega : matrix describing the known (true) or unknown (false) entries.
% - alpha : (only for tensors) weights on the different tensor dimensions.
% - betaMult : parameter controling convergence speed vs accuraty tradeoff.
% It must be higher than 1. (using betaMult=1 may never converge).
% - maxIter : maximum number of iterations.
% - epsilon : global noise tolerance parameter. The parameter may be given as a vector of 2-elements : 
%            - 1st element = epsilon value used for stopping criterion /
%            - 2nd element = epsilon value used in the error constraint.
%            Note that the input epsilon values are normalized (multiplied by the norm of the input matrix taking only the known elements).
% - TraceNormToRankParam : 0 (default) => trace norm minimization (soft thresholding of singular values) /
%                          inf => rank minimization /
%                          Other positive values should give intermediate results.
% - CorruptCols : List of indices of the columns in the matrix T for which even the known entries may contain errors.
% - maxRank : maximum rank (default=inf). The algorithm stops when maxRank is reached.
% - X : User defined initialization of the approximation matrix. (default = T_ij if W_ij>0 / mean(T) otherwise).

% Outputs:
% - X_final : Completed matrix without any modification of the entries considered as known in the input.
% - X : Completed matrix with possible modifications of the entries considered as known in the input.
% - M : Cell containing matrix M (In case of tensors, the cell contains all the matrices unfolded along each tensor dimension).
%       A way to verify the convergence is to check that X is sufficiently close to M.
%       M has exactly the rank obtained at the last iteration, but may not satisfty the constraint on weighted norm if convergence is not reached.
% - rank : values of the matrix rank (for tensors: list of ranks of the matrices unfolded along each tensor dimension).
% - errList : list of error at each iteration (if the algrithm converged errList(end) should be less than epsilon).


function [X_final, X, M, rank, errList] = LRTCADMM(T, Omega, alpha, betaMult, maxIter, epsilon, TraceNormToRankParam, CorruptCols, maxRank, X)

DoCenter = false; %center the data before completion.
legacyThreshold = false; %true  : threshold value determined for nuclear norm minimization.
                         %false : threshold value determined for rank minimization.
                         %This parameter is automatically set to true when TraceNormToRankParam=0 (nuclear norm minimization).


Omega = logical(Omega);
meanT = sum(T(:).*Omega(:)) / sum(Omega(:));

if(~exist('CorruptCols','var'))
    if(isa(T,'single')), CorruptCols = zeros(0,'uint32');
    else, CorruptCols = [];end
end

if(DoCenter)
    OmegaNan = double(Omega); OmegaNan(~Omega) = nan;
    meanTCol = nanmean(T.*OmegaNan,2);
    clear OmegaNan
    T = bsxfun(@minus,T,meanTCol);
    
    normT2 = sum((T(:).*Omega(:)).^2);
else
    normT2 = sum(((T(:)-meanT).*Omega(:)).^2);
end

%numKnownCorrupt=sum(vec(Omega(:,CorruptCols,:)));
%maxVal2 = 255^2;

if(isscalar(epsilon))
    epsilonStop = epsilon;
    epsilonNoise = epsilon;
else
    epsilonStop = epsilon(1);
    epsilonNoise = epsilon(2);
end
if(epsilonNoise==0||isempty(CorruptCols)||~any(ismember(CorruptCols,1:size(T,2))))
    NoisyData=false;
    CorruptCols = zeros(0,class(CorruptCols));
    NoCorruptCols = cast(1:size(T,2),'like',CorruptCols);
else
    NoisyData=true;
    CorruptCols = intersect(CorruptCols,1:size(T,2));
    NoCorruptCols = setdiff(cast(1:size(T,2),'like',CorruptCols),CorruptCols);
    if(isempty(NoCorruptCols))
        AllCorrupt = true;
        normTCorrupt2 = normT2;
    else
        AllCorrupt = false;
        if(DoCenter)
            normTCorrupt2 = sum(vec((T(:,CorruptCols,:).*Omega(:,CorruptCols,:))).^2);
        else
            meanTCorrupt = sum(vec((T(:,CorruptCols,:).*Omega(:,CorruptCols,:)))) / sum(vec(Omega(:,CorruptCols,:)));
            normTCorrupt2 = sum(vec(((T(:,CorruptCols,:)-meanTCorrupt).*Omega(:,CorruptCols,:))).^2);
        end
    end
end

if (~exist('X','var')||isempty(X))
    X = T;
    if(DoCenter)
        X(~Omega) = 0;
    else
        X(~Omega) = meanT;%X(logical(1-Omega)) = 256;
    end
    %OmegaNan=double(Omega);OmegaNan(~Omega)=nan;
    %X = T.*Omega + (1-Omega).*repmat(nanmean(T.*OmegaNan),size(T,1),1);
    %clear OmegaNan
end
if (~exist('TraceNormToRankParam','var')|| isempty(TraceNormToRankParam))
    TraceNormToRankParam=0;
end
if(~exist('maxRank','var')||isempty(maxRank))
    maxRank=inf;
end

if(TraceNormToRankParam==0)
    legacyThreshold = true;
end


%Check the valid dimensions for tensor unfolding
TensorDim = ndims(T);
if(TensorDim==2)
    TensorDim=1;
    alpha=1;
    ValidDims=1;
    numValidDims=1;
else
    ValidDims = 1:TensorDim; ValidDims(alpha(ValidDims)==0)=[];
    numValidDims = length(ValidDims);
end

dim = size(T);
M = cell(TensorDim, 1);
Y = M;
rank = zeros(TensorDim,1);
beta = zeros(1,TensorDim, class(T));

%trMult=@(x)x'*x;
for i = ValidDims
    M{i} = X;
    Y{i} = zeros(dim, class(T));
    %beta(i) = 1.25 / lansvd(Unfold(X,dim,i), 1, 'L');
    %SV12 = sqrt(svds(trMult(Unfold(X,dim,i)), 2));
    SV12 = firstSVs(Unfold(double(X),dim,i),2);
    if(legacyThreshold)
        beta(i) = 2 / (SV12(1)+SV12(2));    %initial threshold = 1/beta = mean(SV1+SV2);
        beta(i) = beta(i)*(alpha(i).^2);
    else
        beta(i) = 8 / (SV12(1)+SV12(2)).^2; %initial threshold = sqrt(2/beta) = mean(SV1+SV2);
        beta(i) = beta(i)*alpha(i);
    end
end
%beta(ValidDims) = max(beta)*ones(1,numValidDims);


Ysum = zeros(dim, class(T));
Msum = zeros(dim, class(T));
errList=[];
%equalityErrorPrev=-1;
fprintf('Iteration:                  ');
for k = 1: maxIter
    
    % update M
    Ysum = 0*Ysum;
    Msum = 0*Msum;
    for i = ValidDims
        if(legacyThreshold)
            [M{i}, rank(i), ~] = Pro2TraceNorm(Unfold(X-Y{i}/beta(i), dim, i), alpha(i)/beta(i), TraceNormToRankParam);
        else
            [M{i}, rank(i), ~] = Pro2TraceNorm(Unfold(X-Y{i}/beta(i), dim, i), alpha(i)*sqrt(2/beta(i)), TraceNormToRankParam);
        end
        M{i} = Fold(M{i}, dim, i);
        Ysum = Ysum + Y{i};
        Msum = Msum + beta(i)*M{i};
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%5i (rank =%4i)',k,rank(i));
    end
    
    if(any(rank > maxRank))
        break;
    end
    
    % update X (outside of Omega)
    X = (Msum + Ysum) / sum(beta);
    
    
    % compute the error and check stop criterion
    if(isa(T,'single')),[errList(end+1),errCorrupt2] = SQErrorsLight(X,T,Omega,CorruptCols);
    else,[errList(end+1),errCorrupt2] = SQErrors(X,T,Omega,CorruptCols);
    end
    errList(end) = errList(end) / normT2;
    %errList(end+1) = norm((X-T).*Omega,'fro')^2 / normT2;
    %errList(end+1) = sum(((X(:)-T(:)).*Omega(:)).^2) / normT2;%/ numKnown / maxVal2;
    if (errList(end) < epsilonStop || k == maxIter )
    %if k == maxIter
        X_final = X;
        X_final(Omega) = T(Omega);
        if(DoCenter)
            X_final = bsxfun(@plus,X_final,meanTCol);
        end
        fprintf('\n');
        break;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update X (in Omega)
if(NoisyData) %l2 tolerance:
    % Allow a small error epsilon (in l2 norm) on the known elements of the matrix.
    if(AllCorrupt)
        lambda = max(sqrt( errList(end) / (epsilonNoise) ) - 1, 0);
        X = X.*(1-Omega) + ((X + lambda*T)/(1+lambda)).*Omega;
        %X(Omega) = ( X(Omega) + lambda * T(Omega) ) / (1 + lambda); %Much slower by indexing with Omega.
    else
        errCorrupt2 = errCorrupt2 / normTCorrupt2;
        %errCorrupt2 = norm((X(:,CorruptCols,:)-T(:,CorruptCols,:)).*Omega(:,CorruptCols,:),'fro')^2 / normTCorrupt2;
        %errCorrupt2 = sum(vec((X(:,CorruptCols,:)-T(:,CorruptCols,:)).*Omega(:,CorruptCols,:)).^2)/ normTCorrupt2;%numKnownCorrupt / maxVal2;
        lambda = max(sqrt( errCorrupt2 / (epsilonNoise) ) - 1, 0);
        
        if(isa(T,'single')), UpdateXOmegaLight(X,T,Omega,lambda,CorruptCols); %faster with mex.
        else, UpdateXOmega(X,T,Omega,lambda,CorruptCols);end
        %X(:,CorruptCols,:) = X(:,CorruptCols,:).*(1-Omega(:,CorruptCols,:)) + ((X(:,CorruptCols,:) + lambda*T(:,CorruptCols,:))/(1+lambda)).*Omega(:,CorruptCols,:);
        %X(:,NoCorruptCols,:) = X(:,NoCorruptCols,:).*(1-Omega(:,NoCorruptCols,:)) + T(:,NoCorruptCols,:).*Omega(:,NoCorruptCols,:);
    end
else
    X = X.*(1-Omega) + T.*Omega;
    %X(Omega) = T(Omega);%Much slower by indexing with Omega.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Y
    for i = ValidDims
        Y{i} = Y{i} + beta(i)*(M{i} - X);
    end
    
    % update beta
    beta  = beta  * betaMult;
    %{
    %Adaptive update strategy, tuned for the inpainting problem (faster convergence but hard to justify...)
    equalityError = sqrt(sum( (Ysum(:)/sum(beta)).^2 ));
	%mult = .5+min(max(equalityError/equalityErrorPrev,.5),10);
    th0=.5; th1=10; multmin=1; multmax=10;
    mult = (min(max(equalityError/equalityErrorPrev,th0),th1) - th0) * (multmax-multmin) / (th1-th0) + multmin;
    equalityErrorPrev = equalityError;
    %if(increasePenalty)
        beta  = beta  * mult;
    %end
    %}
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SVs=firstSVs(X,numSVs)
    [m, n] = size(X);
    if( m > 2*n )
        SVs = sqrt(svds(X'*X,numSVs));
    elseif( n > 2*m )
        SVs = sqrt(svds(X*X',numSVs));
    else
        SVs = svds(X,numSVs);
    end
end
  
