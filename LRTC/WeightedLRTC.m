%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Low Rank Tensor/Matrix Completion
%
% min_M: rank(M)
% s.t.:  ||W o (M - T)||_F^2 <= epsilon
%
% Resolution is performed using ADMM by introducing an intermediate matrix X and by solving:
% min_M: rank(M)
% s.t.:  M = X
%        ||W o (X - T)||_F^2 <= epsilon
%
% May still be compatible for tensors but not tested. (minimizing the averaged ranks of
% matrices unfolded along the different tensors dimensions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inputs:
% - T : input matrix to complete/approximate.
% - Weights : element-wise weights. Input weights should be in the range [0 1]. Note: the
% actual weights values used in the minimization are quantized and remapped.
% - numWeightVals : number of quantization steps for the Weights (integer between 2 and 255).
% - alpha : (only for tensors) weights on the different tensor dimensions.
% - betaMult : parameter controling convergence speed vs accuraty tradeoff.
% It must be higher than 1. (using betaMult=1 may never converge).
% - maxIter : maximum number of iterations.
% - epsilon : global noise tolerance parameter. The parameter may be given as a vector of 2-elements : 
%            - 1st element = epsilon value used for stopping criterion /
%            - 2nd element = epsilon value used in the error constraint.
%            Note that the input epsilon values are normalized (divided by the input matrix weighted norm).
% - TraceNormToRankParam : 0 (default) => trace norm minimization (soft thresholding of singular values) /
%                          inf => rank minimization /
%                          Other positive values should give intermediate results.
% - maxRank : maximum rank (default=inf). The algorithm stops when maxRank is reached.
% - X : User defined initialization of the approximation matrix. (default = T_ij if W_ij>0 / mean(T) otherwise).

% Outputs:
% - X : Weighted Low Rank Approximation result.
% - M : Cell containing matrix M (In case of tensors, the cell contains all the matrices unfolded along each tensor dimension).
%       A way to verify the convergence is to check that X is sufficiently close to M.
%       M has exactly the rank obtained at the last iteration, but may not satisfty the constraint on weighted norm if convergence is not reached.
% - rank : values of the matrix rank (for tensors: list of ranks of the matrices unfolded along each tensor dimension).
% - errList : list of error at each iteration (if the algrithm converged errList(end) should be less than epsilon).

function [X, M, rank, errList] = WeightedLRTC(T, Weights, numWeightVals, alpha, betaMult, maxIter, epsilon, TraceNormToRankParam, maxRank, X)

DoCenter = false; %center the data of each column before completion.
legacyThreshold = false; %true  : threshold value determined for nuclear norm minimization.
                         %false : threshold value determined for rank minimization.
                         %This parameter is automatically set to true when TraceNormToRankParam=0 (nuclear norm minimization).
betaMult0 = betaMult;

%%%%%%%%%%%%%%%%%%%%%%%%  Check Input Parameters  %%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('TraceNormToRankParam','var')|| isempty(TraceNormToRankParam)), TraceNormToRankParam=0; end
if(~exist('maxRank','var')||isempty(maxRank)), maxRank=inf; end
if(isscalar(epsilon)), epsilon=[epsilon,epsilon];end
if(epsilon(1)<0||epsilon(2)<0),error('epsilon parameter should not be negative!');end


%%%%%%%%%%%%%%%%%%%%%%%  Prepare other parameters  %%%%%%%%%%%%%%%%%%%%%%%%
if(TraceNormToRankParam==0), legacyThreshold = true; end
HasNoiseTolerance = epsilon(2)>0;

%Find Valid dimensions for tensor unfolding
TensorDim = ndims(T);
if(TensorDim==2)%2-way tensor = matrix => no unfolding necessary.
    TensorDim=1;
    alpha=1;
    ValidDims=1;
else
    ValidDims = 1:TensorDim; ValidDims(alpha(ValidDims)==0)=[];
end

%Quantize the Weights so that the number of different weight values is at most numWeightVals (for efficient computation of lagrangian parameter lambda).
if(HasNoiseTolerance)
    numWeightVals = max(min(numWeightVals,255),2);
    WLabels = uint8(min(Weights,1)*(numWeightVals-1));
    Weights = double(WLabels)/(numWeightVals-1);
    AllWeights2=([0:1/(numWeightVals-1):1]').^2;%list of squared weights (only the squared weight values are needed in the iterations).
    maxLabel = uint8(numWeightVals-1);
    
%    tol = .01;
%    Weights = (1+tol)./(1+tol-Weights) - 1;
%    AllWeights2 = ((1+tol)./(1+tol-([0:1/(numWeightVals-1):1]')) - 1).^2;
    
    %epsilon = epsilon*tol.^2;
%    epsilon(2) = epsilon(2)*tol.^2;
    %epsilon(1) = epsilon(1)/tol.^2;
    
    %%%%%%%%%%%%%%%%%%%%%%
    %{
    Weights = WLabels==maxLabel;
    WLabels = uint8(Weights);
    AllWeights2 = [0;1];
    maxLabel=uint8(1);
    numWeightVals=2;
    %}
    %%%%%%%%%%%%%%%%%%%%%%
    
else %if there is no tolerance to noise, we are back to the non-weigthed matrix completion case (binary weights).
    Weights = Weights>0;
end

%Compute mean, centering and norm of input data.
meanT = sum(T(:).*Weights(:)) / sum(Weights(:));
if(DoCenter)
    meanTCol = sum(T.*Weights,1) ./ sum(Weights,1);
    T = bsxfun(@minus,T,meanTCol);
    
    normT2 = sum((T(:).*Weights(:)).^2);
else
    normT2 = sum(((T(:)-meanT).*Weights(:)).^2);
end

errorStop = epsilon(1)*normT2;  %Error used as a stopping criterion (should not be lower than the error tolerance errorTol, otherwise the criterion is never met).
errorTol = epsilon(2)*normT2;   %Error tolerance on the weighted norm.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        ADMM Initialization        %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('X','var')||isempty(X))
    X = T;
    if(DoCenter)
        X(Weights==0) = 0;
    else
        X(Weights==0) = meanT;
    end
end

dim = size(T);
M = cell(TensorDim, 1);
Y = M;
rank = zeros(TensorDim,1);
beta = zeros(1,TensorDim);

for i = ValidDims
    M{i} = X;
    Y{i} = zeros(dim, class(T));
    SV12 = firstSVs(Unfold(X,dim,i),2);
    if(legacyThreshold)
        beta(i) = 2 / (SV12(1)+SV12(2));    %initial threshold = 1/beta = mean(SV1+SV2);
        beta(i) = beta(i)*(alpha(i).^2);
    else
        beta(i) = 8 / (SV12(1)+SV12(2)).^2; %initial threshold = sqrt(2/beta) = mean(SV1+SV2);
        beta(i) = beta(i)*alpha(i);
    end
end
%beta(ValidDims) = max(beta)*ones(1,length(ValidDims)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%            ADMM solving           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(HasNoiseTolerance)
    Weights = Weights.^2; %Only the squared weights values are needed within the iterations.
    %clear Weights
end
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
    
    % update X (1st step => without taking error constraint into account).
    X = (Msum + Ysum) / sum(beta);
    
    % compute the error and check stop criterion
    if(HasNoiseTolerance)
        SSDPerLabel = SQErrorsPerLabel(X,T,WLabels,maxLabel);
        errList(end+1) = SSDPerLabel*AllWeights2;
    else
        errList(end+1) = SQErrors(X,T,Weights);
    end
    
    %Stop criterion : Here the error used for the stop criterion is computed on the partially updated matrix X because next update step
    %requires computation of sum of squared errors, which is also useful for stop criterion evaluation. But in theory the matrix M (exactly
    %with rank displayed) should be used for the stop criterion (we simplified to avoid computing SSD for both M and X).
    if (errList(end) < errorStop || k == maxIter )
        errList = errList/normT2;%return normalized error so that it can be compared to the input epsilon parameter.
        if(DoCenter)
            X = bsxfun(@plus,X,meanTCol);
        end
        fprintf('\n');
        break;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update X (2nd step => satisfy constraint on the error).
if(HasNoiseTolerance) %l2 tolerance:
    % Allow an error epsilon (in l2 norm) on the known elements of the matrix.
    if(numWeightVals==2) %Particular case with closed form solution. Assumes that AllWeights2(1)==0 and AllWeights2(2)>0 (enforced in the definition of AllWeights2).
        lambda = max( (sqrt( errList(end) / (errorTol) ) - 1)/AllWeights2(2) , 0);
    else %find lambda with Newton's method.
        lambda = ComputeLambdaWeightedLRC(AllWeights2, SSDPerLabel, errorTol, 1e-4, 10000);
    end
    %X = (X + lambda*T.*Weights.^2) ./ (1+lambda*Weights.^2);
    %X = (X + lambda*T.*AllWeights2(1+WLabels)) ./ (1+lambda*AllWeights2(1+WLabels));
%    X = (X + lambda*T.*(double(WLabels)/(numWeightVals-1)).^2) ./ (1+lambda*(double(WLabels)/(numWeightVals-1)).^2);%Assuming Weights=double(WLabels)/(numWeightVals-1);
    X = (X + lambda*T.*Weights) ./ (1+lambda*Weights); %Assuming 'Weights' stores the squared weights;
else
    X = X.*(1-Weights) + T.*Weights;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Y
    for i = ValidDims
        Y{i} = Y{i} + beta(i)*(M{i} - X);
    end
    
    % update beta
%    betaMult = betaMult0 * (1+rank(i)/size(T,2)).^2;
    %betaMult = betaMult0 + 4*(errorTol/errList(end)).^2;
    beta  = beta  * betaMult;
    %{
    %Adaptive update strategy, tentative of tuning for the inpainting problem (faster convergence but hard to justify=> not used for the article...)
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
  
