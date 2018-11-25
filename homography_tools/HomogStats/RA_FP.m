function [V, b, m, alpha] = RA_FP(X,b_fix,varargin)
% Function that estimates the parameters of a multivariate generalized
% Gaussian distribution (MGGD). 
% INPUT
%   - X : pxN observations drawn from a zero mean MGGD
%   - b_fix : Fixed value of beta. If b_fix == 0 then estimate beta, else
%   not.
% OUTPUT
%   - V : Normalized covariance matrix
%   - b : Shape parameter. If b_fix==0 then beta is the estimated value of
%   the shape parameter.
%   - m : Scale parameter.
%   - alpha : adaptive choice of step
%
% Coded by Lionel Bombrun (lionel.bombrun at u-bordeaux.fr);
%          & Zois Boukouvalas (zb1 at umbc.edu); 
% 
% Machine Learning for Signal Processing Lab (MLSP-Lab)
% http://mlsp.umbc.edu
%
% Reference:
%
% [1] Boukouvalas, Z. and Said, S. and Bombrun, L. and Berthoumieu, Y. and
% Adali, T., "A New Riemannian Averaged Fixed-Point Algorithm for MGGD
% Parameter Estimation," Signal Processing Letters, IEEE, vol.22, no.12, pp.2314-2318, December 2015.


%% Call internal test function when no arguments supplied
if nargin==0
   help RA_FP
   test_RA_FP
   return
end
 
%% Gather Options
Params=struct('NMaxIter',100,'ThresHold',0.01);

% load in user supplied options
Params=getopt(Params,varargin{:});

% Get the size of the data
[p, n] = size(X);

% Use MoM to get an initial value.
[V, b, Fail] = Moments_MGGD_MEP_X(X, b_fix);

if p == 1
    S = X.^2 * inv(V);
else
    S = dot(X,V\X);
end

tol = 100 * ones(p^2 + 1, 1);
kk = 0;
NIter = 0;

while (any(tol > Params.ThresHold) && (NIter < Params.NMaxIter))  % Repeat until tolerance is under 1 percent for every parameter (or NIter >= 100).
    
    NIter=NIter+1;
    kk = kk + 1;
    alpha = 1/(kk);
    
    if p == 1
        V_new = (b/n) * sum(S.^(b - 1) .* X.^2);
    else
        V_new = zeros(p);
        for i = 1:n
            V_new = V_new + ((S(i))^(b - 1) * X(:, i) * X(:, i)');
        end
        V_new = (V^((1/2))) * (( (V^(-0.5)) * (V_new) * (V^(-0.5)) )^(alpha)) * (V^((1/2))); %Riemannian Averages
    end
    
    V_new = p*V_new/trace(V_new);
    
    if p == 1
        S = X.^2 * inv(V);
    else
        S = dot(X,V_new\X);
    end
    
    if b_fix ~= 0
        b_new = b_fix;
    else
        
        b_new = Newton_Raphson(b, S, p, n);
    end
    
    tol = 100 * [reshape(abs((V_new - V)./V), p^2, 1); abs((b_new - b)/b)];  % Relative change in percent
    V = V_new;
    b = b_new;
    
end

m = calc_m(b, S, p, n);

return;

function test_RA_FP
disp('Executing internal test function for RA-FP.')

%% Generate samples from an MGGD
N = 1000; % number of samples
p = 3; % dimension of the MGGD (scatter matrix)
rho = 0.5; %Correlation
beta =.5; %Shape Parameter

%Scatter matrix
Sigma = eye(p);
for i=1:p
    for j=1:p
        Sigma(i,j) = rho^(abs(i-j));
    end
end

X = MGGD_generation(N, p, Sigma, beta);

b_fix =0; % Put 0 if we want to estimate beta

[V_Mom, beta_hat, Fail] = Moments_MGGD_MEP_X(X, b_fix); %This is the MoM
errorMoM = norm(V_Mom - Sigma,'fro');

[V_RA_FP, beta_hat, m, alpha] = RA_FP(X, b_fix); %This is RA-FP
V_RA_FP = m*V_RA_FP;
errorRA_FP = norm(V_RA_FP - Sigma,'fro');

disp(['The error of estimating the scatter matrix using MoM is ' num2str(errorMoM)])
disp(['The error of estimating the scatter matrix using RA-FP is ' num2str(errorRA_FP)])

return

%%%%%%% More functions %%%%%%%

function m = calc_m(b, S, p, n)

F1 = sum(S.^b);
m = ((b*F1)/(p*n))^(1/b);

return

function [V, b, Fail] = Moments_MGGD_MEP_X(X, b_fix)
% Same as Moments_MMD, but using an input matrix X (p x n) for the data.
% 
% Calculates estimates of the parameters beta and V of a
% multivariate generalized Gaussian pdf fit to the data in WIm, which is
% 3-band colour data, using the method of moments. The image with index
% NIm is chosen for a single wavelet subband (NS). If b_fix ~= 0 then a
% fixed shape parameter beta is used (1 for Gaussian, 0.5 for Laplace).
% Otherwise, to solve the equation for beta, the half-interval method is
% used. Col is a vector of colour indices (1 = r, 2 = g, 3 = b) indicating
% the colour bands that one actually wishes to model. This also determines
% the dimensionality of the probability space p.
%
% ATTENTION: MEP version!!

    Fail = 0;
    [p, n] = size(X);    
    Var = zeros(p);  % Sample variance
    for i = 1:n
        Var = Var + X(:, i) * X(:, i)';
    end 
    Var = Var/n;
    
    
    if b_fix == 0
        IV = inv(Var);
        gamma2 = 0;  % Sample kurtosis
        for i = 1:n
            gamma2 = gamma2 + (X(:, i)' * IV * X(:, i))^2;
        end
        gamma2 = (1/n) * gamma2 - p * (p + 2);
        
        [b, IsMaxIter] = Zero_b(gamma2, p, n);
    else
        b = b_fix;
        IsMaxIter = 0;  % No Maximum Likelihood necessary.
    end
    
    if IsMaxIter == 1
        Fail = 1;
    end     
    
    V = (p * gamma(p/(2*b)) / (2^(1/b) * gamma((p + 2)/(2*b)))) * Var;                           
return


function [b, IsMaxIter] = Zero_b(gamma2, p, n)

    MaxIter = 200;
    IsMaxIter = 0;    
    R = [0.1, 5];  % Initial interval for b   
    Err = 0;
    sa = sign(g(R(1), gamma2, p, n));
    sb = sign(g(R(2), gamma2, p, n));  
    if sa == sb
        Err = 1;
        b = 1;
    else
        gb = 1;
        NIter = 0;        
        while (abs(gb) > 1e-3) && (NIter < MaxIter)
            NIter = NIter + 1;
            b = mean(R);  % Half-point
            gb = g(b, gamma2, p, n);
            if sb * gb > 0
                R(2) = b;
            elseif sb * gb < 0
                R(1) = b;
            elseif gb == 0
                break;  % Found an exact root.
            end
        end
        if NIter == MaxIter
%             disp('Moments_MGGD: maximum number of iterations reached. Results may be unreliable.');
            IsMaxIter = 1;
%            b = 0.5;
        end        
    end
return


function Res = g(b, gamma2, p, n)
% b = beta
   
    Res = p^2 * gamma(p./(2*b)) .* gamma((p + 4)./(2*b)) ...
        - (gamma((p + 2)./(2*b))).^2 * (p * (p + 2) + gamma2);
return


function b_new = Newton_Raphson(b, S, p, n)
% NEWTON_RAPHSON for estimating the shape parameter
% Input :
%   - b : Initial value for the shape parameter.
%   - S : Quadratic term.
%   - p : Dimension of the space.
%   - n : Number of obsevations.
% Output
%   - b_new : Estimated value of the shape parameter.

    b_new = b - f(b, S, p, n)/gi(b, S, p, n);
    if(b_new<=0)
        b_new = b;
    end
    if(imag(b_new)~=0)
        b_new = b;
    end
return

function Res = f(b, S, p, n)   
    F1 = sum(log(S) .* S.^b);
    F2 = sum(S.^b);
    Res = p*n*F1/(2*F2) - p*n/(2*b)*(psi(p/(2*b))+log(2)) - n - p*n/(2*b)*log(b*F2/(p*n));
    Res = Res/n;
return

function Res = gi(b, S, p, n)
    F1 = sum(S.^b);
    F2 = sum((log(S)) .* S.^b);
    F3 = sum((log(S)).^2 .* S.^b);
    
    Res = p*n/2 * ((F1*F3 - F2^2)/(F1^2)) + p*n/(2*b^2)*(psi(p/(2*b))+log(2)) - p*n/(2*b)*(-p*psi(1,p/(2*b))/(2*b^2));
    Res = Res + p*n/(2*b^2)*log(b/(p*n)) - p*n/(2*b^2) + p*n/(2*b^2)*log(F1) - p*n/(2*b)*F2/F1;
    Res = Res/n;
return


% MGGD_GENERATION
%
% This matlab script contains three functions that can be used to generate 
% multivariate generalized Gaussian (MGGD) distributed sources. MGGD
% depends on two parameters, the scatter matrix (which is symmetric
% positive definite) and the shape parameter. If shape parameter is less
% than 1 the distribution of the marginals is more peaky, with heavier
% tails, and if it is greater than 1, it is less peaky with lighter tails.
% If shape parameter is equal to one then we generate Gaussian sources.
%
% EXAMPLE:
% 
% X = MGGD_generation(1000, 3, 0.5, 1);
%
% Coded by Lionel Bombrun (lionel.bombrun at u-bordeaux.fr);
%          Zois Boukouvalas (zb1 at umbc.edu);
%
% References: E. Gomez, M. Gomez-Viilegas, and J. Marin, "A multivariate generalization of the power exponential family of distributions," Communications in Statistics-Theory and Methods, vol. 27, no. 3, pp. 589{600}, 1998.


function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.
%%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties =
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})

% Function from
% http://mathforum.org/epigone/comp.soft-sys.matlab/sloasmirsmon/bp0ndp$crq5@cui1.lmms.lmco.com

% dgleich
% 2003-11-19
% Added ability to pass a cell array of properties

if ~isempty(varargin) && (iscell(varargin{1}))
   varargin = varargin{1};
end;

% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
   arg = varargin{ii};
   if isempty(TargetField)
      if ~ischar(arg)
         error('Property names must be character strings');
      end
      %f = find(strcmp(prop_names, arg));
      if isempty(find(strcmp(prop_names, arg),1)) %length(f) == 0
         error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
      end
      TargetField = arg;
   else
      properties.(TargetField) = arg;
      TargetField = '';
   end
end
if ~isempty(TargetField)
   error('Property names and values must be specified in pairs.');
end
return


