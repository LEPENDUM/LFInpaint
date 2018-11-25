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
% link : http://www.academia.edu/9019268/A_multivariate_generalization_of_the_power_exponential_family_of_distributions

function X = MGGD_generation(N, p, Sigma, beta)
% GENERATE_MGGD
%
% Main function that generates MGGD distributed sources.
%
% INPUT: 
%   - N : Number of realizations
%   - p : Dimension 
%   - rho : Correlation parameter
%   - beta : Shape parameter 
% 
% OUTPUT:
%   - X : pxN matrix containing MGGD sources 

X = transpose(RandSphere(N,p));

X = (Sigma)^(0.5)*X;

tau = (gamrnd(p/(2*beta),2,1,N)).^(1/(2*beta));
tau = repmat(tau, p, 1);

X = tau.*X;


return


function X=RandSphere(N,dim)
% RANDSPHERE
%
% RandSphere generates uniform random points on the surface of a unit radius
% N-dim sphere centered in the origin . This script uses differents algorithms
% according to the dimensions of points:
%
%    -2D:  random generation of theta [0 2*pi]
%    -3D:  the "trig method".
%    -nD:  Gaussian distribution
%
%
% INPUT:
%
%    N: integer number representing the number of points to be generated
%    dim: dimension of points, if omitted 3D is assumed as default
%
% OUTPUT:
%
%   X: Nxdim double matrix representing the coordinates of random points
%   generated
%
% EXAMPLE:
% 
%   N=1000;
%   X=RandSphere(N);
%   hold on
%   title('RandSphere')
%   plot3(X(:,1),X(:,2),X(:,3),'.k');
%   axis equal
%
%Authors: Luigi Giaccari,Ed Hoyle

switch dim
    case 3 %3D
        
        %trig method
        X=zeros(N,dim);%preallocate
        X(:,3)=rand(N,1)*2-1;%z
        t=rand(N,1)*2*pi;
        r=sqrt(1-X(:,3).^2);
        X(:,1)=r.*cos(t);%x
        X(:,2)=r.*sin(t);%y
    case 2 %2D
        
        %just a random generation of theta
        X(:,2)=rand(N,1)*2*pi;%theta: use y as temp value
        X(:,1)=cos(X(:,2));%x
        X(:,2)=sin(X(:,2));%y
    
    otherwise %nD
        
        %use gaussian distribution
        X=randn(N,dim);
        X=bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
end 
return