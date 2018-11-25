%% Apply homography transform h to the input image I0.
% Inputs:
% - H : 3x3 homography matrix such that cst*[x_target, y_target, 1] = [x_o, y_o, 1]*H
% - I0 : input image to warp.
% - type : (optional) interpolation method (default = 'cubic').
% - defaultVal : default value for extrapolation on image borders. Can be
% set to 'auto' instead of a scalar value to extrapolate the image border
% automatically.
%
% Outputs:
% - I_target : transformed image.
% - mask : binary image indicating unknown borders (as false) of the transformed image I_target.

function [I_target, mask] = homographyWarp(H, I0, type, defaultVal)

if(nargin<3)
    type='cubic';
end
if(nargin<4)
    defaultVal=0;
end

doFill = strcmp(defaultVal,'auto');

resY = size(I0,1);
resX = size(I0,2);
nChan = size(I0,3);

%% original coordinates. Attention: x and y are inversed 
coord_o = ones(resY*resX,3);
[X, Y] = meshgrid(1:resX, 1:resY);
coord_o(:,1) = vec(X);
coord_o(:,2) = vec(Y);

%% derive the target coordinates.
coord_t = coord_o * H;

%% remove the third dimension
coord_t = [coord_t(:,1)./coord_t(:,3), coord_t(:,2)./coord_t(:,3)];

%% grid interpolation
X_Target = reshape(coord_t(:,1),[resY,resX]);
Y_Target = reshape(coord_t(:,2),[resY,resX]);

if(nargout>1)
    mask = X_Target<=resX & X_Target>=1 & Y_Target<=resY & Y_Target>=1;
end

if(doFill)
    X_Target = min(max(X_Target,1),resX);
    Y_Target = min(max(Y_Target,1),resY);
end

I_target = zeros(resY,resX,nChan);
for ch=1:nChan
    I_target(:,:,ch) = interp2(X, Y, I0(:,:,ch), X_Target, Y_Target, type);
end

if(~doFill)
    I_target(isnan(I_target)) = defaultVal;
    %mask = ~isnan(I_target(:,:,1));
end

end