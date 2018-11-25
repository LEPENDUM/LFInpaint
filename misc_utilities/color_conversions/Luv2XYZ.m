%  Conversion function : Lu'v' to XYZ (L=Y / CIE (u' v') system)
%  Lu'v' is the same space as in LogLuv without the log encoding of luminance L.
function [XYZ] = Luv2XYZ(Luv)


XYZ(:,:,2) =  Luv(:,:,1);   % Y = L

x = 9 * Luv(:,:,2) ./ ( 6 * Luv(:,:,2) - 16 * Luv(:,:,3) + 12 );
y = 4 * Luv(:,:,3) ./ ( 6 * Luv(:,:,2) - 16 * Luv(:,:,3) + 12 );

XYZ(:,:,1) = x .* XYZ(:,:,2) ./ y;       % X = x*(Y/y) 
XYZ(:,:,3) = (1-x-y) .* XYZ(:,:,2) ./ y; % Z = (1-x-y)*(Y/y)