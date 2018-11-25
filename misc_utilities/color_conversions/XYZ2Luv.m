%  Conversion function : XYZ to Lu'v' (L=Y / CIE (u' v') system)
%  Lu'v' is the same space as in LogLuv without the log encoding of luminance L.
function [Luv] = XYZ2Luv(XYZ)

Luv(:,:,1) =  XYZ(:,:,2);
u = 4 * XYZ(:,:,1) ./ ( XYZ(:,:,1) + 15 * XYZ(:,:,2) + 3 * XYZ(:,:,3) );
v = 9 * XYZ(:,:,2) ./ ( XYZ(:,:,1) + 15 * XYZ(:,:,2) + 3 * XYZ(:,:,3) );

u(isnan(u)) = 0.1978;
v(isnan(v)) = 0.4683;

Luv(:,:,2) = u;
Luv(:,:,3) = v;