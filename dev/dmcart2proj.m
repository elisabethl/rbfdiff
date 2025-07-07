function [DPx,DPy,DPz] = dmcart2proj(x,y,z,Dx,Dy,Dz)
% Create projected differentiation matrices from Cartesian ones
n = length(x);
spd = @(v) spdiags(v,0,n,n);
DPx = spd(1-x.^2)*Dx+spd(-x.*y)*Dy+spd(-x.*z)*Dz;
DPy = spd(-x.*y)*Dx+spd(1-y.^2)*Dy+spd(-y.*z)*Dz;
DPz = spd(-x.*z)*Dx+spd(-y.*z)*Dy+spd(1-z.^2)*Dz;