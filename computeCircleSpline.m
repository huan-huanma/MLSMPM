
function [w,dwdx,dwdxx,dwdxxx] = computeCircleSpline(x,xI,d,form)
% Compute cubic and quartic spline function
% Inputs:
% x (1x2)  : coordinate of point at which w is to be evaluated
% xI (1x2) : coord of node I
% d        : size of the support

sigma = 7;  % 7

r = sqrt( (x  - xI ) .* (x - xI ) )/d;

switch form
  case 'cubic_spline' 
     [w,dwdr,dwdrr,dwdrrr] = cubicSpline(r);
  case 'quartic_spline'
     [w,dwdr,dwdrr,dwdrrr] = quarticSpline(r);
  case 'exp_spline'
     [w,dwdr,dwdrr,dwdrrr] = expSpline(r,sigma);
  otherwise 
     error('Grr. Unknown functional form');
end

if (r ~= 0)
    drdx = (x - xI)/(r*d*d) ; 
    drdxx =  1/(r*d*d) - (x - xI)^2/(r*r*r*d*d*d*d) ;
    drdxxx =  - 3*(x - xI )/(r*r*r*d*d*d*d) + 3*(x  - xI )^3/(r*r*r*r*r*d*d*d*d*d*d) ;  
else
    drdx = 0 ; 
    drdxx = 0;
    drdxxx = 0;
end

dwdx = dwdr * drdx ; 
dwdxx = dwdrr * drdx*drdx + dwdr * drdxx;
dwdxxx = dwdrrr * drdx*drdx*drdx + 3*dwdrr *drdx * drdxx + dwdr*drdxxx;