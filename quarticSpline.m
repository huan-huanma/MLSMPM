function [w,dwdr,dwdrr,dwdrrr] = quarticSpline(r)
% Compute cubic spline function

if (r <= 1.0)
   w    = 1 - 6*r*r + 8*r*r*r - 3*r*r*r*r;
   dwdr = -12*r + 24*r*r - 12*r*r*r ;
   dwdrr = -12 + 48*r - 36*r*r;
   dwdrrr = 48 - 72*r;
else
   w    = 0.0 ;
   dwdr = 0.0 ;
   dwdrr = 0.0;
   dwdrrr = 0.0;
end
