function [phi,dphidx,dphidxdx] = mlsQuadricBasis1D(pt,index,node,di,form)
% Compute the MLS shape function at point pt for all nodes within the
% support of this point pt.
% Basis used is linear basis pT = [1 x]
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

% --------------------------------------
%      compute the moment matrix A
% --------------------------------------
ni   = length(index);
A    = zeros(3,3);
dAdx    = zeros(3,3);
dAdxdx  = zeros(3,3) ;
phi  = zeros(ni,1);
dphidx  = zeros(ni,1);
dphidxdx = zeros(ni,1);
w    = zeros(ni,1);
dwdx    = zeros(ni,1);
dwdxdx =  zeros(ni,1);
for m = 1 : ni
  idx  = index(m);
  xi   = node(idx);
  [wi,dwidx,dwidxdx]  = computeCircleSpline(pt,xi,di(idx),form);
  pTp  = [1 xi (xi)^2]'*[1 xi (xi)^2] ;
  A    = A    + wi*pTp ;
  dAdx = dAdx + dwidx*pTp ;
  dAdxdx = dAdxdx + dwidxdx*pTp ;
  w(m) = wi ;
  dwdx(m) = dwidx ;
  dwdxdx(m) = dwidxdx ;
end

p  = [1; pt; (pt)^2];
% --------------------------------------
%         compute  matrix c(x)
% --------------------------------------
% A(x)c(x)   = p(x)
% Backward substitutions, two times for c(x)
c  = zeros(3,3);
% Using LU factorization for A
[L,U,PERM] = lu(A);

for i = 1 : 3
    if i == 1         % backward substitution for c(x)
        C = PERM*p;
    elseif i == 2     % backward substitution for c,x(x)
        C = PERM*([0 1 2*pt]' - dAdx*c(1:3,1)); 
    elseif i == 3     % backward substitution for c,x(x)
        C = PERM*([0 0 2]' - (2*dAdx*c(1:3,2) + dAdxdx*c(1:3,1) ) ); 
    end

    D1 = C(1);
    D2 = C(2) - L(2,1)*D1; 
    D3 = C(3) - L(3,1)*D1 - L(3,2)*D2 ;

    c(3,i) = D3/U(3,3);
    c(2,i) = (D2 - U(2,3)*c(3,i))/(U(2,2)); 
    c(1,i) = (D1 - U(1,2)*c(2,i) - U(1,3)*c(3,i))/(U(1,1));
end

for m = 1 : ni
  xi  = node(index(m));
  piT = [1;xi;xi^2];
  phi(m) = c(:,1)'* piT*w(m);
  dphidx(m) = c(:,2)'*piT*w(m) +  c(:,1)'*piT*dwdx(m);
  dphidxdx(m) = c(:,3)'* piT*w(m) + 2*c(:,2)'*piT*dwdx(m) + c(:,1)'*piT*dwdxdx(m) ;
end





