function [index] = defineSupport(node,x,di)
% find nodes in neighbouring of point x
% Inpputs:
%  node : numnode x 1, nodal coordinates
%  x    : scalar, coordinate of point 
%  di   : 1 x numnode, size of support of nodes
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

numnode = size(node,1) ;
dif     = node - ones(numnode,1)*x;
for i = 1 : numnode
    r(i) = abs(dif(i,:));
end 
index   = find(r' - di <= 0.00001);
