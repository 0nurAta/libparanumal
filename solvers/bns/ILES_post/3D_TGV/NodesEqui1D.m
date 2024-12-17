function [r] =NodesEqui1D(N)

% total number of nodes
Np = (N+1);

zer=zeros(Np,1);
dummy=zer;
% Create equispaced nodes on a straight line.
for n=1:Np
    dummy(n) = -1 + (2*(n-1)) / N;
 end

r = dummy;
end
