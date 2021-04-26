function [ L1, L2 ] = diffop(y)
% [ L1, L2 ] = diffop( y )
% This routine sets up arrays L1 and L2 which contain coefficients
% of a standard second order FD scheme for discretization of first
% and second derivatives, respectively.
%
% (c) Ardeshir Hanifi, 2011.
%
n=length(y);
L1=spalloc(n,n,3*n);
L2=spalloc(n,n,3*n);
%
% set up 2nd order FD coeff. at lower boundary
%
d1 = y(2)-y(1); d2 = y(3)-y(1);
L1(1,1:3) = ([ -(d1+d2)/d1/d2, -d2/d1/(d1-d2),  d1/d2/(d1-d2)  ]);
L2(1,1:3) = ([ 2.0/d1/d2,      2.0/d1/(d1-d2), -2.0/d2/(d1-d2) ]);

%
%     set up 2nd order FD coeff. at iner points
%
for j=2:n-1
    
    d1 = y(j-1)-y(j);
    d2 = y(j+1)-y(j);
    
    L1(j,j-1) = -d2/d1/(d1-d2);
    L1(j,j) = -(d1+d2)/d1/d2;
    L1(j,j+1) =  d1/d2/(d1-d2);
    
    L2(j,j-1) =  2.0/d1/(d1-d2);
    L2(j,j) =  2.0/d1/d2;
    L2(j,j+1) = -2.0/d2/(d1-d2);
    
end
%
%     set up 2nd order FD coeff. at upper boundary
%
d1 = y(n-1)-y(n); d2 = y(n-2)-y(n);
L1(n,n-2:n) = ([ d1/d2/(d1-d2),   -d2/d1/(d1-d2), -(d1+d2)/d1/d2 ]);
L2(n,n-2:n) = ([ -2.0/d2/(d1-d2), 2.0/d1/(d1-d2), 2.0/d1/d2      ]);
end


