function [ curv ] = curvature( x,y )
%[ curv ] = curvature( x,y ) 
%   computes the curvature of curve y(x)
x=x(:); y=y(:);
s=[0; cumsum(sqrt(diff(x).^2+diff(y).^2))];
[L1,L2]=diffop(s);
xt=L1*x;
yt=L1*y;
xtt=L2*x;
ytt=L2*y;
curv=abs(xt.*ytt-yt.*xtt)./(xt.*xt+yt.*yt).^1.5;
end

