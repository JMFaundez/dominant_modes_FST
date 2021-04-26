function [y,D1,D2,W] = chebmat_trans(N,ymax,yi)
% [y,D1,D2,W] = chebmat(N,ymax)
% This function sets up the integration weights and differential operatores
% corresponding to the spectral discretization using Chenychev polynomials.
% The differential operators are also transormed using a mapping from
% the spectral domain [-1,1] to [0,ymax].
%
% N: order of Chebychem polynomial is N-1
% ymax: Height of last node in wall-normal direction
%
% The output variable are
% y: coordinates of Gauss-Lobbato nodes in physical space [0,ymax]
% D1, D2: first and second derivative operatires
% W: spectral integraation weights
%
%
% (c) Ardeshir Hanifi, 2020
%

[x, D] = chebdif(N, 2);
D1 = D(:,:,1);
D2 = D(:,:,2);

% Transformation of grid points to physical space
a=ymax*yi/(ymax-2*yi);
b=1+2*a/ymax;
y=a*(1+x)./(b-x);
tranfac= (b-x).^2/(b+1)/a; % dx/dy

% Map onto interval [0 ymax]
D1 = diag(tranfac)*D1;
D2 = D1*D1;

% Get integration weights
W = iwt(N);

N1= N-1;
n = 0:1:N1;
j = 0:1:N1;
b = ones(1, N);
b([1 N]) = 0.5;
c = 2*b;
b = b/N1;
S = cos(n'*j*pi/N1);
dW= diag(W);
F = S*(dW./tranfac);
W = diag(b.*((c.*F')*S));

end
