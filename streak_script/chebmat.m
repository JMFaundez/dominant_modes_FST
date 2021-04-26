function [y,D1,D2,W] = chebmat(N,ymax)
% [z,D1,D2,W] = chebmat(N,zmax)
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
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

[yc, D] = chebdif(N, 2);
D1 = D(:,:,1);
D2 = D(:,:,2);

% Define the mapping to the physical domain
y=(yc+1)*ymax/2;

% Map Derivative operators onto interval [0 ymax]
D1 = D1*(2/ymax);
D2 = D2*(2/ymax)^2;

% Get integration weights
W = iwt(N)*ymax/2;

end
