function [eulco] = compeulco(order, dx1, dx2)
%  [eulco] = compeulco(order, dx1, dx2) computes the coefficients for backward Euler FD scheme
%
% Input
%   order: 'first' or 'second' for first- or second-order backward Euler scheme
%   dx1: x(i)-x(i-1)
%   dx2: x(i-1)-x(i-2)
%
% Output:
%   eulco: an array containing the discretization coefficients such that 
%          df/dx=eulco(1)*f(i-2)+eulco(2)*f(i-1)+eulco(3)*f(i)
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

if strcmp(order, 'first')
    eulco=([0, -1, 1])/dx1;
elseif strcmp(order, 'second')
    d=dx2+dx1;
    eulco = ([dx1/(dx2*d), -d/(dx1*dx2), (2*dx1+dx2)/(dx1*d)]);
end
