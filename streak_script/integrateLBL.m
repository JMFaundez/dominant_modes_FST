function [q] = integrateLBL(q0, x, y, U, metric, D1, D2, beta, omega, Re, Nstation, FDorder)
% [q,alphax,ncalc] = integratePSE(q0, x, y, U, D1, D2, beta, omega, Re, Nstation, FDorder)
% Thuis function integrates the Parabolised Stability Equations (PSE).
%
% Input
%    q0 : complex eigenfunction 
%    x, y: arrays containing streamwise and normal coordinates of the grid
%    U: structure array containing the meanflow quantities
%    D1, D2: first- and second derivative operators
%    beta, omega, Re: spanwise wavenumber, frequency and Reynolds number
%    Nstation: number of stations in x-direction to compute
%    FDorder: order of backward Euler discretization scheme (can be 'first' or 'second')
%
% Output
%    q: complex array contatining eigen functions (u,v,w,p) at computed x-stations
%    ncalc: number of sucessfully computed stations
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

%initialise variables
N=length(y);
q=zeros(N*4,Nstation);        
q(:,1) = q0;

dx= [0 diff(x)];

% integration loop
n=1;
fprintf('\n\t Integrate state equation. Station no.: ');
while (n < Nstation)
    n=n+1;
    if (mod(n,floor(Nstation/10))==0); fprintf('%i ',n); end
    
    %decide which FDorder for backward scheme
    if (n == 2 || strcmp(FDorder, 'first'))
        scheme_order = 'first';
    elseif (n ~= 2 && strcmp(FDorder, 'second'))
        scheme_order = 'second';
    end
    
    %get coefficients for backward scheme
    [eulco] = compeulco(scheme_order, dx(n), dx(n-1));
    
    %set up operator matrices with boundary conditions
    [RHS, LHS] = operator(U(n), metric(n), D1, D2,...
        beta,omega,Re,eulco,q(:,max(1,n-2)),q(:,n-1));
    
    %solve equation system
    q(:,n) = LHS\RHS;
end
end

