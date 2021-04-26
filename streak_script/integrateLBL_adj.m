function [qadj] = integrateLBL_adj(qadjN, x, y, U, metric, D1, D2,...
    beta, omega, Re, Nstation, FDorder)
% [q,alphax,ncalc] = integratePSE(q0, alpha0, x, y, U, D1, D2, W,...
%                           beta, omega, Re, Nstation, stab, dpdx, FDorder, iaux)
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
qadj=zeros(N*4,Nstation);        
qadj(:,Nstation) = qadjN;

dx= [-diff(x) 0];

%integration loop
n=Nstation;
fprintf('\n\t Integrate adjoint equation. Station no.: ');
while (n >1)
    n=n-1;
    if (mod(n,floor(Nstation/10))==0); fprintf('%i ',n); end
    
    %decide which FDorder for backward scheme
    if (n == Nstation-1 || strcmp(FDorder, 'first'))
        scheme_order = 'first';
    elseif (n ~= Nstation-1 && strcmp(FDorder, 'second'))
        scheme_order = 'second';
    end
    
    %get coefficients for backward scheme
    [eulco] = compeulco(scheme_order, dx(n), dx(n+1));
    
    %set up operator matrices with boundary conditions
    [RHS, LHS] = operator_adj(U(n), metric(n), D1, D2,...
        beta,omega,Re,eulco,qadj(:,min(Nstation,n+2)),qadj(:,n+1));
    
    %solve equation system
    qadj(:,n) = LHS\RHS;
end

end

