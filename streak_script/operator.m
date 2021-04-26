function [RHS,LHS] = operator(U, metric, D1, D2, beta, omega, Re, eulco, qm2, qm1)

% OPERATOR sets up the operator corrsponding to local and nonlocal stability
% equations.
%
% [RHS,LHS]=OPERATOR(U, D1, D2, beta, omega, Re, eulco, q1, q2)
%
% Input
%	U: structure containing meanflow field and its derivatives
%	alpha: streamwise wavenumber
%	D1 and D2: matrices corresponding to first- and second derivative operators
%	beta: spanwise wavenumber
%	omega: angular frequency
%	Re: Reynolds number
%	eulco: coefficients of Euler backward discretization scheme
%	qm1 and qm2: amplitudefunctions at x(i-1) and x(i-2)
%
% Output
%	LHS and RHS: matrices corresponding to PSE, LHS.q=RHS
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%
%------------------------------------------------------------------
N=length(U.u);
%------------------------------------------------------------------
% Set up operator matrices
% A.q + B.q_y + C.q_yy + D.q_x=0

i = sqrt(-1);
I = eye(size(D1));
O = zeros(size(D1));
xi=diag(-i*omega + (beta^2)/Re + i*(beta*U.w));

h1 =metric.h1;
m13=metric.m13;
c0 = 1; % set to 0 for restrict Görtler problem

A12 = diag(m13)*c0;
A13 = i*beta*I;
A21 = xi+diag(U.ux);
A22 = diag(U.uy)+diag(m13)*diag(U.u)*c0;
A31 = diag(U.vx)-2*diag(m13)*diag(U.u);
A32 = xi+diag(U.vy);
A41 = diag(U.wx);
A42 = diag(U.wy);
A43 = xi;
A44 = i*beta*I;

A = [O,   A12, A13, O  ; ...
     A21, A22, O  , O  ; ...
     A31, A32, O  , O  ;
     A41, A42, A43, A44]; 
    
%------------------------------------------------------------------

VD1 = diag(U.v)*D1;

B = [O,   D1,  O,   O;  ...
     VD1, O,   O,   O;  ...
     O,   VD1, O,   D1; ...
     O,   O,   VD1, O];

%------------------------------------------------------------------

C = [O,      O,      O,      O; ...
     -D2/Re, O,      O,      O; ...
     O,      -D2/Re, O,      O; ...
     O,      O,      -D2/Re, O];

%------------------------------------------------------------------    
    u = diag(U.u);
    
    D = [I, O, O, O; ...
         u, O, O, O; ...
         O, u, O, O; ...
         O, O, u, O];
%------------------------------------------------------------------

%Left-hand-side and Right-hand-side of equation system
invh1=blkdiag( diag(1./h1), diag(1./h1), diag(1./h1), diag(1./h1) );
D=invh1*D;

LHS = A + B + C + D*eulco(3);
RHS = -D*(eulco(1)*qm2 + eulco(2)*qm1);

%---------------------------------------------------------------------
%Boundary conditions
Ov = zeros(1,4*N);

% u(0)=0
LHS(2*N,:) = Ov;
LHS(2*N, N) = 1;
RHS(2*N) = 0;

% v(0)=0
LHS(3*N,:) = Ov;
LHS(3*N, 2*N) = 1;
RHS(3*N) = 0;

% w(0)=0
LHS(4*N,:) = Ov;
LHS(4*N, 3*N) = 1;
RHS(4*N) = 0;

bcase=0; % 1 for dq/dy=0, 0 for q=0
if bcase==0
    % u(ymax)=0
    LHS(N+1,:) = Ov;
    LHS(N+1,1) = 1;
    RHS(N+1) = 0;
    
    % v(ymax)=0
    LHS(2*N+1,:) = Ov;
    LHS(2*N+1, N+1) = 1;
    RHS(2*N+1) = 0;
    
    % w(ymax)=0
    LHS(3*N+1,:) = Ov;
    LHS(3*N+1, 2*N+1) = 1;
    RHS(3*N+1) = 0;
    
else
    % du/dy(ymax)=0
    LHS(N+1,:) = Ov;
    LHS(N+1,1:N) = D1(1,:);
    RHS(N+1) = 0;
    
    % dv/dy (ymax)=0
    LHS(2*N+1,:) = Ov;
    LHS(2*N+1,N+1:2*N) = D1(1,:);
    RHS(2*N+1) = 0;
    
    % dw/dy(ymax)=0
    LHS(3*N+1,:) = Ov;
    LHS(3*N+1,2*N+1:3*N) = D1(1,:);
    RHS(3*N+1) = 0;
    
end

end
