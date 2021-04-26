function [RHS,LHS] = operator_adj(U, metric, D1, D2, beta, omega, Re, eulco, qm2, qm1)

% OPERATOR sets up the operator corrsponding to local and nonlocal stability
% equations.
%
% [RHS,LHS]=OPERATOR(U, D1, D2, beta, omega, Re, eulco, q1, q2)
%
% Input
%	U: structure containing meanflow field and its derivatives
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
xi=diag(i*omega + (beta^2)/Re - i*(beta*U.w));

h1 =metric.h1;
m13=metric.m13;
c0 = 1; % set to 0 for restrict Görtler problem

A12 = xi+diag(U.ux);
A13 = diag(U.vx)-2*diag(m13)*diag(U.u);
A14 = diag(U.wx);
A21 = diag(m13)*c0;
A22 = diag(U.uy)+diag(m13)*diag(U.u)*c0;
A23 = xi+diag(U.vy);
A24 = diag(U.wy);
A31 = -i*beta*I;
A34 = xi;
A44 = -i*beta*I;

A = [O,   A12,  A13, A14; ...
     A21, A22,  A23, A24; ...
     A31, O,    O,   A34;
     O,   O,    O,   A44];
%------------------------------------------------------------------
V   = diag(U.v);
VD1 = diag(U.v)*D1;
Vy  = diag(U.vy);

B = [O,   VD1, O,    O;  ...
     D1,  O,   VD1,  O;  ...
     O,   O,   O,    VD1; ...
     O,   O,   D1,   O ];
 
B0 = [O,  V,   O,  O;  ...
      I,  O,   V,  O;  ...
      O,  O,   O,  V; ...
      O,  O,   I,  O];

By = [O,  Vy, O,  O;  ...
      O,  O,  Vy, O;  ...
      O,  O,  O,  Vy; ...
      O,  O,  O,  O ];

%------------------------------------------------------------------

C = [O, -D2/Re,   O,    O; ...
     O,   O,   -D2/Re,  O; ...
     O,   O,      O,  -D2/Re; ...
     O,   O,      O,    O];
 
C1 = [O, -D1/Re,   O,    O; ...
      O,   O,   -D1/Re,  O; ...
      O,   O,      O,  -D1/Re; ...
      O,   O,      O,    O];

%------------------------------------------------------------------
u = diag(U.u);
ux = diag(U.ux);

D = [I, u, O, O ; ...
     O, O, u, O ; ...
     O, O, O, u ; ...
     O, O, O, O ];
Dx= [O, ux, O, O ; ...
     O, O, ux, O ; ...
     O, O,  O, ux ; ...
     O, O,  O, O ];
%---------------------------------------------------------------------
M13 = blkdiag( diag(m13), diag(m13), diag(m13), diag(m13) );

A=A-By-Dx-M13*B0;
B=-B+2*M13*C1;
D=-D;

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
LHS(2*N, N*2) = 1;
RHS(2*N) = 0;

% v(0)=0
LHS(3*N,:) = Ov;
LHS(3*N, 3*N) = 1;
RHS(3*N) = 0;

% w(0)=0
LHS(4*N,:) = Ov;
LHS(4*N, 4*N) = 1;
RHS(4*N) = 0;

% u(ymax)=0
LHS(N+1,:) = Ov;
LHS(N+1,N+1) = 1;
RHS(N+1) = 0;

% v(ymax)=0
LHS(2*N+1,:) = Ov;
%LHS(2*N+1, 2*N+1) = 1;
LHS(2*N+1, 2*N+1) = U.v(N);
LHS(2*N+1, 1) = 1;
RHS(2*N+1) = 0;

% w(ymax)=0
LHS(3*N+1,:) = Ov;
LHS(3*N+1, 3*N+1) = 1;
RHS(3*N+1) = 0;

end
