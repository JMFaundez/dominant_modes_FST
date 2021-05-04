function [xw,q,yp,energy,beta,x] = streak_readBL(xf,x0,beta,F,if_plot,N,q0 )
addpath('./streak_script')
%
% streak.m is the main program for solving the Linearized Boundary Laayer equations
%
% Input and output energy norm correspond to those in Andersson et al., 1999.
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%
% Note on this version: The global length scale is L=sqrt(x_1*nu/Ue_1),
% where subscript '1' denotes location where the energy is optimised.
%
%clear all
%% --------------------------------------------------------------------------
% Flow parameters
Re=5.33e5;    % Reynolds number 
L = load('meanflow_BL.mat');
Uin = L.U; % u,v,w,ux,vx,wx,uy,vy,wy (u: tangent vel, v: normal vel)
X = L.X;
x = X.x;    % arclength along surface
xc= X.xc;   % x/c positions
L2 = load('Base_Flow.mat');
%Xw = L2.Xw; % Airfoil coordinates
%xw = Xw.x';
%yw = Xw.y';
yin = L2.X.y;


A = 0.2969;                                                                                       
B = 0.1260;
C = 0.3516;
D = 0.2843;                                                                                       
E = 0.1015;                                                                                       
naca = @(x) 5*0.08*(A*sqrt(x) - B*x - C*x.^2 + D*x.^3 - E*x.^4);                                  
dnaca = @(x) 5*0.08*(0.5*A*1./sqrt(x) - B - 2*C*x + 3*D*x.^2 - 4*E*x.^3);      
xw = xc;
yw = naca(xw);

%% --------------------------------------------------------------------------
% Numerical Parameters
%x0=0.003;    % x_0
%xf=0.34;   % x at final station

x0_ind=find(xw>=x0,1,'first');
xf_ind=find(xw<=xf,1,'last');

dstar=trapz(Uin(xf_ind).y,1-Uin(xf_ind).u/Uin(xf_ind).u(end)); % dispalcement thickness at last station

%N=150;             % Number of nodes in wall-normal direction
%ymax=40*dstar;     % Height of last node in wall-normal direction
FDorder='first';   % order of backward Euler discretization scheme (can be 'first' or 'second')

%% --------------------------------------------------------------------------
% Spectral discretisation operators (Chebycheff)

%[y,D1,D2,W] = chebmat(N,ymax);

if beta>=1100
    ymax = 1*max(yin);
    ymid = 0.3*ymax;
elseif beta<=100
    ymax = 12*max(yin);
    ymid = 0.2*ymax;
else
    ymax = 3*max(yin);
    ymid = 0.3*ymax;
end
ymax = 0.081;
if abs(beta)<100
    ymax = 0.22;
end
ymid = 0.2*ymax;
[y,D1,D2,W] = chebmat_trans(N,ymax,ymid);
yp = y;
%% --------------------------------------------------------------------------
% Interpolate meanflow

disp('Interpolate mean flow')
U = MeanFlow_BL(y,Uin);


%U(1).vx=(U(2).v-U(1).v)/(x(2)-x(1));
%for i = 2:length(U)
%    U(i).vx=(U(i).v-U(i-1).v)/(x(i)-x(i-1));
%end

% Metrics
%xw=x; yw=xw*0;     % should be replaced with geometry of the airfoil
curv=curvature(xw,yw);
metric=struct;
for k=1:length(U)
    h1=1+y*curv(k);
    m13=curv(k)./h1;
    metric(k).h1=h1;
    metric(k).m13=m13;

    U(k).ux=U(k).ux./h1;
    U(k).vx=U(k).vx./h1;
    U(k).wx=U(k).wx./h1;
end
%% --------------------------------------------------------------------------
x=x(x0_ind:xf_ind);
xw=xw(x0_ind:xf_ind);
U=U(x0_ind:xf_ind);
metric=metric(x0_ind:xf_ind);
Nstation=length(x);       % Number of x-stations

disp(['x0= ',num2str(x(1)),'   xf= ',num2str(x(end))]);


%% --------------------------------------------------------------------------
% Perturbation parameters
% if if_blasius
% 	%beta_blasius=0.45; % Spanwise wavenumber based on the Blasius lengthscale at xf
% 	beta_scale=sqrt(Re*U(end).u(1)/x(end));
% 	beta= beta_blasius*beta_scale;  % Spanwise wavenumber based on the lengthscale used in the meanflow
% end
%F=0;           % Reduced frequency
omega=F*2*pi;    % Angular frequency


%% --------------------------------------------------------------------------
% initial guess
%q0=zeros(4*N,1)+1;
E0=0.5*q0(1:3*N)'*blkdiag(W,W,W)*q0(1:3*N);
q0=q0/sqrt(E0);

[q] = integrateLBL(q0, x, y, U, metric, D1, D2, beta, omega, Re, Nstation, FDorder);

Ef0=0.5*q(1:3*N,Nstation)'*blkdiag(W,W,W)*q(1:3*N,Nstation); % energ at xf
dE=1.0;
fprintf('\n\n E/E0= %f\n',Ef0);

%% --------------------------------------------------------------------------
% Integrate state and adjoint equations untill convergence reached

while dE>1e-3
    % Initial condition for adjoint quantities
    qf=q(:,Nstation);
    qadjN=[qf(1:N); qf(1:N)*0; qf(N+1:2*N)*0; qf(2*N+1:3*N)*0];
    
    [qadj] = integrateLBL_adj(qadjN, x, y, U, metric, D1, D2, beta, omega, Re, Nstation, FDorder);
    
    % Initial condition for state quantities
    U1=U(1); qadj1=qadj(:,1);
    
    q0 = [qadj1(N+1:2*N).*U1.u*0;  qadj1(2*N+1:3*N).*U1.u; qadj1(3*N+1:4*N).*U1.u; qadj1(1:N)*0;];
    
    E0 = 0.5*q0(1:3*N)'*blkdiag(W,W,W)*q0(1:3*N);
    q0 = q0/sqrt(E0); % normalize the energy at x=x_0
    
    [q] = integrateLBL(q0, x, y, U, metric, D1, D2, beta, omega, Re, Nstation, FDorder);
    
    Ef=abs(0.5*q(1:3*N,Nstation)'*blkdiag(W,W,W)*q(1:3*N,Nstation)); % energ at xf
    
    dE=abs(Ef-Ef0)/Ef0;
    Ef0=Ef;
    fprintf('\n\n E/E0= %f\t dE/E= %f\n',Ef,dE);
end

%% --------------------------------------------------------------------------
% Compute amplitude
energy=zeros(size(x));
for k=1:length(x)
    energy(k)=abs(0.5*q(1:3*N,k)'*blkdiag(W,W,W)*q(1:3*N,k));
end

if if_plot
	figure(111)
	plot(x,energy)
	ylabel('E/E_0')
	xlabel('x')
	q = q;
	%%
	% Plot amplitude functions
	figure(222)
	plot(abs(q(1:N,1)),y,abs(q(N+1:N*2,1)),y,'r-',abs(q(2*N+1:N*3,1)),y,'r--')
	legend('abs(u)','abs(v)','abs(w)')
	title('Optimal initial perturbations')
	
	figure(333)
	plot(abs(q(1:N,end)),y,abs(q(N+1:N*2,end)),y,'r-',abs(q(2*N+1:N*3,end)),y,'r--')
	legend('abs(u)','abs(v)','abs(w)')
	title('optimal response')
end
%%

end
