clear all
casen = 3;

switch casen
    case 1
        L = load('dominant_modes_05_sL_han_N.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.01+0*[0.01,0.015,0.01,0.005,0.01,0.01,0.0075,0.015]; % small L
        ymaxp = 3e-3;
        maxvp=1e-3;
        colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    ];
        indm = [1:8];
    case 2
        L = load('dominant_modes_05_lL_han_N.mat');
        ydm = 0.045;
        ymid = 5e-3;
        xiv =0.02+0*[0.05,0.04,0.02,0.01,0.05,0.025,0.02,0.015,0.02,0.02,0.02,0.02]; % large L
        ymaxp = 5e-3;
        maxvp=1.2e-3;
        colr=[[0.4660,0.6740,0.1880];[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.310,0.745,0.933] ;...
    [0.4660,0.6740,0.1880];[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.310,0.745,0.933]...
    ];

        indm = [1:5,7:11];
        %indm = [6,12];
        colr = colr(indm,:);
    case 3
        L = load('dominant_modes_3_sL_han_N.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.01+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
        ymaxp = 16.5e-3;
        maxvp=2.2e-3;
        colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    ];
        indm = [1:8];
    case 4
        L = load('dominant_modes_3_lL_han_N.mat');
        ydm = 0.045;
        ymid = 5e-3;
        xiv =0.02+0*[0.05,0.04,0.02,0.01,0.05,0.025,0.02,0.015,0.02,0.02]; % large L
        ymaxp = 18e-3;
        maxvp=2.2e-3;
         colr=[[0.4660,0.6740,0.1880];[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560] ;...
    [0.4660,0.6740,0.1880];[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    ];

        indm = [1:10];
    case 5
        L = load('dominant_modes_05_sL_N2.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.01+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
        ymaxp = 3e-3;
        maxvp=1e-3;
    case 6
        L = load('dominant_modes_05_sL_N34.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.01+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
        ymaxp = 3e-3;
        maxvp=1e-3;
end
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m,n,p,[0.12 0.04],[0.1 0.05], [0.05 0.02]);
if ~make_it_tight, clear subplot;end

uh = L.uh(indm);
vh = L.vh(indm);
wh = L.wh(indm);
x = L.X(:,1,1);
n = L.N(1,:,1);
ft = L.ft*2*pi;
fz = L.fz*2*pi;
modes = L.mo(indm);
dth = L.dth; 
xdth = L.xa;
[x0dth, inx0] = min(abs(xdth));
xdth = xdth(inx0:end);
dth = real(dth(inx0:end));

Nmo = length(modes);

% Bending function


ystar = 1-(n-ydm)/(max(n)-ydm);
S = 1./(1+exp(1./(ystar-1) +1./ystar));
S(ystar<=0) = 0;
S(ystar>=1)=1;
fig1 = figure(101);
fig1.Position = [500 500 1500 500];
%%
count =0;
for i = 1:Nmo
    count= count +1;
    xind = find(x>=xiv(i),1,'first');
xi = x(xind);
N = length(n);
q0=zeros(4*N,1);
u0 = squeeze(uh{i}(xind,:)).*S*1;
v0 = squeeze(vh{i}(xind,:)).*S*1;
w0 = squeeze(wh{i}(xind,:)).*S*1;
q0(1:N,1) = flip(u0);
q0(N+1:2*N,1) = flip(v0);
q0(2*N+1:3*N,1) = flip(w0);
mi = modes{i};
beta = fz(mi(1));
freq = -ft(mi(2))/(2*pi);
[xw,q,y_opt,En,bbb2,xarc2] = streak_readBL_eval_mode(0.34,xi,beta,freq,N,q0,max(n),ymid);

q0(1:N,1) = flip(u0)*0;
q0(N+1:2*N,1) = flip(v0);
q0(2*N+1:3*N,1) = flip(w0);
[xw2,q2,y_opt2,En,bbb2,xarc2] = streak_readBL_eval_mode(0.34,xi,beta,freq,N,q0,max(n),ymid);


%%

umax = zeros(size(x)-1);


for j=1:length(x)-1
    dthj = interp1(xdth,dth,x(j+1));
    nind = find(n>=3*dthj,1,'first');
    umax(j) = max(abs(uh{i}(j+1,1:nind)));
end


[nn,Nstations] = size(q);
maxop = zeros(Nstations,1);
maxop2 = maxop;
for j=1:Nstations
    dthj = interp1(xdth,dth,xw(j));
    nind = find(n>=3*dthj,1,'first');
    maxop(j) = max(abs(q(N-nind:N,j)));
    maxop2(j) = max(abs(q2(N-nind:N,j)));
end
%maxop = max(abs(q(1:N,:)),[],1);
%maxop2 = max(abs(q2(1:N,:)),[],1);

%%

figure(101)
subplot(2,Nmo/2,count)
hold on
plot(x(2:end),umax,'Color',colr(i,:),'LineWidth',1.5)
plot(xw,maxop,'Color','k');%[0.5,0.5,0.5])
plot(xw,maxop2,'--','Color','k')%[0.5,0.5,0.5])
if i>Nmo/2
xlabel('$x$','Interpreter','latex','FontSize',22)
end
if i==1 || i==Nmo/2+1
ylabel('$|\hat{u}|_{max}$','Interpreter','latex','FontSize',22)
end
box on
grid on
title(num2str(i))
xlim([0,0.33])
ylim([0,ymaxp])


yop = y_opt;
figure(201)
subplot(2,Nmo/2,count)
hold on
plot(abs(q(1:N,1)),yop,'r')
plot(abs(q(N+1:2*N,1)),yop,'b')
plot(abs(q(2*N+1:3*N,1)),yop,'k')
legend('|u|','|v|','|w|')
xlabel('$|u|,|v|,|w|$','Interpreter','latex','FontSize',16)
ylabel('$n$','Interpreter','latex','FontSize',16)
box on
grid on
ylim([0,0.002])
xlim([0,maxvp])
title(num2str(i))

figure(301)
subplot(2,Nmo/2,count)
hold on
plot(wrapTo2Pi(angle(q(1:N,1))),yop,'r')
plot(wrapTo2Pi(angle(q(N+1:2*N,1))),yop,'b')
plot(wrapTo2Pi(angle(q(2*N+1:3*N,1))),yop,'k')
legend('|u|','|v|','|w|')
xlabel('$phase$','Interpreter','latex','FontSize',16)
ylabel('$n$','Interpreter','latex','FontSize',16)
box on
grid on
ylim([0,0.01])
%xlim([0,1.0]*1e-3)
title(num2str(i))

end
%%

% yop = y_opt;
% figure()
% hold on
% plot(abs(q(1:N,1)),yop)
% plot(abs(q(N+1:2*N,1)),yop)
% plot(abs(q(2*N+1:3*N,1)),yop)
% legend('|u|','|v|','|w|')
% box on
% grid on
% ylim([0,0.01])
% figure()
% hold on
% plot(angle(q(1:N,1)),yop)
% plot(angle(q(N+1:2*N,1)),yop)
% plot(angle(q(2*N+1:3*N,1)),yop)
% legend('|u|','|v|','|w|')
% box on
% grid on
% ylim([0,0.01])
% 
% [uhmax,inn] = max(abs(u0));
% nmax = n(inn);
% nratio =nmax./real(dth);
% nratio(isnan(nratio))=0;
% figure()
% plot(xdth,nratio)
% xlim([0.005,0.35])
% 
% figure()
% plot(xdth,dth)
% xlim([0.005,0.35])




