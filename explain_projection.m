clear all
close all

addpath('./matlab_script')
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m,n,p,[0.12 0.08],[0.15 0.08], [0.08 0.04]);
if ~make_it_tight, clear subplot;end
casen = 1;

switch casen
    case 1
        L = load('dominant_modes_05_sL_hann.mat');
        ininame = 'opt_mode';
         ylimp200 = [0,0.25];
        ylimp300 = [0,2e-10];
    case 2
        L = load('dominant_modes_05_lL_han_N34.mat');
        ininame = 'optw_mode';
        ylimp200 = [0,0.2];
        ylimp300 = [0,6e-10];
    case 3
        L = load('dominant_modes_3_sL_hann.mat');
        ininame = 'opt_mode';
         ylimp200 = [0,0.28];
        ylimp300 = [0,6e-9];
        
    case 4
        L = load('dominant_modes_3_lL_han.mat');
        ininame = 'optw_mode';
         ylimp200 = [0,0.3];
        ylimp300 = [0,3e-8];
end
uh = L.uh;
vh = L.vh;
wh = L.wh;
x = L.X(:,1,1);
n = L.N(1,:,1);
ft = L.ft*2*pi;
fz = L.fz*2*pi;
modes = L.mo;


dth = L.dth; 
xdth = L.xa;

[x0dth, inx0] = min(abs(xdth));
xdth = xdth(inx0:end);
dth = real(dth(inx0:end));

% dthin = interp1(xdth,real(dth),x,'pchip');
% figure(200)
% XX = squeeze(L.X(:,:,1));
% NN = squeeze(L.N(:,:,1));
% Var = abs(uh{1});
% pcolor(XX,NN,Var)
% hold on
% plot(x,3*dthin,'k--')
% plot(xdth,3*dth,'r--')
% view(2)
% shading interp


ny = length(n);
nx = length(x);

colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];


%ki = [length(xfvec),3,5,2,11,3,length(xfvec),length(xfvec)];
clear lableg
fig1 = figure(100);
fig1.Position = [500 500 700 400];
fig2 = figure(200);
fig2.Position = [500 500 300 400];
fig3 = figure(300);
fig3.Position = [500 500 300 400];
fig4 = figure(400);
fig4.Position = [500 500 300 400];
hold on

fig5 = figure(500);
fig5.Position = [500 500 1200 400];

i=3;

O = load([ininame,num2str(i),'_difx.mat']);

xop = O.xw{1};
xff = xop(1);
xi = find(x>=xff,1,'first');
xfvec = linspace(xop(1)*1.5,0.33,30);

% DNS

ui = squeeze(uh{i}(xi,:));
vi = squeeze(vh{i}(xi,:));
wi = squeeze(wh{i}(xi,:));

uenvo = 0*O.xw{1};
Edifv = zeros(size(xfvec));

count = 0;
count2 = 0;

for k = 1:length(xfvec)%[10,20]
count = count +1;
% Opt
[nn2,Nstations2] = size(O.q{k});
N = nn2/4;
yop = O.y_opt{k};
q0 = O.q{k}(:,1);
uop = (q0(1:N,1));
wop = (q0(2*N+1:N*3,1));
vop = (q0(N+1:N*2,1));

% interp Op
uo = interp1(yop,uop,n);
vo = interp1(yop,vop,n);
wo = interp1(yop,wop,n);


aint =  conj(vi).*vo + conj(wi).*wo;
a = 0.5*trapz(n,aint);

Eo =0.5*(abs(vo).^2+abs(wo).^2);

scl = abs(a);
%Ei = 0.5*(abs(ui).^2);
Ei = 0.5*(1*abs(ui).^2+0*(abs(vi)).^2+0*(abs(wi)).^2);



[Emax,nEpmax] = max(Eo);
nEp0 = find(Eo(nEpmax:end)<1e-5*max(Eo),1,'first');
nEp0 = nEp0 + nEpmax;
dthi = interp1(xdth,dth,xff);
nind = find(n>=3*dthi,1,'first');
%nEp0 =nind;
errorE(k) = abs(a).^2/trapz(n(1:nEp0),Ei(1:nEp0));
qmax2 = zeros(Nstations2,1);
for j=1:Nstations2
    dthi = interp1(xdth,dth,O.xw{k}(j));
    nind = find(flip(O.y_opt{k})>=6*dthi,1,'first');
    %qmax2(j) = max(abs(O.q{k}(length(O.y_opt2{k})-nind:end,j)));
    qmax2(j) = max(abs(O.q{k}(1:N,j)));
end


Edifv(k) = abs(a).^2;

for j=1:length(uenvo)
    if qmax2(j)*scl>uenvo(j)
        uenvo(j) = qmax2(j)*scl;
    end
end

figure(100)
hold on
plot(O.xw{k},qmax2*scl,'-','Color',[0.5 0.5 0.5])

xio = find(O.xw{k}>=xfvec(k),1,'first');
plot(O.xw{k}(xio),qmax2(xio)*scl,'*','Color','k')


umax = zeros(nx,1);

for j=1:nx
    dthi = interp1(xdth,dth,x(j));
    if isnan(dthi)
        dthi=0;
    end
    nind = find(n>=3*dthi,1,'first');
    umax(j) = sqrt(1)*max(abs(squeeze(uh{i}(j,1:nind))));
end



figure(100)
plot(x,umax-0*umax(xi),'Color',colr(i,:),'LineWidth',1.5)
%plot(xenv,uenv,'--','Color',colr(count,:),'LineWidth',1.5)
if count==1
%plot(O.xw{k}(1),0,'.','MarkerSize',15,'Color',[0.5,0.5,0.5])
%text(O.xw{k}(1)*1.2,0.0001,'$x_0$','Interpreter','latex','Color',[0.5,0.5,0.5],...
%    'FontSize',18)
end
%plot(O.xw{k}(xio),0,'.','MarkerSize',15,'Color',[0.5,0.5,0.5]);
%text(O.xw{k}(xio)*1.01,0.0001,'$x_f$','Interpreter','latex','Color',[0.5,0.5,0.5],...
%    'FontSize',18)
box on
grid on
xlim([0,0.33])
%plot(x,umax-1*umax(xi),'r')
ylabel('$|\hat{u}|_{max}$','Interpreter','latex','FontSize',22)
xlabel('$x$','Interpreter','latex','FontSize',22)
if count==1
figure(200)
hold on
plot(abs(uo),n,'r-','DisplayName','$|\hat{u}|$')
plot(abs(vo),n,'b-','DisplayName','$|\hat{v}|$')
plot(abs(wo),n,'k-','DisplayName','$|\hat{w}|$')
xlabel('$|\hat{u}|,|\hat{v}|,|\hat{w}|$','Interpreter','latex','FontSize',22)
ylabel('$n$','Interpreter','latex','FontSize',22)
title('Optimal','Interpreter','latex','FontSize',20)
%xlim([0,0.33])
ylim([0,5e-3])
box on
grid on
legend('Interpreter','latex','FontSize',16)
figure(300)
hold on
plot(abs(uo)*abs(a),n,'r-','DisplayName','$|\hat{u}|$')
plot(abs(vo)*abs(a),n,'b-','DisplayName','$|\hat{v}|$')
plot(abs(wo)*abs(a),n,'k-','DisplayName','$|\hat{w}|$')
xlabel('$|\hat{u}|,|\hat{v}|,|\hat{w}|$','Interpreter','latex','FontSize',22)
ylabel('$n$','Interpreter','latex','FontSize',22)
title('Optimal$\cdot |a|$','Interpreter','latex','FontSize',20)
xlim([0,1.2e-3])
ylim([0,5e-3])
box on
grid on
legend('Interpreter','latex','FontSize',16)

figure(400)
hold on
plot(abs(ui),n,'r-','DisplayName','$|\hat{u}|$')
plot(abs(vi),n,'b-','DisplayName','$|\hat{v}|$')
plot(abs(wi),n,'k-','DisplayName','$|\hat{w}|$')
xlabel('$|\hat{u}|,|\hat{v}|,|\hat{w}|$','Interpreter','latex','FontSize',22)
ylabel('$n$','Interpreter','latex','FontSize',22)
title('DNS','Interpreter','latex','FontSize',20)
xlim([0,1.2e-3])
ylim([0,5e-3])
box on
grid on
legend('Interpreter','latex','FontSize',16)



end

if k==5 || k==15 || k==20 || k==25
count2 = count2+1;
    figure(500)
    subplot(1,4,count2)
    hold on
    xidns = find(x>xfvec(k),1,'first');
    ufdns = squeeze(uh{i}(xidns,:));
    pl1 =plot(abs(ufdns),n,'Color',colr(i,:),'DisplayName','DNS');
    pl2 = plot(abs(O.q{k}(1:N,xio))*scl,yop,'Color',[0.5,0.5,0.5],'DisplayName','Optimal');
    [umaxo,indmax] = max(abs(O.q{k}(1:N,xio)));
    plot(umaxo*scl,yop(indmax),'k*')
    xlabel('$|\hat{u}|$','Interpreter','latex','FontSize',22)
    ylabel('$n$','Interpreter','latex','FontSize',22)
    title(['$x=',num2str(xfvec(k),'%1.2f'),'$'],'Interpreter','latex','FontSize',20)
    ylim([0,5e-3])
    box on
grid on
legend([pl1,pl2],'Interpreter','latex')
end
end 
%yline(1,'--')










