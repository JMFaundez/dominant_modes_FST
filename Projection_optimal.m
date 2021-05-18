clear all
close all
addpath('./matlab_script')
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m,n,p,[0.12 0.08],[0.15 0.08], [0.08 0.04]);
if ~make_it_tight, clear subplot;end
casen = 2;

switch casen
    case 1
        L = load('dominant_modes_05_sL_hann.mat');
        ininame = 'opt2_mode';
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
%%

colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];


%ki = [length(xfvec),3,5,2,11,3,length(xfvec),length(xfvec)];
clear lableg
fig1 = figure(100);
fig1.Position = [500 500 1500 500];
fig2 = figure(200);
fig2.Position = [500 500 900 300];
fig3 = figure(300);
fig3.Position = [500 500 900 300];
fig4 = figure(400);
fig4.Position = [500 500 900 400];
hold on
count = 0;

for i=1:8
    count = count+1;



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

for k=1:length(xfvec)

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
subplot(2,4,count)
hold on
plot(O.xw{k},qmax2*scl,'-','Color',[0.5 0.5 0.5])

xio = find(O.xw{k}>=xfvec(k),1,'first');
plot(xfvec(k),qmax2(xio)*scl,'*','Color','k')

if k==1
    xenv = O.xw{k}(1:xio);
    uenv = qmax2(1:xio)*scl;
else
    xenv = [xenv,O.xw{k}(xio)];
    uenv = [uenv;qmax2(xio)*scl];
end

%if k==ki(i)%length(xfvec)
    %plot(O.xw{k},qmax2*scl,'--','Color',colr(count,:),'LineWidth',2)
    %figure(4)
    %subplot(1,2,2)
    %hold on
    %plot(O.xw{k},qmax2*scl,'--','Color',colr(count,:),'LineWidth',1.5)
    %xlabel('$x$','Interpreter','latex','FontSize',16)
    %ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
    %xlim([0,0.33])
    %ylim([0,3e-3])
    %box on
    %grid on
%end

end


umax = zeros(nx,1);

for j=1:nx
    dthi = interp1(xdth,dth,x(j));
    if isnan(dthi)
        dthi=0;
    end
    nind = find(n>=3*dthi,1,'first');
    umax(j) = sqrt(1)*max(abs(squeeze(uh{i}(j,1:nind))));
end
if i<=4
    ipp = 1;
else
    ipp = 2;
end

figure(200)

if ipp==1
    subplot(1,2,1)
    semilogy(xfvec,Edifv,'-o','DisplayName',num2str(i),'Color',colr(count,:))
else
    subplot(1,2,2)
    semilogy(xfvec,Edifv,'-o','DisplayName',num2str(i),'Color',colr(count,:))
end
hold on
box on
grid on
xlabel('$x_f$','Interpreter','latex','FontSize',14)
ylabel('$|a|^2$','Interpreter','latex','FontSize',14)
legend()
ylim(ylimp300)
xlim([0,0.33])

figure(100)
subplot(2,4,count)
plot(x,umax-0*umax(xi),'Color',colr(count,:),'LineWidth',1.5)
%plot(xenv,uenv,'--','Color',colr(count,:),'LineWidth',1.5)
box on
grid on
xlim([0,0.33])
%plot(x,umax-1*umax(xi),'r')
title(num2str(i))
xlabel('$x$','Interpreter','latex','FontSize',14)
ylabel('$u_{max}$','Interpreter','latex','FontSize',14)



figure(300)


if ipp==1
    subplot(1,2,1)
 plot(xfvec,errorE,'-o','Color',colr(count,:),'DisplayName',num2str(i))
else
    subplot(1,2,2)
    plot(xfvec,errorE,'-o','Color',colr(count,:),'DisplayName',num2str(i))
end
hold on
xlabel('$x_f$','Interpreter','latex','FontSize',14)
ylabel('$|a|^2/E_{dns}$','Interpreter','latex','FontSize',14)
box on
grid on
legend()
ylim(ylimp200)
xlim([0,0.33])
%plot(xfvec,Edifv,'<-','Color',colr(count,:))


figure(400)
subplot(1,2,ipp)
%subplot(1,1,1)
hold on
l1 = plot(x,umax-0*umax(xi),'Color',colr(count,:),'LineWidth',1.5,'DisplayName',num2str(i));
plot(xenv,uenv,'k--')%,'Color',colr(count,:),'LineWidth',1.5)
lableg(count) = l1;
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
xlim([0,0.33])
%ylim([0,3e-3])
box on
grid on
legend(lableg)
%yline(1,'--')

end









