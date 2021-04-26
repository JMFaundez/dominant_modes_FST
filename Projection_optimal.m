clear all
%% Base Flow
bflL = to_ut('data_sL/data_BF.mat'); 

%% Data FST
lL05 = to_ut('data_sL/data_05p_sL.mat'); 
%% Commpute pert

%u05 = lL05.ut- bflL.ut;
u05 = lL05.ut(:,:,:,1:end);%- bflL.ut;
u05 = u05 - mean(mean(u05,1),4);

%v05 = lL05.un - bflL.un;
v05 = lL05.un(:,:,:,1:end);% - bflL.un;
v05 = v05 - mean(mean(v05,1),4);

%w05 = lL05.uz - bflL.uz;
w05 = lL05.uz(:,:,:,1:end);% - bflL.uz;
w05 = w05 - mean(mean(w05,1),4);

[nt,nx,ny,nz] = size(u05);
zpad = nt;
uhat = fft(u05,[],4);
uhat = fft(uhat,zpad,1);
fftut = uhat/(zpad*nz);
vhat = fft(v05,[],4);
vhat = fft(vhat,zpad,1);
fftun = vhat/(zpad*nz);
what = fft(w05,[],4);
what = fft(what,zpad,1);
fftuz = what/(zpad*nz);

uhat = abs(uhat);
vhat = abs(vhat);
what = abs(what);

[ft,ftp] = freq_fft(zpad,lL05.time(end));
[fz,fzp] = freq_fft(nz,0.02-0.02/nz);
[W,B] = meshgrid(fftshift(ft),fftshift(fz));
W = W*2*pi; B = B*2*pi;


%%
base = base_case('fringe_m90.f00008',300,22);
dth = base.dth;
xa = base.xx(1,:);

%%
xp = 0.15;
x = lL05.xx(:,1,1);
n = lL05.nn(1,:,1);
xind = find(x>=xp,1,'first');
dthi = interp1(xa,dth,x(xind));
nind = find(n>=1.5*dthi,1,'first');
U = abs(squeeze(uhat(:,xind,nind,:)));
figure()
pcolor(W,B,fftshift(U)')
view(2)
shading interp
xlim([0,40])
ylim(([-1000,1000]))


%%
mo = {[2,6],[2,9],[3,12],[4,15],[nz,6],[nz,9],[nz-1,12],[nz-2,15]};
uh = {}; vh={};wh={};
for i=1:length(mo)
    uh{i} = (squeeze(fftut(mo{i}(2),:,:,mo{i}(1)))+0/(nz*nt));
    vh{i} = (squeeze(fftun(mo{i}(2),:,:,mo{i}(1)))+0/(nz*nt));
    wh{i} = (squeeze(fftuz(mo{i}(2),:,:,mo{i}(1)))+0/(nz*nt));
end
X = lL05.xx;
N = lL05.nn;

save('dominant_modes_05_sL.mat','mo','uh','vh','wh','X','N','fz','ft','X','N','xa','dth')



%%
[usts, xsts,nu, data] = urmsmax('data_sL/stsINT05_sL.mat');


%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m,n,p,[0.1 0.08],[0.15 0.05], [0.1 0.05]);
if ~make_it_tight, clear subplot;end

mo = {[2,5],[2,7],[3,9],[4,11],[nz,5],[nz,7],[nz-1,9],[nz-2,11]};

i=3;
O = load(['Opt_data/opt_m',num2str(i),'.mat']);
x = fut.xx(:,1,1);
n = fut.nn(1,:,1);
xff = 0.005;
xi = find(x>=xff,1,'first');

% DNS

% vi = abs(2*squeeze(fftun(mo{i}(2),xi,:,mo{i}(1))));
% wi = abs(2*squeeze(fftuz(mo{i}(2),xi,:,mo{i}(1))));
% ui = abs(2*squeeze(fftut(mo{i}(2),xi,:,mo{i}(1))));

ui = squeeze(fftut(mo{i}(2),xi,:,mo{i}(1)))/(nz*nt);
vi = squeeze(fftun(mo{i}(2),xi,:,mo{i}(1)))/(nz*nt);
wi = squeeze(fftuz(mo{i}(2),xi,:,mo{i}(1)))/(nz*nt);

% Opt
[nn2,Nstations2] = size(O.q);
N = nn2/4;
yop = O.y_opt2;
q0 = O.q(:,1);
uop = (q0(1:N,1));
wop = (q0(2*N+1:N*3,1));
vop = (q0(N+1:N*2,1));
%uop = (q0(1:N,1));
%wop = (q0(2*N+1:N*3,1));
%vop = (q0(N+1:N*2,1));

% interp Op
uo = interp1(yop,uop,n);
vo = interp1(yop,vop,n);
wo = interp1(yop,wop,n);


% project over optimal
up = uo*0;
vp = vo*0;
wp = wo*0;
% A = [ui,vi,wi];
% B = [uo',vo',wo'];
% prod = dot(A,A,2).^0.5;
% C = sqrt(sum(cross(A,B,2).^2,2));
% dotAB = dot(A,B,2);
% angle = atan2(C,dotAB);
% prod = prod.*cos(angle);
% uabs = sqrt(wo.^2 + vo.^2);
% 
% prod = dotAB;
% bnorm = uabs;
% bnorm(bnorm<=1e-2) = 1;
% prod = dotAB./bnorm';
% 
% 
% vp = prod'.*(vo./uabs);
% wp = prod'.*(wo./uabs);
tol1 = 1e-1;
%vo(vo<tol1)=0;
%wo(wo<tol1)=0;
tol = 10e10;
Bmax=0;
for j=1:ny
    A =[ vi(j),wi(j)];
    B = [vo(j), wo(j)];
    %prod = dot(A,B);
    %prod = dot(A,A).^0.5;
    %C = sqrt(sum(cross(A,B).^2));
    %dotAB = dot(A,B,2);
    %angle = atan2(C,dotAB);
    Bnorm = sqrt(sum(B.^2));
    uabs = Bnorm;
    %if Bnorm<tol
    %    Bnorm =1;
    %end
    %prod = prod*cos(angle);
    prod = dot(A,B);
    if prod>=tol
        prod = prod/Bnorm;
    end
    if Bnorm>Bmax
        Bmax = Bnorm;
    end
    vp(j) = prod*(vo(j)/uabs); 
    wp(j)= prod*(wo(j)/uabs);
end

vp = vp/Bmax;
wp = wp/Bmax;


figure('Position',[500 500 1000, 400])
subplot(1,3,1)
hold on
plot(vi,n)
plot(wi,n)
xlabel('$|v|,|w|$','Interpreter','latex','FontSize',16)
ylabel('$y$','Interpreter','latex','FontSize',16)
box on
ylim([0,0.01])
title('DNS')
subplot(1,3,2)
hold on
plot(vo,n)
plot(wo,n)
ylim([0,0.01])
xlabel('$|v|,|w|$','Interpreter','latex','FontSize',16)
ylabel('$y$','Interpreter','latex','FontSize',16)
box on
title('Optimal')
subplot(1,3,3)
hold on
plot(abs(vp),n)
plot(abs(wp),n)
ylim([0,0.01])
xlabel('$|v|,|w|$','Interpreter','latex','FontSize',16)
ylabel('$y$','Interpreter','latex','FontSize',16)
box on
title('Projection')

%%
Ei =0.5*(vp.^2+wp.^2);
Ei(isnan(Ei))=0;
scl = sqrt(trapz(n,Ei));

qmax2 = zeros(Nstations2,1);
for j=1:Nstations2
    qmax2(j) = max(abs(O.q(1:N,j)));
end
umax = zeros(nx,1);
for j=1:nx
   umax(j) = max( 2*abs(squeeze(fftut(mo{i}(2),j,:,mo{i}(1)))));
end

figure()
hold on
plot(O.xw2,qmax2*scl)
plot(x,umax)

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m,n,p,[0.05 0.04],[0.12 0.05], [0.05 0.02]);
if ~make_it_tight, clear subplot;end

mo = {[2,6],[2,9],[3,12],[4,15],[nz,5],[nz,7],[nz-1,9],[nz-2,11]};
xfvec = linspace(0.03,0.33,30);
colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];

errorv = zeros(size(xfvec));
errorw = zeros(size(xfvec));
errorT = zeros(size(xfvec));
errorE = zeros(size(xfvec));

ki = [length(xfvec),3,5,2,11,3,length(xfvec),length(xfvec)];
clear lableg
fig1 = figure(1);
fig1.Position = [500 500 1600 400];
figure(2)
fig4 = figure(4);
fig4.Position = [500 500 900 400];
hold on
count = 0;

qto = {};
qdns = {};
uto = {};
tyo = {};
xto = {};

for i=1:4
    count = count+1;
x = lL05.xx(:,1,1);
n = lL05.nn(1,:,1);
xff = 0.005;
xi = find(x>=xff,1,'first');

% DNS

%vi = abs(sqrt(2)*squeeze(fftun(mo{i}(2),xi,:,mo{i}(1))));
%wi = abs(sqrt(2)*squeeze(fftuz(mo{i}(2),xi,:,mo{i}(1))));
%ui = abs(sqrt(2)*squeeze(fftut(mo{i}(2),xi,:,mo{i}(1))));

ui = sqrt(1)*(squeeze(fftut(mo{i}(2),xi,:,mo{i}(1)))+0/(nz*nt));
vi = sqrt(1)*(squeeze(fftun(mo{i}(2),xi,:,mo{i}(1)))+0/(nz*nt));
wi = sqrt(1)*(squeeze(fftuz(mo{i}(2),xi,:,mo{i}(1)))+0/(nz*nt));

%ui = 2*(squeeze(fftut(mo{i}(2),xi,:,mo{i}(1))+fftut(mo{i+4}(2),xi,:,mo{i+4}(1)))+0/(nz*nt));
%vi = 2*(squeeze(fftun(mo{i}(2),xi,:,mo{i}(1))+fftun(mo{i+4}(2),xi,:,mo{i+4}(1)))+0/(nz*nt));
%wi = 2*(squeeze(fftuz(mo{i}(2),xi,:,mo{i}(1))+fftuz(mo{i+4}(2),xi,:,mo{i+4}(1)))+0/(nz*nt));


% ui = (squeeze(abs(fftut(mo{i}(2),xi,:,mo{i}(1))).^2+abs(fftut(mo{i+4}(2),xi,:,mo{i+4}(1))).^2));
% ui = sqrt(2*ui);
% vi = (squeeze(abs(fftun(mo{i}(2),xi,:,mo{i}(1))).^2+abs(fftun(mo{i+4}(2),xi,:,mo{i+4}(1))).^2));
% vi = sqrt(2*vi);
% wi = (squeeze(abs(fftuz(mo{i}(2),xi,:,mo{i}(1))).^2+abs(fftuz(mo{i+4}(2),xi,:,mo{i+4}(1))).^2));
% wi = sqrt(2*wi);

vk = vi;
wk = wi;
uk = ui;
O = load(['Opt_data/opt_m',num2str(i),'_difx.mat']);
Edifv = zeros(size(xfvec));
ifplot = 1;
vpold = vi*0;
wpold = wi*0;
Vtot = zeros(length(n),2,length(xfvec));

uenvo = 0*O.xw2{1};


for k=1:length(xfvec)

% Opt
[nn2,Nstations2] = size(O.q{k});
N = nn2/4;
yop = O.y_opt2{k};
q0 = O.q{k}(:,1);
uop = (q0(1:N,1));
wop = (q0(2*N+1:N*3,1));
vop = (q0(N+1:N*2,1));

% interp Op
uo = interp1(yop,uop,n);
vo = interp1(yop,vop,n);
wo = interp1(yop,wop,n);


% project over optimal
up = uo*0;
vp = vo*0;
wp = wo*0;

Bmax=0;Amax=0;
for j=1:ny
    A =[ vk(j),wk(j)];
    B = [vo(j), wo(j)];
    Bnorm = sqrt(sum(abs(B).^2));
    uabs = Bnorm;
    Anorm = sqrt(sum(A.^2));
    uabsdns = Anorm;
    prod = dot(A,B);
    if Bnorm>Bmax
        Bmax = Bnorm;
    end
    if Anorm>Amax
        Amax = Anorm;
    end
    vp(j) = prod*(vo(j)/uabs); 
    wp(j)= prod*(wo(j)/uabs);
%    vk(j) = vk(j) - prod*(vk(j)/uabsdns);
%    wk(j) = wk(j) - prod*(wk(j)/uabsdns);
end


vp = vp/Bmax;
wp = wp/Bmax;

%vk = vk -vp;
%wk = wk - wp;

aint =  conj(vi).*transpose(vo) + conj(wi).*transpose(wo);
a = 0.5*trapz(n,aint);

Ep =0.5*(abs(vp).^2+abs(wp).^2);

Eo =0.5*(abs(vo).^2+abs(wo).^2);
Ep(isnan(Ep))=0;
%scl = sqrt(trapz(n,Ep));
scl = abs(a);
Ei = 0.5*((abs(vi)).^2+(abs(wi)).^2);

[Emax,nEpmax] = max(Eo);
nEp0 = find(Eo(nEpmax:end)<1e-5*max(Eo),1,'first') + nEpmax;


errorE(k) = abs(a).^2/trapz(n(nEp0:end),Ei(nEp0:end));
qmax2 = zeros(Nstations2,1);
for j=1:Nstations2
    dthi = interp1(xa,dth,O.xw2{k}(j));
    nind = find(flip(O.y_opt2{k})>=6*dthi,1,'first');
    %qmax2(j) = max(abs(O.q{k}(length(O.y_opt2{k})-nind:end,j)));
    qmax2(j) = max(abs(O.q{k}(1:N,j)));
end
vp(isnan(vp))=0;
wp(isnan(wp))=0;


qto{count,k,1} = ui;
qto{count,k,2} = vp;
qto{count,k,3} = wp;

qdns{count,k,1} = ui;
qdns{count,k,2} = vi;
qdns{count,k,3} = wi;

if k>1
    ampp = sqrt((sum([vp',wp'].^2,2)));
prod = dot([vp',wp'],[vpold',wpold'],2)./ampp;
Edif = 0.5*(ampp-prod).^2;
Vtot(:,1,k) = vp;%- prod'.*vp./sqrt((sum([vp',wp'].^2,2))');
Vtot(:,2,k) = wp;% - prod'.*wp./sqrt((sum([vp',wp'].^2,2))');
else
    Edif = 0.5*(abs(vp).^2+abs(wp).^2);
    Vtot(:,1,k) = vp;
    Vtot(:,2,k) = wp;
end

%vpdif = vp - vpold;
%wpdif = wp - wpold;
%Edif = 0.5*(abs(vpdif).^2+abs(wpdif).^2);
Edif(isnan(Edif)) =0;
Edifv(k) = abs(a).^2;trapz(n,Edif);
vpold = vp;
wpold = wp;
% errorv(k) = immse(vp,wo*scl)./mean((vo*scl).^2);
% errorw(k) = immse(wp,wo*scl)./mean((wo*scl).^2);
%errorv(k) = immse(vp,vo*scl)./mean((vo*scl).^2);
%errorw(k) = immse(wp,wo*scl)./mean((wo*scl).^2);

%errorT(k) = mean((vp-vo*scl).^2 + (wp-wo*scl).^2)/(mean((vp).^2+(wp).^2));
wii = abs(wi);
vii = abs(vi);
vii(nEp0:end) = 0;
wii(nEp0:end) = 0;
errorT(k) = mean(((abs(vp)-transpose(vii)).^2 + (abs(wp)-transpose(wii)).^2))./mean(abs(vii).^2+abs(wii).^2);
%errorT(k) = mean(abs(vp).^2 + abs(wp).^2)/(mean(abs(vi).^2+abs(wi).^2));


% figure(2)
% subplot(i,8,k+(i-1)*8)
% hold on
% plot(abs(vp),n)
% plot(abs(wp),n)
% plot(abs(vo)*scl,n,'k')
% plot(abs(wo)*scl,n,'k--')
% ylim([0,0.01])
% title(['x_f=',num2str(xfvec(k))])
if sum(Edifv)/trapz(n,Ei)>=1
    ifplot = 1;
end


for j=1:length(uenvo)
    if qmax2(j)*scl>uenvo(j)
        uenvo(j) = qmax2(j)*scl;
    end
end

if ifplot
figure(1)
subplot(1,4,count)
hold on
plot(O.xw2{k},qmax2*scl,'-','Color',[0.5 0.5 0.5])

xio = find(O.xw2{k}>=xfvec(k),1,'first');
plot(xfvec(k),qmax2(xio)*scl,'k*')

if k==ki(i)%length(xfvec)
    plot(O.xw2{k},qmax2*scl,'--','Color',colr(count,:),'LineWidth',2)
    figure(4)
    subplot(1,2,2)
    hold on
    plot(O.xw2{k},qmax2*scl,'-','Color',colr(count,:),'LineWidth',1.5)
    xlabel('$x$','Interpreter','latex','FontSize',16)
    ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
    xlim([0,0.33])
    ylim([0,4.0e-3])
    box on
    grid on
end
 end


end

for k=1:length(xfvec)
    vinit = squeeze(Vtot(:,:,end));
    vf =  squeeze(Vtot(:,:,k));
    ampp = sqrt((sum(vinit.^2,2)));
    prod = dot(vinit,vf,2)./ampp;
    Edif = 0.5*(prod).^2;
    Edif(isnan(Edif)) =0;
    Edifv(k) =trapz(n,Edif)/trapz(n,0.5*ampp.^2);
end
% 
% figure(5)
% subplot(4,2,1+(count-1)*2)
% plot(n,squeeze(abs(Vtot(:,1,:))),'k-')
% xlim([0,0.015])
% subplot(4,2,2+(count-1)*2)
% plot(n,squeeze(abs(Vtot(:,2,:))),'k--')
% xlim([0,0.015])

umax = zeros(nx,1);

for j=1:nx
    dthi = interp1(xa,dth,x(j));
    if isnan(dthi)
        dthi=0;
    end
    nind = find(n>=3*dthi,1,'first');
    umax(j) = sqrt(1)*max(abs(squeeze(fftut(mo{i}(2),j,1:nind,mo{i}(1)))))+0/(nz*nt);
    %umaxrms = (squeeze(abs(fftut(mo{i}(2),j,1:nind,mo{i}(1))).^2+abs(fftut(mo{i+4}(2),j,1:nind,mo{i+4}(1))).^2));
    %umax(j) = max(sqrt(2*umaxrms));
    %  umaxTrms = abs(squeeze(fftut(mo{i}(2),j,1:nind,mo{i}(1)))).^2 + abs(squeeze(fftut(mo{i+4}(2),j,1:nind,mo{i+4}(1)))).^2;
    %  umax(j) = 2*(max(sqrt(umaxTrms)))+0/(nz*nt);
end
yto{count} = n;
uto{count} = umax;
xto{count} = x;
Evor =x(1:end)*0;

for jj=1:length(Evor)
    j = jj;
    %v_2 = 2*(squeeze(abs(fftun(mo{i}(2),j,1:nind,mo{i}(1))).^2+abs(fftun(mo{i+4}(2),j,1:nind,mo{i+4}(1))).^2));
    %w_2 = 2*(squeeze(abs(fftuz(mo{i}(2),j,1:nind,mo{i}(1))).^2+abs(fftuz(mo{i+4}(2),j,1:nind,mo{i+4}(1))).^2));
%     v_v = 2*(squeeze(abs(fftun(mo{i}(2),j,:,mo{i}(1))).^2+abs(fftun(mo{i+4}(2),j,:,mo{i+4}(1))).^2));
%     v_v = sqrt(v_v);
%     w_v = 2*(squeeze(abs(fftuz(mo{i}(2),j,:,mo{i}(1))).^2+abs(fftuz(mo{i+4}(2),j,:,mo{i+4}(1))).^2));
%     w_v = sqrt(w_v);
%     vf = squeeze(Vtot(:,:,15));
%     ampp = sqrt((sum(vf.^2,2)));
%     prod = dot(vf,[v_v,w_v],2)./max(ampp);
%     %Edif = 0.5*(v_v.^2 + w_v.^2);
%     Edif = 0.5*prod.^2;
%     Edif(isnan(Edif)) =0;
%     Evor(jj) = trapz(n,Edif)./trapz(n,0.5*ampp.^2);
end

save('projected_vel.mat','qto','uto','xto','yto','qdns');

figure(6)
%semilogy(x(1:end),Evor)
plot(xfvec,Edifv)
hold on

figure(1)
subplot(1,4,count)
plot(x,umax-0*umax(xi),'Color',colr(count,:),'LineWidth',1.5)
box on
grid on
xlim([0,0.33])
%plot(x,umax-1*umax(xi),'r')
xlabel('$x$','Interpreter','latex','FontSize',14)
ylabel('$u_{max}$','Interpreter','latex','FontSize',14)


figure(3)
%plot(xfvec,errorv);
hold on
%plot(xfvec,errorw);
plot(xfvec,errorT,'o-','Color',colr(count,:),'LineWidth',1.);
xlabel('$x_f$','Interpreter','latex','FontSize',14)
ylabel('$MSE(\mathbf{u}_{dns},\mathbf{u}_{proj})$','Interpreter','latex','FontSize',14)
box on
grid on
%plot(xfvec,errorE,'o-');


figure(2)
hold on
plot(xfvec,errorE,'*-','Color',colr(count,:))
xlabel('$x_f$','Interpreter','latex','FontSize',14)
ylabel('$E_{proj}/E_{dns}$','Interpreter','latex','FontSize',14)
box on
grid on
%plot(xfvec,Edifv,'<-','Color',colr(count,:))


figure(4)
subplot(1,2,1)
hold on
dx = xfvec(2)-xfvec(1);
sumE = (cumsum(errorE)*dx);
%sumE= [errorE(1),sumE];
l1 = plot(x,umax-0*umax(xi),'Color',colr(count,:),'LineWidth',1.5,'DisplayName',num2str(i));
lableg(count) = l1;
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
xlim([0,0.33])
box on
grid on
legend(lableg)
%yline(1,'--')

end







