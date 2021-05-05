clear all
%close all
addpath('./matlab_script')
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m,n,p,[0.05 0.04],[0.12 0.05], [0.05 0.02]);
if ~make_it_tight, clear subplot;end

L = load('dominant_modes_05_lL.mat');
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

ny = length(n);
nx = length(x);

figure()
for i=1:8
mi = modes{i};
omega = ft(mi(2));
beta = fz(mi(1));
xff = 0.005;
xi = find(x>=xff,1,'first');

ui = squeeze(uh{i}(xi,:));
vi = squeeze(vh{i}(xi,:));
wi = squeeze(wh{i}(xi,:));

tf = 0.; 
%z = linspace(-0.015,0.015,100);
z = linspace(0,4*pi/abs(beta),100);
u = transpose(ui)*exp(1j*beta*z + 1j*omega*tf);
v = transpose(vi)*exp(1j*beta*z + 1j*omega*tf);
w = transpose(wi)*exp(1j*beta*z + 1j*omega*tf);

yf = 0.006;
Nr = 60; 
zstart = linspace(min(z),max(z),Nr);max(z)*rand(Nr,1); 
nstart = yf*ones(size(zstart));max(n)*rand(Nr,1); 
[Z,N] = meshgrid(z,n); 
subplot(2,4,i)
hold on
pcolor(Z,N,real(u))
[Zs,Ns,Xs] = meshgrid(z,n,1:2);
ws = Xs*0;
vs = Xs*0;
ws(:,:,1) = real(w);
vs(:,:,1) = real(v);
us = ws*10;
st = streamline(stream2(Z,N,real(w),real(v),zstart,nstart));
%st = streamslice(Zs,Ns,Xs,ws,vs,us,[],[],[1]);
set(st,'color','k')
set(st,'LineWidth',1.2)
%set(st,'AutoScale','on','AutoScaleFactor',0.5)
view(2)
shading interp
%quiver(Z,N,real(w),real(v)*0,'k')
ylim([0,yf])
colorbar()
end


figure()
plot(xdth,1.5*dth)




