clear all
L = load('dominant_modes_05_sL.mat');
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
ymid = 3e-3; % value used to generate dns mesh

colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];

% Bending function
ydm = 0.0145;
ystar = 1-(n-ydm)/(max(n)-ydm);
S = 1./(1+exp(1./(ystar-1) +1./ystar));
S(ystar<=0) = 0;
S(ystar>=1)=1;

xi = 0.005;
xind = find(x>=0.005,1,'first');


%%
count =0;
for i =4:4
    count = count +1;
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
maxop = max(abs(q(1:N,:)),[],1);
maxop2 = max(abs(q2(1:N,:)),[],1);

%%
figure(100)
%subplot(1,4,count)
hold on
plot(x(2:end),umax,'Color',colr(i,:),'LineWidth',1.5)
plot(xw,maxop,'Color','k');%[0.5,0.5,0.5])
plot(xw,maxop2,'--','Color','k')%[0.5,0.5,0.5])
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
box on
grid on
xlim([0,0.33])

end
%%

yop = y_opt;
figure()
hold on
plot(abs(q(1:N,1)),yop)
plot(abs(q(N+1:2*N,1)),yop)
plot(abs(q(2*N+1:3*N,1)),yop)
legend('|u|','|v|','|w|')
box on
grid on
ylim([0,0.01])

[uhmax,inn] = max(abs(u0));
nmax = n(inn);
nratio =nmax./real(dth);
nratio(isnan(nratio))=0;
figure()
plot(xdth,nratio)
xlim([0.005,0.35])







