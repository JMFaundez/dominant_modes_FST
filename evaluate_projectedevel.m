clear all
L = load('projected_vel.mat');
umax= L.uto;
xmax = L.xto;
qp = L.qdns;
[nm,Nxf,~] = size(qp);
betav = [1,1,2,3]*2*pi/0.02;
freqv = -[10 15 19 25]/(2*pi);
i = 4;
colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];
beta = betav(i);
freq = freqv(i);
N = length(qp{1,1,1});
q0=zeros(4*N,1)+1;
y = L.yto{i};
ydm = 0.018;
ystar = 1-(y-ydm)/(max(y)-ydm);
S = 1./(1+exp(1./(ystar-1) +1./ystar));
S(ystar<=0) = 0;
S(ystar>=1)=1;

jj=0;
for j=1:1
    jj=jj+1;
y = L.yto{i};
u0 = qp{i,j,1}.*S'*1;
v0 = qp{i,j,2}.*S';
w0 = qp{i,j,3}.*S';
q0(1:N,1) = flip(u0);
q0(N+1:2*N,1) = flip(v0);
q0(2*N+1:3*N,1) = flip(w0);
[xwj,qj,y_optj,Enj,bbb2,xarc2] = streak_readBL_no_opt(0.34,0.005,beta,freq,0,N,q0,flip(y));
xw2{jj} = xwj;
q{jj} = qj;
En2{jj} = Enj;
y_opt2{jj} = y_optj;

end
%%


figure()
hold on
jj=0;
for j=1:1
    jj=jj+1;
    [nn,Nstations] = size(q{jj});
    maxop = zeros(Nstations,1);
    for k=1:Nstations
        maxop(k) = max(abs(q{jj}(1:N,k)));
    end
    xop = xw2{jj};
    if j==11
        plot(xop,maxop,'Color','r','Linewidth',1.5)
    else
    plot(xop,maxop,'Color',[0.5,0.5,0.5])
    end
end
plot(xmax{i},umax{i},'Color',colr(i,:),'LineWidth',1.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
box on
grid on
xlim([0,0.33])




k=1;
yop = y_opt2{k};
figure()
hold on
plot(abs(q{k}(1:N)).*flip(S),yop)
plot(abs(q{k}(N+1:2*N)),yop)
plot(abs(q{k}(2*N+1:3*N)),yop)
%ylim([0,0.01])





