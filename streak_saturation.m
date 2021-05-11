clear all
freq = 0;
beta = 150;
xi = 0.008;
N=150;

q0=zeros(4*N,1)+1;
[xw1,q1,y_opt1,En1,bbb,xar] = streak_readBL(0.33,xi,beta,freq,0,N,q0);

q0 = q1(:,1);
q0(1:N) = q1(1:N,40);
[xw2,q2,y_opt2,Enj,bbb2,xarc2] = streak_readBL_no_opt(0.33,xi,beta,freq,0,N,q0,y_opt1);

%%
[NN,Nstations] = size(q2);
qmax1 = zeros(Nstations,1);
qmax2 = zeros(Nstations,1);
for j=1:Nstations
    qmax1(j) = max(abs(q1(1:N,j)));
    qmax2(j) = max(abs(q2(1:N,j)));
end


figure()
hold on
plot(xw1,qmax1)
plot(xw2,qmax2)
