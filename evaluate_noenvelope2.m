clear all


beta = -1*2*pi/0.02;
freq = 15/(2*pi);
xfvec = linspace(0.03,0.33,30);
N =150;
q0=zeros(4*N,1)+1;

for j=1:length(xfvec)
xf = xfvec(j);
[xw,q1,y_opt,En,bbb,xar] = streak_readBL(xf,0.005,beta,freq,0,N,q0);


q0 = q1(:,1);
[xwj,qj,y_optj,Enj,bbb2,xarc2] = streak_readBL_no_opt(0.34,0.005,beta,freq,0,N,q0,y_opt);
xw2{j} = xwj;
q{j} = qj;
En2{j} = Enj;
y_opt2{j} = y_optj;

end
save(['opt_m6_difx.mat'],'xw2','q','En2','y_opt2');


