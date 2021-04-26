clear all
L = load('dominant_modes_05_sL.mat');
ft = L.ft*2*pi;
fz = L.fz*2*pi;
modes = L.mo;

xfvec = linspace(0.03,0.33,30);
N =150;
q0=zeros(4*N,1)+1;

for i=1:length(modes)
    mi = modes{i};
    beta = fz(mi(1));
    freq = -ft(mi(2))/(2*pi);
    for j=1:length(xfvec)
        xf = xfvec(j);
        [xw1,q1,y_opt1,En1,bbb,xar] = streak_readBL(xf,0.005,beta,freq,0,N,q0);

        q0 = q1(:,1);
        [xwj,qj,y_optj,Enj,bbb2,xarc2] = streak_readBL_no_opt(0.34,0.005,beta,freq,0,N,q0,y_opt1);
        xw{i,j} = xwj;
        q{i,j} = qj;
        En{i,j} = Enj;
        y_opt{i,j} = y_optj;

    end
end
save(['opt_modessL_difx.mat'],'xw','q','En','y_opt');


