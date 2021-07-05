clear all
L = load('dominant_modes_05_sL_han_N.mat');
ft = L.ft*2*pi;
fz = L.fz*2*pi;
modes = L.mo;
xiv = 0.01 + zeros(length(modes),1) ;[0.05,0.04,0.02,0.01,0.05,0.025,0.02,0.015,0.02,0.02]; % large L
%xiv =+0.01+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
N =150;
q0=zeros(4*N,1)+1;
%%
for i=1:length(modes)
    mi = modes{i};
    beta = fz(mi(1));
    freq = -ft(mi(2))/(2*pi);
    xfvec = linspace(xiv(i)*1.5,0.33,30);
    for j=1:length(xfvec)
        xf = xfvec(j);
        [xw1,q1,y_opt1,En1,bbb,xar] = streak_readBL(xf,xiv(i),beta,freq,0,N,q0);

        q0 = q1(:,1);
        [xwj,qj,y_optj,Enj,bbb2,xarc2] = streak_readBL_no_opt(0.34,xiv(i),beta,freq,0,N,q0,y_opt1);
        xw{j} = xwj;
        q{j} = qj;
        En{j} = Enj;
        y_opt{j} = y_optj;

    end
    save(['opt_mode',num2str(i),'_difx.mat'],'xw','q','En','y_opt');
end



