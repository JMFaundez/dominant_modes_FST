clear all

D = load('modes_dns_largeL_3p.mat');
O = load('modes_opt_largeL_3p.mat');

%D = load('modes_dns_smallL_3p.mat');
%O = load('modes_opt_smallL_3p.mat');

ui = D.ui;
vi = D.vi;
wi = D.wi;
om = D.omega;
beta = D.beta;
ydns = D.ydns;
xi = D.xidns;


uo = O.uopt;
vo = O.vopt;
wo = O.wopt;
yopt = O.yopt;
xiopt = O.xiopt;
xfopt = O.xfopt;


figure('Position',[500,500,1000,500])
for i=1:8
    subplot(2,4,i)
    hold on
    plot(abs(ui{i}),ydns{i},'r')
    plot(abs(vi{i}),ydns{i},'b')
    plot(abs(wi{i}),ydns{i},'k')
    ylim([0,5e-3])
    title(['\omega=',num2str(om{i},'%1.1f'),',\beta=',num2str(beta{i},'%1.1f'),...
        'xi=',num2str(xi{i},'%1.2f')])
end

k = 1;
figure('Position',[500,500,1000,500])
for i=1:8
    subplot(2,4,i)
    hold on
    plot(abs(uo{i,k}),yopt{i},'r')
    plot(abs(vo{i,k}),yopt{i},'b')
    plot(abs(wo{i,k}),yopt{i},'k')
    ylim([0,5e-3])
    title(['\omega=',num2str(om{i},'%1.1f'),',\beta=',num2str(beta{i},'%1.1f'),...
        'xi=',num2str(xiopt{i},'%1.2f')])
end

sgtitle(['x_f=',num2str(xfopt(k))]) 