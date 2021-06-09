clear all
close all

addpath('./matlab_script')
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m,n,p,[0.12 0.08],[0.15 0.08], [0.08 0.04]);
if ~make_it_tight, clear subplot;end
casen = 1;

switch casen
    case 1
        L = load('dominant_modes_05_sL_hann.mat');
        ininame = 'opt_mode';
      fnamedns = 'modes_dns_smallL_05p.mat';
        fnameopt = 'modes_opt_smallL_05p.mat';
    case 2
        L = load('dominant_modes_05_lL_han_N34.mat');
        ininame = 'optw_mode';
        ylimp200 = [0,0.2];
        ylimp300 = [0,6e-10];
    case 3
        L = load('dominant_modes_3_sL_hann.mat');
        ininame = 'opt_mode';
        fnamedns = 'modes_dns_smallL_3p.mat';
        fnameopt = 'modes_opt_smallL_3p.mat';
        
    case 4
        L = load('dominant_modes_3_lL_han.mat');
        ininame = 'optw_mode';
        fnamedns = 'modes_dns_largeL_3p.mat';
        fnameopt = 'modes_opt_largeL_3p.mat';
end
uh = L.uh;
vh = L.vh;
wh = L.wh;
x = L.X(:,1,1);
n = L.N(1,:,1);
Ndns = length(n);
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

colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];


%ki = [length(xfvec),3,5,2,11,3,length(xfvec),length(xfvec)];
omega ={};
beta = {};
ui = {};
vi = {};
wi = {};
uic = {};
vic = {};
wic = {};
xidns = {};
ydns = {};
uopt = {};
vopt = {};
wopt = {};
yopt = {};
xiopt = {};

uosum = {};
vosum = {};
wosum = {};
Umatrix = {}

for i=1:8
    O = load([ininame,num2str(i),'_difx.mat']);
    
    xop = O.xw{1};
    xff = xop(1);
    xi = find(x>=xff,1,'first');
    xfvec = linspace(xop(1)*1.5,0.33,30);
    x(xi);
    % DNS
    ui{i} = squeeze(uh{i}(xi,:));
    vi{i} = squeeze(vh{i}(xi,:));
    wi{i} = squeeze(wh{i}(xi,:));
    uic{i} = ui{i};
    vic{i} = vi{i};
    wic{i} = wi{i};
    nk = length(xfvec);
    Umatrix{i} = zeros(2*length(n),nk);
    
    uosum{i} = 0;
    vosum{i} = 0;
    wosum{i} = 0;
    
    omega{i} = ft(modes{i}(2));
    beta{i} = fz(modes{i}(1));
    xidns{i} = x(xi);
    ydns{i} = n;
    
    uenvo = 0*O.xw{1};
    Edifv = zeros(size(xfvec));
    
    
    for k =1:length(xfvec)%[10,20]
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
        
        
        aint =  conj(vi{i}).*vo + conj(wi{i}).*wo;
        a = 0.5*trapz(n,aint);
        scale = a;
        %ui{i} = ui{i}-uo*scale;
        %vi{i} = vi{i}-vo*scale;
        %wi{i} = wi{i}-wo*scale;
        
        uopt{i,k} = uop*scale;
        vopt{i,k} = vop*scale;
        wopt{i,k} = wop*scale;
        yopt{i} = yop;
        xiopt{i} = xff;
        Umatrix{i}(1:Ndns,k) = vo*scale;
        Umatrix{i}(Ndns+1:end,k) = wo*scale;
        uosum{i} = uosum{i} + uop*scale;
        vosum{i} = vosum{i} + vop*scale;
        wosum{i} = wosum{i} + wop*scale;
   
    end
end
xfopt = xfvec;

save(fnamedns,'ui','vi','wi','omega','beta','ydns','xidns');
save(fnameopt,'uopt','vopt','wopt','omega','beta','yopt','xiopt','xfopt','uosum','vosum','wosum');
%%
for i=1:8
    ran = rank(Umatrix{i}(1:Ndns,:))
end

for i=1:8
    nk = length(xfvec);
    uosum{i} = uopt{i,nk};
    vosum{i} = vopt{i,nk};
    wosum{i} = wopt{i,nk};
    for k= nk:-1:1
        aint =  conj(vopt{i,k}).*vosum{i} + conj(wopt{i,k}).*wosum{i};
        a = 0.5*trapz(n,aint);
        a = a;
        difu =  uopt{i,k} - uosum{i}*a;
        difv =  vopt{i,k} - vosum{i}*a;
        difw =  wopt{i,k} - wosum{i}*a;
        uosum{i} = uosum{i} + difu;
        vosum{i} = vosum{i} + difv;
        wosum{i} = wosum{i} + difw;
    end
end

for i=1:8
    subplot(2,4,i)
    hold on
    plot(abs(vosum{i}),yopt{i},'r--')
    plot(abs(wosum{i}),yopt{i},'b--')
    plot(abs(vic{i}),ydns{i},'r')
    plot(abs(wic{i}),ydns{i},'b')
    ylim([0,8e-3])
    Esum = -trapz(yopt{i},abs(vosum{i}).^2 +  abs(wosum{i}).^2)
    Edns = trapz(ydns{i},abs(vic{i}).^2 +  abs(wic{i}).^2)
    
end


