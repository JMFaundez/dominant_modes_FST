clear all
casen = 1;

switch casen
    case 1
        L = load('dominant_modes_05_sL_hann.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.01+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
        ymaxp = 3e-3;
        maxvp=1e-3;
    case 2
        L = load('dominant_modes_05_lL_han.mat');
        ydm = 0.045;
        ymid = 5e-3;
        xiv =0.02+0*[0.05,0.04,0.02,0.01,0.05,0.025,0.02,0.015]; % large L
        ymaxp = 6.5e-3;
        maxvp=1.2e-3;
    case 3
        L = load('dominant_modes_3_sL_hann.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.0+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
        ymaxp = 15e-3;
        maxvp=2.2e-3;
    case 4
        L = load('dominant_modes_3_lL_han.mat');
        ydm = 0.045;
        ymid = 5e-3;
        xiv =0.02+0*[0.05,0.04,0.02,0.01,0.05,0.025,0.02,0.015]; % large L
        ymaxp = 18e-3;
        maxvp=2.2e-3;
    case 5
        L = load('dominant_modes_05_sL_N2.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.01+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
        ymaxp = 3e-3;
        maxvp=1e-3;
    case 6
        L = load('dominant_modes_05_sL_N34.mat');
        ydm = 0.015;
        ymid = 3e-3; % value used to generate dns mesh 3e-3 small L, 5e-3 large L
        xiv =0.01+0*[0.01,0.015,0.024,0.02,0.01,0.01,0.018,0.025]; % small L
        ymaxp = 3e-3;
        maxvp=1e-3;
end
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

colr=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];

% Bending function


ystar = 1-(n-ydm)/(max(n)-ydm);
S = 1./(1+exp(1./(ystar-1) +1./ystar));
S(ystar<=0) = 0;
S(ystar>=1)=1;


%%
count =0;

figure(1)
for i =[1:8]
    count= count +1;
    xind = find(x>=xiv(i),1,'first');
xi = x(xind);
N = length(n);
q0=zeros(4*N,1);
u0 = squeeze(uh{i}(xind,:)).*S*1;
v0 = squeeze(vh{i}(xind,:)).*S*1;
w0 = squeeze(wh{i}(xind,:)).*S*1;

figure(1)
subplot(2,4,count)
hold on
plot(abs(u0),n,'r')
plot(abs(v0),n,'b')
plot(abs(w0),n,'k')
legend('|u|','|v|','|w|')
xlabel('$|u|,|v|,|w|$','Interpreter','latex','FontSize',16)
ylabel('$n$','Interpreter','latex','FontSize',16)
box on
grid on
ylim([0,0.01])
xlim([0,maxvp])
title(num2str(i))


end