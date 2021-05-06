clear all

L = load('dominant_modes_05_sL.mat');
uh = L.uh;
vh = L.vh;
wh = L.wh;
x = L.X(:,1,1);
n = L.N(1,:,1);
ft = L.ft*2*pi;
fz = L.fz*2*pi;
nz = length(fz);
modes = L.mo;
nmi =length(modes);
modes{nmi+1}(1) = 3;
modes{nmi+1}(2) = 11;
modes{nmi+2}(1) = nz-1;
modes{nmi+2}(2) = 11;
modes{nmi+3}(1) = 2;
modes{nmi+3}(2) = 33;
modes{nmi+4}(1) = nz;
modes{nmi+4}(2) = 33;
% modes{nmi+5}(1) = 3;
% modes{nmi+5}(2) = 38;
% modes{nmi+6}(1) = nz-1;
% modes{nmi+6}(2) = 38;

nm = length(modes);
omega = zeros(nm*2,1);
beta = zeros(nm*2,1);

for i=1:nm
    omega(i) = ft(modes{i}(2));
    beta(i) = fz(modes{i}(1));
end

beta(nm+1:end) = beta(1:nm);
omega(nm+1:end) = -omega(1:nm);
%openfig('spectrum_3_sL_x002.fig');
openfig('spectrum_3_sL_x015.fig');
hold on
scatter(omega,beta,'ko','LineWidth',1.5);
for i =1:2*nm-1
    for j=i+1:2*nm
        db = beta(i)+beta(j);
        do = omega(i)+omega(j);
        plot(do,db,'rs')
    end
end
xlim([0,65])
openfig('spectrum_3_sL_x002.fig');
hold on
scatter(omega,beta,'ko','LineWidth',1.5);
xlim([0,65])
