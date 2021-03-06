clear all
casen = 2;
switch casen
    case 1
        L1 = load('dominant_modes_05_sL_hann.mat');
        L2 = load('dominant_modes_05_sL_hann_N34.mat');
        L3 = load('dominant_modes_05_sL_hann_N2.mat');
        ymaxp = 3e-3;
    case 2
        L1 = load('dominant_modes_05_lL_han.mat');
        L2 = load('dominant_modes_05_lL_han_N34.mat');
        L3 = load('dominant_modes_05_lL_han_N2.mat');
        ymaxp = 6.5e-3;
end

x = L1.X(:,1,1);
n = L1.N(1,:,1);
dth = L1.dth; 
xdth = L1.xa;
[x0dth, inx0] = min(abs(xdth));
xdth = xdth(inx0:end);
dth = real(dth(inx0:end));
uh1 = L1.uh;
uh2 = L2.uh;
uh3 = L3.uh;



for i=1:8
    umax1 = zeros(size(x)-1);
    umax2 = zeros(size(x)-1);
    umax3 = zeros(size(x)-1);
    for j=1:length(x)-1
        dthj = interp1(xdth,dth,x(j+1));
        nind = find(n>=3*dthj,1,'first');
        umax1(j) = max(abs(uh1{i}(j+1,1:nind)));
        umax2(j) = max(abs(uh2{i}(j+1,1:nind)));
        umax3(j) = max(abs(uh3{i}(j+1,1:nind)));
    end
    
figure(101)
subplot(2,4,i)
hold on
plot(x(2:end),umax1,'LineWidth',1.5,'DisplayName','N')
plot(x(2:end),umax2,'LineWidth',1.5,'DisplayName','3N/4')
plot(x(2:end),umax3,'LineWidth',1.5,'DisplayName','N/2')



xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
if i==4
legend()
end
box on
grid on
title(num2str(i))
xlim([0,0.33])
ylim([0,ymaxp])

figure(201)
err1 = abs((umax1-umax2)./(umax1));
err2 = abs((umax1-umax3)./(umax1));
subplot(2,4,i)
hold on
plot(x(2:end),err1,'LineWidth',1.5)
plot(x(2:end),err2,'LineWidth',1.5)

xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u_{max}$','Interpreter','latex','FontSize',16)
box on
grid on
title(num2str(i))
xlim([0.01,0.33])
    
end
