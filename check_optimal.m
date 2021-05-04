clear all
%close all
figure()

j = 0;
for i=1:4
    j=j+1;
    O = load(['opt_mode',num2str(i),'_difx.mat']);
for k=1:15

% Opt
[nn2,Nstations2] = size(O.q{k});
N = nn2/4;
yop = O.y_opt{k};
q0 = O.q{k}(:,1);
uop = (q0(1:N,1));
wop = (q0(2*N+1:N*3,1));
vop = (q0(N+1:N*2,1));
subplot(2,4,j)
ylim([0,0.08])
hold on
plot(abs(vop),yop,'k')
subplot(2,4,j+4)

hold on
plot(abs(wop),yop,'k')
ylim([0,0.08])
end
end