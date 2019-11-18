% Simulate the repressilator model
clc; clear all;
close all

alpha=50;

alpha0=0;
beta=.2;
n=2;
gamma = -2;


%for numreps=1:4 

p = [alpha,alpha0,beta,n, gamma];
x0 = 30*rand(12,1) ;

Tmax=100;

[T,Y] = ode45(@repress,[0 Tmax],x0,[],p);

figure(1)
plot(T,Y(:,1:3),'LineWidth',3) ; hold on;
plot(T,Y(:,4:6),':','LineWidth',3) ; hold on;
plot(T,Y(:,7:9),':','LineWidth',3) ; hold on;
plot(T,Y(:,10:12),':','LineWidth',3) ; hold on;

legend('m lalcl','m tetR','m cl','p lacl','p tetR','p cl')
xlabel('t') ; 
set(gca,'FontSize',16)


% figure(2)
% plot3(Y(:,1),Y(:,2),Y(:,3),'LineWidth',3) ; hold on 
% xlabel('m lalcl');ylabel('m tetR');zlabel('m cl')
% set(gca,'FontSize',16)
%end

% So it would seem even very small gamma values cause the coupled protein
% values to converge to one another so that there are 6 distinct
% trajectories, corresponding to 6 sets of concentrations. The overall
% reason for this is because the gamma term creates a coupling that is
% symmetric, in that each concentration now has a damping effect on one
% other (in some cases this relationship is 'twice removed'). The result is
% that if any particular concentration gets larger, it will, circuitously,
% give a bump to another concetration that inhibits it. 