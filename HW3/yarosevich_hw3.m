% Simulate the repressilator model
clc; clear all;
close all

alpha=50;

alpha0=0;
beta=.5;
n=2;
gamma = -.27;


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

legend({'$m_1$', '$m_2$', '$m_3$', '$p_1$', '$p_2$','$p_3$','$n_1$', '$n_2$', '$n_3$', '$q_1$', '$q_2$','$q_3$'}, 'Interpreter', 'latex')
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

%%
clc; clear all; close all;

% Simulate the repressilator model
% Positive Feedforward loop

alpha=2;
beta =2;
gamma = .5;
n = 2;
k=1;

p = [alpha, beta, gamma, n, k];
y0 = [0, 0];
tmax = 20;

[T,Y] = ode45(@posfeed,[0, tmax],y0,[],p);



subplot(2,1,1);
plot(T,Y(:,1), 'k', 'LineWidth', 2)
title('X Behavior')
xlabel('t')
subplot(2,1,2);
plot(T,Y(:,2), 'b',  'LineWidth', 2)
xlabel('t')
title('Y Behavior')

%Seems to work perfectly! Had to tune the strength of the z effect and the
%rate of decay.

%% Double negative feedback
clc; clear all; close all;

alpha=2;
beta =2;
gamma = .5;
n = 2;
k=1;

p = [alpha, beta, gamma, n, k];
y0 = [0, 1];
tmax = 10;

[T,Y] = ode45(@negfeed,[0, tmax],y0,[],p);



subplot(2,1,1);
plot(T,Y(:,1), 'k', 'LineWidth', 2)
title('X Behavior')
subplot(2,1,2);
plot(T,Y(:,2), 'b',  'LineWidth', 2)
title('Y Behavior')

% This works, but it lacks some of the intuition from the prompt, at least
% insofar as I'm basing Y's decay on x in order to represent repression by
% X).
%%
clc; clear all; close all;

alpha=.001;
beta = 5;
gamma =1;
n = 1;
k=1;

p = [alpha, beta, gamma, n, k];
y0 = [0, 0, 0, 0, 0, 0, 0];
tmax = 20;

[T,Y] = ode45(@sharktooth,[0, tmax],y0,[],p);

figure(1)
plot(T,Y(:,2), 'k', 'LineWidth', 2)
title('')
hold on;
plot(T, Y(:,5), 'b', 'LineWidth', 2)
plot(T, Y(:,7), 'm', 'LineWidth', 2)
legend({'$Z_1$', '$Z_2$', '$Z_3$'}, 'Interpreter', 'latex')







