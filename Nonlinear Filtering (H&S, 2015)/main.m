%==========================================================================
%                       DSGE MODEL ESTIMATION:  
%              Particle Filter Approximation of Likelihood 
%
%
% Author: Minsu Chang        minsuc@sas.upenn.edu
% Last modified: 2/24/2016
%==========================================================================


clear
clc
close all
delete *.asv

tic

l = path;

path('Mfiles',path);
path('LRE',path);

% load data and consider parameters in Table 8.1.

yt      = load('us.txt');
param_m = [2.09 0.98 2.25 0.65 0.34 3.16 0.51 0.81 0.98 0.93 0.19 0.65 0.24];
param_l = [3.26 0.89 1.88 0.53 0.19 3.29 0.73 0.76 0.98 0.89 0.20 0.58 0.29];

% check whether likelihood in Table 8.1. is replicated
% dsgeliki_text(param_m)

[T1, ~, T0, ~, ~, ~] = model_solution(param_m);
[A,B,H,R,S2,Phi] = sysmat(T1,T0,param_m);

% Kalman filter result
[liki, measurepredi, statepredi, varstatepredi] = kalman(A,B,H,R,S2,Phi,yt);

ns      = size(B,2);
T       = size(yt,1);

% initialize
x0 = zeros(ns,1);
P0 = nearestSPD(dlyap(Phi, R*S2*R'));  % to make it positive semidefinite

N       = 400; % number of particles (for BSPF: 40000, for COPF: 400)

[lik, all_s_up, Neff] = PF_lik(A, B, H, Phi, R, S2, N, yt, x0, P0, 0);
% Last input denotes the indicator for bootstrap particle filtering.
% if == 1, it performs bootstrap particle filtering. otherwise, it does
% conditionally-optimal particle filtering.

sum(lik)


%=========================================================================
%                  FIGURE 1: Log Likelihood Increments
%=========================================================================


figure('Position',[20,20,900,600],'Name',...
    'Log Likelihood Increments','Color','w')

plot(liki,'LineStyle','-','Color','b','LineWidth',2.5) % from Kalman filter
hold on 
plot(lik,'LineStyle','--','Color','r','LineWidth',2.5) % from CO Particle filter
title('$ln \hat{p} (y_{t}|Y_{1:t-1}, \theta^{m}$) vs. $lnp(y_{t}|Y_{1:t-1}, \theta^{m})$','Interpreter','latex')


%=========================================================================
%                  FIGURE 2: Filtered States 
%=========================================================================

figure('Position',[20,20,1200,800],'Name',...
    'Filtered States','Color','w')

subplot(3,1,1)
plot(statepredi(2:end,5),'LineStyle','-','Color','b','LineWidth',2.5) % g from Kalman filter
hold on 
plot(mean(all_s_up(:,5,:),3),'LineStyle','--','Color','r','LineWidth',2.5) % from CO Particle filter
title('$\hat{E} (\hat{g}_{t} |Y_{1:t}, \theta^{m}$) vs. $E(\hat{g}_{t}|Y_{1:t}, \theta^{m})$','Interpreter','latex')
subplot(3,1,2)
plot(statepredi(2:end,6),'LineStyle','-','Color','b','LineWidth',2.5) % z from Kalman filter
hold on 
plot(mean(all_s_up(:,6,:),3),'LineStyle','--','Color','r','LineWidth',2.5) % from CO Particle filter
title('$\hat{E} (\hat{z}_{t} |Y_{1:t}, \theta^{m}$) vs. $E(\hat{z}_{t}|Y_{1:t}, \theta^{m})$','Interpreter','latex')
subplot(3,1,3)
plot(statepredi(2:end,1),'LineStyle','-','Color','b','LineWidth',2.5) % y from Kalman filter
hold on 
plot(mean(all_s_up(:,1,:),3),'LineStyle','--','Color','r','LineWidth',2.5) % from CO Particle filter
title('$\hat{E} (\hat{y}_{t} |Y_{1:t}, \theta^{m}$) vs. $E(\hat{y}_{t}|Y_{1:t}, \theta^{m})$','Interpreter','latex')


%=========================================================================
%                  FIGURE 3: Effective Sample Size
%=========================================================================


figure('Position',[20,20,900,600],'Name',...
    'Effective Sample Size','Color','w')

plot(Neff,'LineStyle','-','Color','b','LineWidth',2.5) % from CO Particle filter
title('Effective Sample Size')






path(l);
disp(['         ELAPSED TIME:   ', num2str(toc)]);
elapsedtime=toc;



