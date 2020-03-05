%==========================================================================
%                       DSGE MODEL ESTIMATION:  
%              Particle Filter Approximation of Likelihood 
%
%
% Juan Castellanos Silván
% Date: 03/04/2020
%==========================================================================

clear all, clc

%=========================================================================
%                            EXPERIMENTS
%=========================================================================

sample = 1;         % 1 => Old sample 1983:I - 2002:IV
                    % 2 => New sample 1999:IV - 2019:III

accuracy = 1;       % 0 => filters only run once
                    % 1 => filters run multiple times

if sample == 1
    xAxis = 1983:0.25:2002.75;
    a = 1983;
    b = 2002;
else
    xAxis = 1999.75:0.25:2019.5;
    a = 1999;
    b = 2019;
end
    
%=========================================================================
%                           LOAD RESULTS 
%=========================================================================

workpath = '/Users/Castesil/Documents/GitHub/Econ 722 - Schorfheide/PS4/particle_filters';
savepath = [workpath, '/results/'];

for naive = [0,1]
    
    loadfilename = ['PF_naive', num2str(naive),'_sample', num2str(sample),'_accuracy', num2str(accuracy),'.mat'];

    cd(savepath);
    load(loadfilename);

    if naive == 1 
        lik1 = lik;
        all_s_up1 = all_s_up;
        Neff1 = Neff;
        clear lik all_s_up Neff
        
        if accuracy == 1
            Delta11 = Delta1; 
            clear Delta1
        end
        
    else
        lik0 = lik;
        all_s_up0 = all_s_up;
        Neff0 = Neff;
        clear lik all_s_up Neff
        
        if accuracy == 1
            Delta10 = Delta1;
            clear Delta1
        end
        
    end
       
end

cd(workpath);

[density0, x0_den] =ksdensity(mean(Delta10,2));
[density1, x1_den] = ksdensity(mean(Delta11,2));
return
%=========================================================================
%                  FIGURE 1: Log Likelihood Increments
%=========================================================================


figure('Position',[20,20,900,600],'Name',...
    'Log Likelihood Increments','Color','w')

plot(xAxis,liki,'LineStyle','-','Color','k','LineWidth',1.5) % from Kalman filter
hold on 
plot(xAxis,lik0, 'LineStyle',':','Color','k','LineWidth',2) % from CO Particle filter
hold on
plot(xAxis,lik1,'LineStyle','--','Color','k','LineWidth',1.5) % from BS Particle filter
grid on 
title('$ln\:\hat{p}(y_{t}|Y_{1:t-1}, \theta^{m}$) vs. $ln\:p(y_{t}|Y_{1:t-1}, \theta^{m})$','fontsize', 15,'interpreter','latex')
xlim([a,b])

x = 29.7;                  % A4 paper size
y = 21.0;                  % A4 paper size
xMargin = 1;               % left/right margins from page borders
yMargin = 1;               % bottom/top margins from page borders
xSize = x - 2*xMargin;     % figure size on paper (widht & hieght)
ySize = y - 2*yMargin;     % figure size on paper (widht & hieght)

set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[x y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

save = '/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS4/LaTeX/';
filename = strcat('pLikelihood', num2str(sample),'.pdf');
saveas(gcf, strcat(save,filename));


%=========================================================================
%                  FIGURE 2: Filtered States 
%=========================================================================

figure('Position',[20,20,1200,800],'Name',...
    'Filtered States','Color','w')

subplot(3,1,1)
plot(xAxis, statepredi(2:end,5),'LineStyle','-','Color','k','LineWidth',1.5) % g from Kalman filter
hold on 
plot(xAxis, mean(all_s_up0(:,5,:),3),'LineStyle',':','Color','k','LineWidth',2) % from CO Particle filter
hold on 
plot(xAxis, mean(all_s_up1(:,5,:),3),'LineStyle','--','Color','k','LineWidth',1.5) % from BS Particle filter
grid on 
title('$\hat{E} (\hat{g}_{t} |Y_{1:t}, \theta^{m}$) vs. $E(\hat{g}_{t}|Y_{1:t}, \theta^{m})$','fontsize', 15,'interpreter','latex')
xlim([a,b])

subplot(3,1,2)
plot(xAxis,statepredi(2:end,6),'LineStyle','-','Color','k','LineWidth',1.5) % z from Kalman filter
hold on 
plot(xAxis,mean(all_s_up0(:,6,:),3),'LineStyle',':','Color','k','LineWidth',2) % from CO Particle filter
hold on
plot(xAxis,mean(all_s_up1(:,6,:),3),'LineStyle','--','Color','k','LineWidth',1.5) % from BS Particle filter
grid on
title('$\hat{E} (\hat{z}_{t} |Y_{1:t}, \theta^{m}$) vs. $E(\hat{z}_{t}|Y_{1:t}, \theta^{m})$','fontsize', 15,'interpreter','latex')
xlim([a,b])

subplot(3,1,3)
plot(xAxis,statepredi(2:end,1),'LineStyle','-','Color','k','LineWidth',1.5) % y from Kalman filter
hold on 
plot(xAxis,mean(all_s_up0(:,1,:),3),'LineStyle',':','Color','k','LineWidth',2) % from CO Particle filter
hold on
plot(xAxis,mean(all_s_up1(:,1,:),3),'LineStyle','--','Color','k','LineWidth',1.5) % from BS Particle filter
grid on
title('$\hat{E} (\hat{y}_{t} |Y_{1:t}, \theta^{m}$) vs. $E(\hat{y}_{t}|Y_{1:t}, \theta^{m})$','fontsize', 15,'interpreter','latex')
xlim([a,b])

set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[x y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

save = '/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS4/LaTeX/';
filename = strcat('pFilteredStates', num2str(sample),'.pdf');
saveas(gcf, strcat(save,filename));
