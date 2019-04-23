EKF;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotting
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotID = 3001;
figure(plotID);
set(plotID, 'Position', [1 700 600 300], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'PaperPositionMode', 'auto');
hold on;
plot(IEKFitcount, 'b');
title('IEKF iterations at each sample');
if (printfigs == 1)
    fpath = sprintf('fig_demoKFStatesEstimatesMeasurements');
    savefname = strcat(figpath, fpath);
    print(plotID, '-dpng', '-r300', savefname);
end

plotID = 3002;
figure(plotID);
    subplot(3,1,1)    
        plot([0:N]*dt,Z_k(1,:))
        hold on
        plot([0:N]*dt,Z_k1k1(1,:))
        xlim([0 N*dt])
        xlabel('time[s]')
        ylabel('\alpha (AoA) [rad]');
        title('Output');
        hold off
    subplot(3,1,2)
        plot([0:N]*dt,Z_k(2,:))
        hold on
        plot([0:N]*dt,Z_k1k1(2,:))
        xlim([0 N*dt])
        xlabel('time[s]')
        ylabel('\beta (Sideslip) [rad]');
        hold off
    subplot(3,1,3)
        plot([0:N]*dt,Z_k(3,:))
        hold on
        plot([0:N]*dt,Z_k1k1(3,:))
        xlabel('time[s]')
        xlim([0 N*dt])
        ylabel('Velocity V [m/s]');
        hold off


plotID = 3003;
figure(plotID);
title('Error');
subplot(3,1,1);
    plot([0:N]*dt,errrr(1,:))
    xlim([0 N*dt])
    xlabel('time[s]');
    ylabel('\alpha (AoA) [rad]');
    title('Error');
subplot(3,1,2);
    plot([0:N]*dt,errrr(2,:))
    xlim([0 N*dt])
    xlabel('time[s]');
    ylabel('\beta (Sideslip) [rad]');
subplot(3,1,3);
    plot([0:N]*dt,errrr(3,:))
    xlim([0 N*dt])
    xlabel('time[s]');
    ylabel('Velocity V [m/s]');

plotID = 3004;
figure(plotID);
title('State Estimation');
subplot(4,1,1);
    plot([0:N-1]*dt,XX_k1k1(1,:))
    xlim([0 N*dt])
    xlabel('time[s]');
    ylabel('U [m/s]');
    title('State Estimates');
subplot(4,1,2);
    plot([0:N-1]*dt,XX_k1k1(2,:))
    xlim([0 N*dt])
    xlabel('time[s]');
    ylabel('V [m/s]');
subplot(4,1,3);
    plot([0:N-1]*dt,XX_k1k1(3,:))
    xlim([0 N*dt])
    xlabel('time[s]');
    ylabel('W [m/s]');
subplot(4,1,4);
    plot([0:N-1]*dt,XX_k1k1(4,:))
    xlim([0 N*dt])
    xlabel('time[s]');
    ylabel('C_{a}');
Ints=100;
alpha = Z_k1k1(1, Ints:N);
beta = Z_k1k1(2, Ints:N);
xinp=[alpha;beta];
plotID = 3005;
figure(plotID);
subplot(1,2,1);
    plot(alpha,beta)
    title('Filtered')
    xlabel('\alpha (AoA)');
    ylabel('\beta (sideslip)');
subplot(1,2,2);
    scatter(Z_k(1,:),Z_k(2,:),'.');
    title('Raw data');
    xlabel('\alpha (AoA)');
    ylabel('\beta (sideslip)');
    %%
    x=alpha';
    y=beta';
    z=Cm(Ints:N);
    plotID = 3006;
figure(plotID);
    k=1;
    for i=1:3
        for j=1:3
            subplot(3,3,k)
            k=k+1;
            s=sprintf('poly%x%x',i,j)
            sf = fit([x, y],z,s);
            plot(sf,[x,y],z)
            title("Polynomial Fit "+string(s));
            xlabel('\alpha (AoA)');
            ylabel('\beta (Sideslip)');
            zlabel('C_{m}');
        end
    end