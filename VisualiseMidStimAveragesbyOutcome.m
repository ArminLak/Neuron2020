rewardColor = [0 0.6 0
    0 1 0
    0.8 0 0];
color = 3;

downSample = 1200;

start = -0.2 % s this should be -1 or less
stop = 0.8   % s

load('midStimAveragesbyOutcome.mat', 'a');

figure;
hold on 


    plot(nanmean(a(1:353,:)), 'color', rewardColor(3,:),'LineWidth',2)
    hold on
    plot(nanmean(a(354:1375,:)), 'color', rewardColor(2,:),'LineWidth',2)
    hold on 
    plot(nanmean(a(1376:2521,:)), 'color', rewardColor(1,:),'LineWidth',2)
    
    

legend('No reward','Small reward', 'Large reward','location','best')
     title('stimulus aligned for 0.25 contrast')
     xlim([3700-abs(start*downSample) 3700+abs(stop*downSample)])
%      xlim([3500 7400])
     xticks([3700 3700+abs(stop*downSample)])
     xticklabels({'0','0.8'})
     ylabel('Norm {\Delta} F / F', 'FontWeight', 'bold')

    
     set(gca, 'TickDir', 'out')
     xlabel('Time (s)', 'FontWeight', 'bold')
  
    