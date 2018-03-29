clear all
close all

% list of animals
Animals = [48 50 51]

load('BehPhotoM_Exp23')

 

color = [
    1 0 0         % red
    1 0.5 0       % orange
    1 1 0         % yellow
    0.5 1 0.5     % light green
    0 1 1         % light blue
    0  0.5 1      % medium blue
    0 0 1];       % blue

colorGray = [ 0.8 0.8 0.8
    0.6 0.6 0.6
    0.4 0.4 0.4
    0 0 0
    ];

StimAllowed = [-0.5 -0.25 -0.12 0 0.12 0.25 0.5];

GrandPopAbsStimResp = zeros(4,7100);
GrandPopAbsStimRespMiddleStimCorrect= zeros(1,7100);
GrandPopAbsStimRespMiddleStimError= zeros(1,7100);

GrandPopAbsStimRespMiddleStimCorLarge= zeros(1,7100);
GrandPopAbsStimRespMiddleStimCorSmall= zeros(1,7100);



GrandPopNormBinStimNoFold = zeros(2,7);
GrandPopNormBinRewardNoFold = zeros(2,7);
GrandPopNormBinStimNoFoldCorrError = zeros(2,7);
GrandPopNormBinRewardNoFoldCorrError = zeros(2,7);


GrandPopStim2AwayDACorError = zeros(4,7100);

GrandPopBeep2AwayDACorError = zeros(4,7100);

GrandPopRew2AwayDACorError = zeros(4,7100);

GrandPopStimBin = zeros(4,4);

GrandPopRewBin = zeros(4,4);

c=1;

for iAnimal = Animals
    
    PerBlock1(c,:)=BehPhotoM(iAnimal).GrandSummary.Performance(1,:);
    RTBlock1(c,:)=BehPhotoM(iAnimal).GrandSummary.RT(1,:);
     
    PerBlock2(c,:)=BehPhotoM(iAnimal).GrandSummary.Performance(2,:);
    RTBlock2(c,:)=BehPhotoM(iAnimal).GrandSummary.RT(2,:);
    
    
    % Stim-align rastersCor
    SingleAnimalStimTrace= BehPhotoM(iAnimal).GrandSummary.AbsStimRaster; 
    SingleAnimalNormStimTrace = SingleAnimalStimTrace ./ max(max(SingleAnimalStimTrace));
    GrandPopAbsStimResp = SingleAnimalNormStimTrace + GrandPopAbsStimResp ;
    
      % Stim-align rasters(middle stimuli, correct/error)
    SingleAnimalStimTraceMiddleStimCorrect= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterCorrect(1,:); 
    SingleAnimalNormStimTraceMiddleStimCorrect = SingleAnimalStimTraceMiddleStimCorrect ./ max(max(SingleAnimalStimTraceMiddleStimCorrect));
    GrandPopAbsStimRespMiddleStimCorrect = SingleAnimalNormStimTraceMiddleStimCorrect + GrandPopAbsStimRespMiddleStimCorrect ;
    
    SingleAnimalStimTraceMiddleStimError= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterError(1,:); 
    SingleAnimalNormStimTraceMiddleStimError = SingleAnimalStimTraceMiddleStimError ./ max(max(SingleAnimalStimTraceMiddleStimCorrect));
    GrandPopAbsStimRespMiddleStimError = SingleAnimalNormStimTraceMiddleStimError + GrandPopAbsStimRespMiddleStimError ;
    
        % Stim-align rasters(middle stimuli, correct/ small/large reward)
    SingleAnimalStimTraceMiddleStimLargeCorrect= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterLargeCorrect(1,:); 
    SingleAnimalNormStimTraceMiddleStimLargeCorrect = SingleAnimalStimTraceMiddleStimCorrect ./ max(max(SingleAnimalStimTraceMiddleStimLargeCorrect));
    GrandPopAbsStimRespMiddleStimCorLarge = SingleAnimalNormStimTraceMiddleStimLargeCorrect + GrandPopAbsStimRespMiddleStimCorLarge ;
    
    SingleAnimalStimTraceMiddleStimSmallCorrect= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterSmallCorrect(1,:); 
    SingleAnimalNormStimTraceMiddleStimSmallCorrect = SingleAnimalStimTraceMiddleStimSmallCorrect ./ max(max(SingleAnimalStimTraceMiddleStimLargeCorrect));
    GrandPopAbsStimRespMiddleStimCorSmall = SingleAnimalNormStimTraceMiddleStimSmallCorrect + GrandPopAbsStimRespMiddleStimCorSmall ;
    
    
    % tuning curve for different blocks
    SingleAnimalTunningStim= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimNoFold;
    SingleAnimalNormTunningStim = SingleAnimalTunningStim ./ max(max(SingleAnimalTunningStim));
    GrandPopNormBinStimNoFold = SingleAnimalNormTunningStim + GrandPopNormBinStimNoFold ;
        GrandPopNormBinStimNoFold1(c,:)=SingleAnimalNormTunningStim(1,:);
        GrandPopNormBinStimNoFold2(c,:)=SingleAnimalNormTunningStim(2,:);

    
    SingleAnimalTunningReward= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardNoFold;
    SingleAnimalNormTunningRew = SingleAnimalTunningReward ./ max(max(SingleAnimalTunningReward));
    GrandPopNormBinRewardNoFold = SingleAnimalNormTunningRew + GrandPopNormBinRewardNoFold ;
    GrandPopNormBinRewardNoFold1(c,:)=SingleAnimalNormTunningRew(1,:);
    GrandPopNormBinRewardNoFold2(c,:)=SingleAnimalNormTunningRew(2,:);
  
    
    % tuning curve for corr/error
    SingleAnimalTunningStimCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrectErrorNoFold;
    SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError ./ max(max(SingleAnimalTunningStimCorrError));
    GrandPopNormBinStimNoFoldCorrError = SingleAnimalNormTunningStimCorrError + GrandPopNormBinStimNoFoldCorrError ;
    GrandPopNormBinStimNoFoldCorrError1(c,:)=SingleAnimalNormTunningStimCorrError(1,:);
    GrandPopNormBinStimNoFoldCorrError2(c,:)=SingleAnimalNormTunningStimCorrError(2,:);
    
  
    SingleAnimalTunningRewardCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardCorrectErrorNoFold;
    SingleAnimalNormTunningRewardCorrError = SingleAnimalTunningRewardCorrError ./ max(max(SingleAnimalTunningRewardCorrError));
    GrandPopNormBinRewardNoFoldCorrError = SingleAnimalNormTunningRewardCorrError + GrandPopNormBinRewardNoFoldCorrError ;
  GrandPopNormBinRewardNoFoldCorrError1(c,:)=GrandPopNormBinRewardNoFoldCorrError(1,:);
    GrandPopNormBinRewardNoFoldCorrError2(c,:)=GrandPopNormBinRewardNoFoldCorrError(2,:);
  
    SingleAnimalBeep2AwayDACorrErr(1,:)= BehPhotoM(iAnimal).GrandSummary.Beep2DACorr;
    SingleAnimalBeep2AwayDACorrErr(2,:)= BehPhotoM(iAnimal).GrandSummary.BeepAwayDACorr;
    SingleAnimalBeep2AwayDACorrErr(3,:)= BehPhotoM(iAnimal).GrandSummary.Beep2DAErr;
    SingleAnimalBeep2AwayDACorrErr(4,:)= BehPhotoM(iAnimal).GrandSummary.BeepAwayDAErr;
    
    
    SingleAnimalBeep2AwayDACorrErrNorm = SingleAnimalBeep2AwayDACorrErr ./ max(max(SingleAnimalBeep2AwayDACorrErr));
    
    GrandPopBeep2AwayDACorError = SingleAnimalBeep2AwayDACorrErrNorm + GrandPopBeep2AwayDACorError ;
    
    SingleAnimalStim2AwayDACorrErr(1,:)= BehPhotoM(iAnimal).GrandSummary.Stim2DACorr;
    SingleAnimalStim2AwayDACorrErr(2,:)= BehPhotoM(iAnimal).GrandSummary.StimAwayDACorr;
    SingleAnimalStim2AwayDACorrErr(3,:)= BehPhotoM(iAnimal).GrandSummary.Stim2DAErr;
    SingleAnimalStim2AwayDACorrErr(4,:)= BehPhotoM(iAnimal).GrandSummary.StimAwayDAErr;
    
    
    SingleAnimalStim2AwayDACorrErrNorm = SingleAnimalStim2AwayDACorrErr ./ max(max(SingleAnimalStim2AwayDACorrErr));
     
    GrandPopStim2AwayDACorError = SingleAnimalStim2AwayDACorrErrNorm + GrandPopStim2AwayDACorError ;
    
    
    SingleBinStimCorrError(1:2,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrect;
    SingleBinStimCorrError(3:4,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimError;
    
    SingleBinStimCorrErrorNorm = SingleBinStimCorrError ./ max(max(SingleBinStimCorrError));
    
    
    GrandPopStimBin = SingleBinStimCorrErrorNorm + GrandPopStimBin ;
    
    
    SingleAnimalRew2AwayDACorrErr(1,:)= BehPhotoM(iAnimal).GrandSummary.Rew2DACorr;
    SingleAnimalRew2AwayDACorrErr(2,:)= BehPhotoM(iAnimal).GrandSummary.RewAwayDACorr;
    SingleAnimalRew2AwayDACorrErr(3,:)= BehPhotoM(iAnimal).GrandSummary.Rew2DAErr;
    SingleAnimalRew2AwayDACorrErr(4,:)= BehPhotoM(iAnimal).GrandSummary.RewAwayDAErr;
    
    
    SingleAnimalRew2AwayDACorrErrNorm = SingleAnimalRew2AwayDACorrErr ./ max(max(SingleAnimalRew2AwayDACorrErr));
     
    GrandPopRew2AwayDACorError = SingleAnimalRew2AwayDACorrErrNorm + GrandPopRew2AwayDACorError ;
    
    
    SingleBinRewCorrError(1:2,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardCorrect;
    SingleBinRewCorrError(3:4,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardError;
    
    SingleBinRewCorrErrorNorm = SingleBinRewCorrError ./ max(max(SingleBinRewCorrError));
    
    
    GrandPopRewBin = SingleBinRewCorrErrorNorm + GrandPopRewBin ;
    
    
    
    
    
    c=c+1;
end

%%
subplot(5,3,1); hold on
xlabel('Contrast')
ylabel('P(R)')
title( 'Psychometric curves')
set(gca,'TickDir','out','Box','off');

plot(StimAllowed,nanmean(PerBlock1),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
plot(StimAllowed,nanmean(PerBlock2),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
legend('LargeRew@L','LargeRew@R','Location','southeast')


subplot(5,3,2); hold on
xlabel('Contrast')
ylabel('Norm RT')
title( 'Reaction Time')
set(gca,'TickDir','out','Box','off');


plot(StimAllowed,nanmean(RTBlock1),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)

plot(StimAllowed,nanmean(RTBlock2),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)



% subplot(5,3,4); hold on
% plot(StimAllowed,GrandPopNormBinStimNoFold(1,:)./ length(Animals),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinStimNoFold(2,:)./ length(Animals),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% title('Stimulus Align')
% xlabel('Contrast')

    
    subplot(5,3,4); hold on % errorbars reflecting across sessions
errorbar(StimAllowed,GrandPopNormBinStimNoFold(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFold1) ./ sqrt(12),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinStimNoFold(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFold2) ./ sqrt(12),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
title('Stimulus Align')
xlabel('Contrast')

    


%subplot(5,3,5); hold on
% plot(StimAllowed,GrandPopNormBinRewardNoFold(1,:)./ length(Animals),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinRewardNoFold(2,:)./ length(Animals),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% xlabel('Contrast')
% title('Outcome Align')

subplot(5,3,5); hold on % errorbars reflecting across sessions
errorbar(StimAllowed,GrandPopNormBinRewardNoFold(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFold1) ./ sqrt(12),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinRewardNoFold(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFold2) ./ sqrt(12),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Contrast')


% subplot(5,3,7); hold on
% plot(StimAllowed,GrandPopNormBinStimNoFoldCorrError(1,:)./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinStimNoFoldCorrError(2,:)./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% title('Stimulus Align')
% xlabel('Contrast')

subplot(5,3,7); hold on
errorbar(StimAllowed,GrandPopNormBinStimNoFoldCorrError(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFoldCorrError1) ./ sqrt(14),'color','r','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinStimNoFoldCorrError(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFoldCorrError2) ./ sqrt(14),'color','g','LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
title('Stimulus Align')
xlabel('Contrast')




% subplot(5,3,8); hold on
% plot(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(1,:)./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(2,:)./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% xlabel('Contrast')
% title('Outcome Align')

subplot(5,3,8); hold on
<<<<<<< HEAD

errorbar(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFoldCorrError1) ./ sqrt(14),'color','r','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFoldCorrError2) ./ sqrt(14),'color','g','LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Contrast')


subplot(5,3,10); hold on


plot(unique(abs(StimAllowed)),GrandPopStimBin(1,:)./ length(Animals),'--g','LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(abs(StimAllowed)),GrandPopStimBin(2,:)./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)

plot(unique(abs(StimAllowed)),GrandPopStimBin(3,:)./ length(Animals),'--r','LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(abs(StimAllowed)),GrandPopStimBin(4,:)./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
xlabel('Contrast')
title('Stimulus Align')


subplot(5,3,11); hold on


plot(unique(abs(StimAllowed)),GrandPopRewBin(1,:)./ length(Animals),'--g','LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(abs(StimAllowed)),GrandPopRewBin(2,:)./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)

plot(unique(abs(StimAllowed)),GrandPopRewBin(3,:)./ length(Animals),'--r','LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(abs(StimAllowed)),GrandPopRewBin(4,:)./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)

set(gca,'TickDir','out','Box','off');
xlabel('Contrast')
title('Outcome Align')


subplot(5,3,12); hold on

plot(GrandPopStim2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopStim2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopStim2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopStim2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);

xlim([3500 4900])
ylim([-0.3 1])
set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('Stimulus Align')
xlabel('Time (s)')
ylabel('Norm response')
legend('CorLarge','CorSmall','ErrLarge','ErrSmall')




subplot(5,3,13); hold on

plot(GrandPopRew2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopRew2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopRew2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopRew2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);


xlim([3500 4900])
ylim([-0.8 1])
set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Time (s)')



subplot(5,3,14); hold on

plot(GrandPopBeep2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopBeep2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopBeep2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopBeep2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);

xlim([3000 4900])
ylim([-0.3 0.7])
set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('Beep Align')
xlabel('Time (s)')
ylabel('Norm response')

for c = 1:4
    
    
    subplot(5,3,15); hold on
    plot((GrandPopAbsStimResp(c,:) ./ length(Animals)),'color',colorGray(c,:),'LineWidth',2)
  
end


title('Stimulus Align')

xlim([3500 4900])
ylim([-0.3 1])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')



subplot(5,3,3); hold on
title('middleSim,CorrError')
plot(GrandPopAbsStimRespMiddleStimCorrect./ length(Animals),'g','LineWidth',2);

plot(GrandPopAbsStimRespMiddleStimError./ length(Animals),'r','LineWidth',2);

xlim([3500 4900])
ylim([-0.3 0.5])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')


subplot(5,3,6); hold on
title('middleSim,SmallLarge')

plot(GrandPopAbsStimRespMiddleStimCorLarge./ length(Animals),'g','LineWidth',2);

plot(GrandPopAbsStimRespMiddleStimCorSmall./ length(Animals),'--g','LineWidth',2);

xlim([3500 4900])
ylim([-0.3 0.5])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')

    
    