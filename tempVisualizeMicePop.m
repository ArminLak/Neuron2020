clear all
close all

load('GrandSummaryExp23_2.mat')

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
GrandPopNormBinStimNoFold = zeros(2,7);
GrandPopStim2AwayDACorError = zeros(4,7100);

GrandPopBeep2AwayDACorError = zeros(4,7100);

GrandPopRew2AwayDACorError = zeros(4,7100);



GrandPopStimBin = zeros(4,4);

GrandPopRewBin = zeros(4,4);


Animals = [48 50]

for iAnimal = Animals
    
    SingleAnimalStimTrace= GrandSummary(iAnimal).AbsStimResp;
    
    SingleAnimalNormStimTrace = SingleAnimalStimTrace ./ max(max(SingleAnimalStimTrace));
    
    GrandPopAbsStimResp = SingleAnimalNormStimTrace + GrandPopAbsStimResp ;
    
    
    SingleAnimalTunning= GrandSummary(iAnimal).PopNormBinStimNoFold;
    
    SingleAnimalNormTunning = SingleAnimalTunning ./ max(max(SingleAnimalTunning));
    
    GrandPopNormBinStimNoFold = SingleAnimalNormTunning + GrandPopNormBinStimNoFold ;
    
    
    SingleAnimalBeep2AwayDACorrErr(1,:)= GrandSummary(iAnimal).Beep2DACorr;
    SingleAnimalBeep2AwayDACorrErr(2,:)= GrandSummary(iAnimal).BeepAwayDACorr;
    SingleAnimalBeep2AwayDACorrErr(3,:)= GrandSummary(iAnimal).Beep2DAErr;
    SingleAnimalBeep2AwayDACorrErr(4,:)= GrandSummary(iAnimal).BeepAwayDAErr;
    
    
    SingleAnimalBeep2AwayDACorrErrNorm = SingleAnimalBeep2AwayDACorrErr ./ max(max(SingleAnimalBeep2AwayDACorrErr));
    
    GrandPopBeep2AwayDACorError = SingleAnimalBeep2AwayDACorrErrNorm + GrandPopBeep2AwayDACorError ;
    
    SingleAnimalStim2AwayDACorrErr(1,:)= GrandSummary(iAnimal).Stim2DACorr;
    SingleAnimalStim2AwayDACorrErr(2,:)= GrandSummary(iAnimal).StimAwayDACorr;
    SingleAnimalStim2AwayDACorrErr(3,:)= GrandSummary(iAnimal).Stim2DAErr;
    SingleAnimalStim2AwayDACorrErr(4,:)= GrandSummary(iAnimal).StimAwayDAErr;
    
    
    SingleAnimalStim2AwayDACorrErrNorm = SingleAnimalStim2AwayDACorrErr ./ max(max(SingleAnimalStim2AwayDACorrErr));
    
    
    GrandPopStim2AwayDACorError = SingleAnimalStim2AwayDACorrErrNorm + GrandPopStim2AwayDACorError ;
    
    
    SingleBinStimCorrError(1:2,:)= GrandSummary(iAnimal).PopNormBinStimCorrect;
    SingleBinStimCorrError(3:4,:)= GrandSummary(iAnimal).PopNormBinStimError;
    
    SingleBinStimCorrErrorNorm = SingleBinStimCorrError ./ max(max(SingleBinStimCorrError));
    
    
    GrandPopStimBin = SingleBinStimCorrErrorNorm + GrandPopStimBin ;
    
    
    
    
    SingleAnimalRew2AwayDACorrErr(1,:)= GrandSummary(iAnimal).Rew2DACorr;
    SingleAnimalRew2AwayDACorrErr(2,:)= GrandSummary(iAnimal).RewAwayDACorr;
    SingleAnimalRew2AwayDACorrErr(3,:)= GrandSummary(iAnimal).Rew2DAErr;
    SingleAnimalRew2AwayDACorrErr(4,:)= GrandSummary(iAnimal).RewAwayDAErr;
    
    
    SingleAnimalRew2AwayDACorrErrNorm = SingleAnimalRew2AwayDACorrErr ./ max(max(SingleAnimalRew2AwayDACorrErr));
    
    
    
    GrandPopRew2AwayDACorError = SingleAnimalRew2AwayDACorrErrNorm + GrandPopRew2AwayDACorError ;
    
    
    SingleBinRewCorrError(1:2,:)= GrandSummary(iAnimal).PopNormBinRewardCorrect;
    SingleBinRewCorrError(3:4,:)= GrandSummary(iAnimal).PopNormBinRewardError;
    
    SingleBinRewCorrErrorNorm = SingleBinRewCorrError ./ max(max(SingleBinRewCorrError));
    
    
    GrandPopRewBin = SingleBinRewCorrErrorNorm + GrandPopRewBin ;
    
    
end


for c = 1:4
    
    
    subplot(4,2,1); hold on
    plot((GrandPopAbsStimResp(c,:) ./ length(Animals)),'color',colorGray(c,:),'LineWidth',2)
    
    
end


title('C) Stimulus Align')

xlim([3500 4900])
ylim([-0.3 1])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')


subplot(4,2,2); hold on
plot(StimAllowed,GrandPopNormBinStimNoFold(1,:)./ length(Animals),'color',[0.5 0.2 0.1],'LineWidth',2)
plot(StimAllowed,GrandPopNormBinStimNoFold(2,:)./ length(Animals),'color',[1 0.6 0.2],'LineWidth',2)
set(gca,'TickDir','out','Box','off');


subplot(4,2,3); hold on

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


subplot(4,2,5); hold on

plot(GrandPopStim2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopStim2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopStim2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopStim2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);

xlim([3500 4900])
ylim([-0.3 1])
set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('C) Stimulus Align')
xlabel('Time (s)')
ylabel('Norm response')



subplot(4,2,6); hold on

plot(GrandPopRew2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopRew2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopRew2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopRew2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);


xlim([3500 4900])
ylim([-0.8 1])
set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('C) Outcome Align')
xlabel('Time (s)')
ylabel('Norm response')

subplot(4,2,7); hold on


plot(unique(abs(StimAllowed)),GrandPopStimBin(1,:)./ length(Animals),'--g','LineWidth',2)
plot(unique(abs(StimAllowed)),GrandPopStimBin(2,:)./ length(Animals),'g','LineWidth',2)

plot(unique(abs(StimAllowed)),GrandPopStimBin(3,:)./ length(Animals),'--r','LineWidth',2)
plot(unique(abs(StimAllowed)),GrandPopStimBin(4,:)./ length(Animals),'r','LineWidth',2)
set(gca,'TickDir','out','Box','off');


subplot(4,2,8); hold on


plot(unique(abs(StimAllowed)),GrandPopRewBin(1,:)./ length(Animals),'--g','LineWidth',2)
plot(unique(abs(StimAllowed)),GrandPopRewBin(2,:)./ length(Animals),'g','LineWidth',2)

plot(unique(abs(StimAllowed)),GrandPopRewBin(3,:)./ length(Animals),'--r','LineWidth',2)
plot(unique(abs(StimAllowed)),GrandPopRewBin(4,:)./ length(Animals),'r','LineWidth',2)

set(gca,'TickDir','out','Box','off');




%%


c=1;
for iStim = abzStim
    
    AbsStimResp(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim, :));
    
    
    subplot(4,2,3); hold on
    plot((AbsStimResp(c,:)),'color',colorGray(c,:),'LineWidth',2)
    
    
    
    c=c+1;
end

if length(abzStim)==3
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)))
    
elseif  length(abzStim)==4
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)))
    
elseif length(abzStim)==5
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)),num2str(abzStim(5)))
    
end

title('C) Stimulus Align')

xlim([3500 4900])
ylim([-0.3 2])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')


%% Stimulus figure

subplot(4,2,5); hold on

plot(nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==1,:)),'g','LineWidth',2)
plot(nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==-1,:)),'--g','LineWidth',2)

plot(nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==1,:)),'r','LineWidth',2)
plot(nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==-1,:)),'--r','LineWidth',2)

xlim([3500 4800])
ylim([-0.3 2])
set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('C) Stimulus Align')
xlabel('Time (s)')
ylabel('Norm response')


NormBinStim = mean(StimData(:,4200:5200),2);

c=1;
for iblock= [1 2]
    
    cStim = 1;
    
    for iStimAbs = unique((BehData(:,2)))'
        
        PopNormBinStimNoFold (c,cStim)= nanmean(NormBinStim (BehData(:,9)==1 & BehData(:,2)==iStimAbs & BehData(:,8)==iblock));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end

subplot(4,2,4); hold on

plot(unique(BehData(:,2))',PopNormBinStimNoFold(1,:),'color',[0.5 0.2 0.1],'LineWidth',2)
plot(unique(BehData(:,2))',PopNormBinStimNoFold(2,:),'color',[1 0.6 0.2],'LineWidth',2)


%%

c=1;
for RewardSize= [-1 1]
    
    cStim = 1;
    
    for iStimAbs = abzStim
        
        PopNormBinStimCorrect (c,cStim)= nanmean(NormBinStim (BehData(:,9)==1 & abs(BehData(:,2))==iStimAbs & BehData(:,16)==RewardSize));
        PopNormBinStimError (c,cStim)= nanmean(NormBinStim (BehData(:,9)==0 & abs(BehData(:,2))==iStimAbs & BehData(:,16)==RewardSize));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end


subplot(4,2,6); hold on
plot(unique(abs(BehData(:,2)))',PopNormBinStimCorrect(1,:),'--g','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinStimCorrect(2,:),'g','LineWidth',2)

plot(unique(abs(BehData(:,2)))',PopNormBinStimError(1,:),'--r','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinStimError(2,:),'r','LineWidth',2)

set(gca, 'TickDir','out','Box','off');

title('E) Stimulus Align')

xlabel('Contrast')
ylabel('Norm response')


%% Reward figure



subplot(4,2,7); hold on

plot(nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==1,:)),'g','LineWidth',2)
plot(nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==-1,:)),'--g','LineWidth',2)

plot(nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==1,:)),'r','LineWidth',2)
plot(nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==-1,:)),'--r','LineWidth',2)

xlim([3500 4900])
ylim([-2 2])

set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('D) Reward Align')
xlabel('Time (s)')
ylabel('Norm response')

NormBinReward = mean(RewardData(:,4200:5300),2) - mean(RewardData(:,3500:3800),2);

c=1;
for RewardSize= [-1 1]
    
    cStim = 1;
    
    for iStimAbs = abzStim
        
        PopNormBinRewardCorrect (c,cStim)= nanmean(NormBinReward (BehData(:,9)==1 & abs(BehData(:,2))==iStimAbs & BehData(:,16)==RewardSize));
        PopNormBinRewardError (c,cStim)= nanmean(NormBinReward (BehData(:,9)==0 & abs(BehData(:,2))==iStimAbs & BehData(:,16)==RewardSize));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end


subplot(4,2,8); hold on
plot(unique(abs(BehData(:,2)))',PopNormBinRewardCorrect(1,:),'--g','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinRewardCorrect(2,:),'g','LineWidth',2)

plot(unique(abs(BehData(:,2)))',PopNormBinRewardError(1,:),'--r','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinRewardError(2,:),'r','LineWidth',2)


set(gca, 'TickDir','out','Box','off');

xlabel('Contrast')
ylabel('Norm response')

title('F) Outcome Align')
legend('CorSmall','CorLarge','ErrSmall','ErrLarge')




