%clear all
close all

load('BehPhotoM_Exp23')

animal_ID = 48;

RTLimit = 10;


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
%%
BehData = [];
BeepData = [];
StimData = [];
RewardData = [];

for iSession = 1:length(BehPhotoM(animal_ID).Session)
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStim;
    
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;
    
    BeepData = [BeepData;TempBeepData];
    
    StimData = [StimData;TempStimData];
    
    RewardData = [RewardData;TempRewardData];
    
    BehData = [BehData; TempBehData];
    
    
end

RT = BehData(:,10) - BehData(:,13);

BehData(RT > RTLimit,:) = [];

BeepData(RT > RTLimit,:) = [];
StimData(RT > RTLimit,:) = [];
RewardData(RT > RTLimit,:) = [];




ToLargeR = find((BehData(:,3)==-1 & BehData(:,8)==1)  | ...
    (BehData(:,3)==1 & BehData(:,8)==2));

BehData(ToLargeR,16)=1;

ToSmallR = setdiff(1:size(BehData,1),ToLargeR)';

BehData(ToSmallR, 16)=-1;

abzStim = unique(abs(BehData(:,2)))';


figure; hold on

for iBlock = [1 2]
    
    TempData = BehData(BehData(:,8)==iBlock,:);
    
    c = 1;
    for istim = unique(BehData(:,2))'
        
        performance(iBlock,c) = nanmean (TempData(TempData(:,2)==istim,3));
        RT(iBlock,c) = nanmean (TempData(TempData(:,2)==istim,7));
        
        c=c+1;
    end
    
    
end

subplot(4,2,1); hold on
xlabel('Contrast')
ylabel('P(R)')
title( 'A) Psychometric curves')

plot(unique(BehData(:,2))',performance(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',10)
plot(unique(BehData(:,2))',performance(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',10)
legend('LargeRew@L','LargeRew@R','Location','southeast')


subplot(4,2,2); hold on
xlabel('Contrast')
ylabel('Norm RT')
title( 'B) Reaction Time')


plot(unique(BehData(:,2))',RT(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',10)

plot(unique(BehData(:,2))',RT(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',10)





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

    



    
