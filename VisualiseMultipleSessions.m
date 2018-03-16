clear all
close all

% give animal name and it will plot data averaged across all session

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)


%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 51;

load('BehPhotoM_Exp23')

RTLimit = 10; % in s, excluding trials with RT longer than this 

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
        RTAv(iBlock,c) = nanmean (TempData(TempData(:,2)==istim,7));
        
        c=c+1;
    end
    
end

subplot(5,3,1); hold on
xlabel('Contrast')
ylabel('P(R)')
title( 'Psychometric curves')
set(gca,'TickDir','out','Box','off');

plot(unique(BehData(:,2))',performance(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(BehData(:,2))',performance(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
legend('LargeRew@L','LargeRew@R','Location','southeast')


subplot(5,3,2); hold on
xlabel('Contrast')
ylabel('Norm RT')
title( 'Reaction Time')
set(gca,'TickDir','out','Box','off');


plot(unique(BehData(:,2))',RTAv(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)

plot(unique(BehData(:,2))',RTAv(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)

%% beep figure


subplot(5,3,9); hold on

plot(nanmean(BeepData),'k','LineWidth',2)

title('Beep Align')

xlim([3500 4900])
%ylim([-0.3 2])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')


%% Stimulus fig panels 

c=1;
for iStim = abzStim
    
    AbsStimRaster(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim, :));
      
    subplot(5,3,13); hold on
    plot((AbsStimRaster(c,:)),'color',colorGray(c,:),'LineWidth',2)
     
    c=c+1;
end

if length(abzStim)==3
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)))
    
elseif  length(abzStim)==4
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)))
    
elseif length(abzStim)==5
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)),num2str(abzStim(5)))
    
end

title('Stimulus Align')

xlim([3500 4900])
%ylim([-0.3 2])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')


NormBinStim = mean(StimData(:,4000:5200),2);

% this figure looks at response as a function of contrast separated for
% blocks with large reward on L or R
c=1;
for iblock= [1 2]
    
    cStim = 1;
    
    for iStimAbs = unique((BehData(:,2)))'
        
        % only correct trials
        PopNormBinStimBlocksNoFold (c,cStim)= nanmean(NormBinStim (BehData(:,9)==1 & BehData(:,2)==iStimAbs & BehData(:,8)==iblock));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end

subplot(5,3,4); hold on

plot(unique(BehData(:,2))',PopNormBinStimBlocksNoFold(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(BehData(:,2))',PopNormBinStimBlocksNoFold(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
title('Stimulus Align')
xlabel('Contrast')
ylabel('Norm response')

% this figure plots rasters at the stimulus time separated based on the
% pendding outcome
subplot(5,3,14); hold on

plot(nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==1,:)),'g','LineWidth',2)
plot(nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==-1,:)),'--g','LineWidth',2)

plot(nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==1,:)),'r','LineWidth',2)
plot(nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==-1,:)),'--r','LineWidth',2)

xlim([3500 4800])
%ylim([-0.3 2])
set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('Stimulus Align')
xlabel('Time (s)')
ylabel('Norm response')

% this figure looks at response at stim time separated based on contrast
% and whether choice was correct or error
c=1;
for CorrError= [0 1]
    
    cStim = 1;
    
    for istim = unique(BehData(:,2))'
        
        PopNormBinStimCorrectErrorNoFold (c,cStim)= nanmean(NormBinStim (BehData(:,9)==CorrError & BehData(:,2)==istim ));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end

subplot(5,3,7); hold on
plot(unique((BehData(:,2)))',PopNormBinStimCorrectErrorNoFold(1,:),'r','LineWidth',2)
plot(unique((BehData(:,2)))',PopNormBinStimCorrectErrorNoFold(2,:),'g','LineWidth',2)

set(gca, 'TickDir','out','Box','off');

title('Stimulus Align')

xlabel('Contrast')
ylabel('Norm response')


% thi figure folds contrasts and looks responses at the stimulus time
% separated based on the reward size and correct/error
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

subplot(5,3,10); hold on
plot(unique(abs(BehData(:,2)))',PopNormBinStimCorrect(1,:),'--g','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinStimCorrect(2,:),'g','LineWidth',2)

plot(unique(abs(BehData(:,2)))',PopNormBinStimError(1,:),'--r','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinStimError(2,:),'r','LineWidth',2)

set(gca, 'TickDir','out','Box','off');

title('Stimulus Align')

xlabel('Contrast')
ylabel('Norm response')



%% Reward figure

% this figure plots rasters at the outcome time separated based on the
% pendding outcome

subplot(5,3,15); hold on

plot(nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==1,:)),'g','LineWidth',2)
plot(nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==-1,:)),'--g','LineWidth',2)

plot(nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==1,:)),'r','LineWidth',2)
plot(nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==-1,:)),'--r','LineWidth',2)

xlim([3500 4900])
%ylim([-2 2])

set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Time (s)')
ylabel('Norm response')

NormBinReward = mean(RewardData(:,4000:5300),2) - mean(RewardData(:,3400:3800),2);

% this figure looks at response as a function of contrast separated for
% blocks with large reward on L or R
c=1;
for iblock= [1 2]
    
    cStim = 1;
    
    for iStimAbs = unique((BehData(:,2)))'
        
        % only correct trials
        PopNormBinRewardBlocksNoFold (c,cStim)= nanmean(NormBinReward (BehData(:,9)==1 & BehData(:,2)==iStimAbs & BehData(:,8)==iblock));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end

subplot(5,3,5); hold on

plot(unique(BehData(:,2))',PopNormBinRewardBlocksNoFold(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(BehData(:,2))',PopNormBinRewardBlocksNoFold(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
title('Outcome Align')
set(gca, 'TickDir','out','Box','off');

xlabel('Contrast')
ylabel('Norm response')


% this figure looks at response at outcome time separated based on contrast
% and whether choice was correct or error
c=1;
for CorrError= [0 1]
    
    cStim = 1;
    
    for istim = unique(BehData(:,2))'
        
        PopNormBinRewardCorrectErrorNoFold (c,cStim)= nanmean(NormBinReward (BehData(:,9)==CorrError & BehData(:,2)==istim ));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end

subplot(5,3,8); hold on
plot(unique((BehData(:,2)))',PopNormBinRewardCorrectErrorNoFold(1,:),'r','LineWidth',2)
plot(unique((BehData(:,2)))',PopNormBinRewardCorrectErrorNoFold(2,:),'g','LineWidth',2)

set(gca, 'TickDir','out','Box','off');

title('Outcome Align')

xlabel('Contrast')
ylabel('Norm response')

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


subplot(5,3,11); hold on
plot(unique(abs(BehData(:,2)))',PopNormBinRewardCorrect(1,:),'--g','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinRewardCorrect(2,:),'g','LineWidth',2)

plot(unique(abs(BehData(:,2)))',PopNormBinRewardError(1,:),'--r','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinRewardError(2,:),'r','LineWidth',2)


set(gca, 'TickDir','out','Box','off');

xlabel('Contrast')
ylabel('Norm response')

title('Outcome Align')
legend('CorSmall','CorLarge','ErrSmall','ErrLarge')

%%
% Grand Summary data of the animal



BehPhotoM(animal_ID).GrandSummary.Performance = performance;

BehPhotoM(animal_ID).GrandSummary.RT = RTAv;

BehPhotoM(animal_ID).GrandSummary.AbsStimRaster = AbsStimRaster;

BehPhotoM(animal_ID).GrandSummary.PopNormBinStimCorrectErrorNoFold=PopNormBinStimCorrectErrorNoFold;

BehPhotoM(animal_ID).GrandSummary.PopNormBinStimNoFold = PopNormBinStimBlocksNoFold;

BehPhotoM(animal_ID).GrandSummary.PopNormBinStimCorrect=PopNormBinStimCorrect;

BehPhotoM(animal_ID).GrandSummary.PopNormBinStimError=PopNormBinStimError;

BehPhotoM(animal_ID).GrandSummary.PopNormBinRewardCorrectErrorNoFold=PopNormBinRewardCorrectErrorNoFold;

BehPhotoM(animal_ID).GrandSummary.PopNormBinRewardNoFold = PopNormBinRewardBlocksNoFold;

BehPhotoM(animal_ID).GrandSummary.PopNormBinRewardCorrect=PopNormBinRewardCorrect;

BehPhotoM(animal_ID).GrandSummary.PopNormBinRewardError=PopNormBinRewardError;


BehPhotoM(animal_ID).GrandSummary.Beep2DACorr = nanmean(BeepData(BehData(:,9)==1 & BehData(:,16)==1,:));
BehPhotoM(animal_ID).GrandSummary.BeepAwayDACorr = nanmean(BeepData(BehData(:,9)==1 & BehData(:,16)==-1,:));

BehPhotoM(animal_ID).GrandSummary.Beep2DAErr = nanmean(BeepData(BehData(:,9)==0 & BehData(:,16)==1,:));
BehPhotoM(animal_ID).GrandSummary.BeepAwayDAErr = nanmean(BeepData(BehData(:,9)==0 & BehData(:,16)==-1,:));


BehPhotoM(animal_ID).GrandSummary.Stim2DACorr = nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==1,:));
BehPhotoM(animal_ID).GrandSummary.StimAwayDACorr = nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==-1,:));

BehPhotoM(animal_ID).GrandSummary.Stim2DAErr = nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==1,:));
BehPhotoM(animal_ID).GrandSummary.StimAwayDAErr = nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==-1,:));

BehPhotoM(animal_ID).GrandSummary.Rew2DACorr = nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==1,:));
BehPhotoM(animal_ID).GrandSummary.RewAwayDACorr = nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==-1,:));

BehPhotoM(animal_ID).GrandSummary.Rew2DAErr = nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==1,:));
BehPhotoM(animal_ID).GrandSummary.RewAwayDAErr = nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==-1,:));








