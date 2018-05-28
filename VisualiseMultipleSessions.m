clear all
close all

% give animal name and it will plot data averaged across all session

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)


%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 51
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
    0 0 0];

colorRed = [ 1 0.8 0.8
    1 0.6 0.6
    1 0.4 0.4
    1 0 0];


%%
BehData = [];
BeepData = [];
StimData = [];
RewardData = [];

sessionz = 1:length(BehPhotoM(animal_ID).Session);
%if animal_ID == 51 
%sessionz = [1 2 3 4 5 8]; end % sessions used for ALK071 Exp 23 [1 2 3 4 5 8]

for iSession = sessionz
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStim;
    
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;
    
    BeepData = [BeepData;TempBeepData];
    
    StimData = [StimData;TempStimData];
    
    RewardData = [RewardData;TempRewardData];
    
    BehData = [BehData; TempBehData];
    
    
end

% BeepData = diff(BeepData')';
% 
% StimData = diff(StimData')';
% RewardData = diff(RewardData')';



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


if animal_ID == 48 
NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2)- mean(StimData(:,3400:3800),2);
end

if animal_ID == 50 
NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,3100:3400),2);
end

NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2)- mean(StimData(:,3400:3800),2);



% for derivative
%NormBinStim = mean(StimData(:,3800:4200),2);

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

% for 2 middle stimuli, make raster for large/small or correct/error
c=1;
middleStim =abzStim(2:3);
middleStim =abzStim(1:4);

for iStim = middleStim
    
    AbsStimRasterCorrect(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,9)==1, :));
    AbsStimRasterError(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,9)==0, :));
    
    AbsStimRasterLargeCorrect(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,16)==1 & BehData(:,9)==1, :));
    AbsStimRasterSmallCorrect(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :));
     
    c=c+1;
end

%plot((AbsStimRasterCorrect(2,:)),'g','LineWidth',2)
%plot((AbsStimRasterError(2,:)),'r','LineWidth',2)
for i = 1:4
    
    subplot(5,3,3); hold on

plot((smooth(AbsStimRasterCorrect(i,:),70)),'color',colorGray(i,:),'LineWidth',2)

plot((smooth(AbsStimRasterError(i,:),70)),'color',colorRed(i,:),'LineWidth',2)

end

title('Stimulus Align')

%xlim([3500 4900])
%ylim([-0.3 2])


%set(gca, 'XTick', [3700, 4300, 4900]);
%set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')
subplot(5,3,6); hold on
plot((AbsStimRasterLargeCorrect(2,:)),'b','LineWidth',2)
plot((AbsStimRasterSmallCorrect(2,:)),'--b','LineWidth',2)


title('Stimulus Align')

xlim([3500 4900])
%ylim([-0.3 2])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')


%% Reward figure

% this figure plots rasters at the outcome time separated based on the
% pendding outcome

% subplot(5,3,15); hold on
% 
% plot(nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==1,:)),'g','LineWidth',2)
% plot(nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==-1,:)),'--g','LineWidth',2)
% 
% plot(nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==1,:)),'r','LineWidth',2)
% plot(nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==-1,:)),'--r','LineWidth',2)

c=1;


for iStim = abzStim(1:4)
    
    AbsStimRasterCorrectREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,9)==1, :));
    AbsStimRasterErrorREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,9)==0, :));
    
    AbsStimRasterLargeCorrectREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,16)==1 & BehData(:,9)==1, :));
    AbsStimRasterSmallCorrectREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :));
    
    
    c=c+1;
end
    
for i = 1:4
    
    subplot(5,3,15); hold on

plot(smooth(AbsStimRasterCorrectREw(i,:),70),'color',colorGray(i,:),'LineWidth',2)

plot(smooth(AbsStimRasterErrorREw(i,:),70),'color',colorRed(i,:),'LineWidth',2)

end


%xlim([3500 4900])
%ylim([-2 2])

%set(gca, 'XTick', [3700, 4300, 4900]);
%set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Time (s)')
ylabel('Norm response')

NormBinReward = mean(RewardData(:,4000:5300),2) - mean(RewardData(:,3400:3800),2);

if animal_ID == 48
NormBinReward = mean(RewardData(:,4100:5000),2) - mean(RewardData(:,3800:3900),2);
end

if animal_ID == 50
NormBinReward = mean(RewardData(:,4100:5000),2) - mean(RewardData(:,3800:3900),2);
end


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

% looking at reward responses separated for short and ling RT

RT = BehData(:,10) - BehData(:,13);
c=1;
for RewardSize= 1
    
    cStim = 1;
    
    for iStimAbs = abzStim
        
        shortRT_threshold = prctile(RT(find(abs(BehData(:,2))==iStimAbs)),20);
        longRT_threshold = prctile(RT(find(abs(BehData(:,2))==iStimAbs)),80);
        
        
        PopNormBinRewardshortRT (cStim)= nanmean(NormBinReward (BehData(:,9)==1 & abs(BehData(:,2))==iStimAbs & BehData(:,16)==RewardSize & RT<shortRT_threshold ));
        PopNormBinRewardlongRT (cStim)= nanmean(NormBinReward (BehData(:,9)==1 & abs(BehData(:,2))==iStimAbs & BehData(:,16)==RewardSize & RT>longRT_threshold));
        
        cStim = cStim + 1;
    end
    
    c=c+1;
end


subplot(5,3,12); hold on
plot(unique(abs(BehData(:,2)))',PopNormBinRewardshortRT,'--b','LineWidth',2)
plot(unique(abs(BehData(:,2)))',PopNormBinRewardlongRT,'b','LineWidth',2)


set(gca, 'TickDir','out','Box','off');

xlabel('Contrast')
ylabel('Norm response')

title('Outcome Align')
legend('shortRT','longRT')
%%
% conditional psych function ( condition to DA response)
figure

BehData01 = BehData; % changed response to 0 and 1
BehData01 (BehData01(:,3)==-1,3) = 0;
cz = 1;
for thresholdRange = 80 %[60 70 80 90 95]

    DA_threshold = percentile(NormBinStim, thresholdRange);
    
    
    c = 1;
    for istim = unique(BehData01(:,2))'
        
        performance_lowDA(c) = nanmean (BehData01(BehData01(:,2)==istim & NormBinStim < DA_threshold,3));
        performance_highDA(c) = nanmean (BehData01(BehData01(:,2)==istim & NormBinStim > DA_threshold,3));
        
        c=c+1;
    end
    
    subplot(1,6,cz)
    plot(performance_lowDA,'k')
    hold on
    plot(performance_highDA,'r')
    
    cz = cz +1;

end

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

BehPhotoM(animal_ID).GrandSummary.AbsStimRasterCorrect=AbsStimRasterCorrect;
BehPhotoM(animal_ID).GrandSummary.AbsStimRasterError=AbsStimRasterError;

BehPhotoM(animal_ID).GrandSummary.AbsStimRasterLargeCorrect=AbsStimRasterLargeCorrect;
BehPhotoM(animal_ID).GrandSummary.AbsStimRasterSmallCorrect = AbsStimRasterSmallCorrect;

BehPhotoM(animal_ID).GrandSummary.AbsStimRasterCorrectREw=AbsStimRasterCorrectREw;
BehPhotoM(animal_ID).GrandSummary.AbsStimRasterErrorREw=AbsStimRasterErrorREw;

BehPhotoM(animal_ID).GrandSummary.AbsStimRasterLargeCorrectREw=AbsStimRasterLargeCorrectREw;
BehPhotoM(animal_ID).GrandSummary.AbsStimRasterSmallCorrectREw = AbsStimRasterSmallCorrectREw;

BehPhotoM(animal_ID).GrandSummary.performance_lowDA = performance_lowDA;
BehPhotoM(animal_ID).GrandSummary.performance_highDA = performance_highDA;


% where we save the data
%cd('\\zubjects.cortexlab.net\Lab\Share\Lak\Morgane')



