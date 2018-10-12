clear all
close all

% give animal name and it will plot data averaged across all session

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)
% Armin Feb 2018
% Armin July 2018 added bilateral recoding



%VTA : [48, 50,51]  coresponding to ALK068, 70 and 71
% DMS : [53, 55] coresponding to ALK074(Bi), ALK075(Bi)
% NAc : [56, 57,59] coresponding to  ALK078(Bi), MMM001(Un), MMM002(Un)


% select animal
animal_ID = 53

% select database
% load('BehPhotoM_Exp23')

%load('BehPhotoM_Exp23_NAc')

 load('BehPhotoM_Exp23_DMS')

% define implant
Implant = 'Bi' 


if strcmp(Implant,'Un')
    ChanNum =1;
elseif strcmp(Implant,'Bi')
    ChanNum =[1 2];
end;


RTLimit = 6; % in s, excluding trials with RT longer than this

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

colorGray = [ 0.8 0.8 0.8
    0.6 0.6 0.6
    0.4 0.4 0.4
    0 0 0
     1 0.8 0.8
    1 0.6 0.6
    1 0.4 0.4
    1 0 0];

colorRed = [ 1 0.8 0.8
    1 0.6 0.6
    1 0.4 0.4
    1 0 0];

colorRed = [ 0.8 0.8 0.8
    0.6 0.6 0.6
    0.4 0.4 0.4
    0 0 0
     1 0.8 0.8
    1 0.6 0.6
    1 0.4 0.4
    1 0 0];


%%
BehData = [];
BeepData = [];
StimData = [];
ActionData = [];
RewardData = [];

BeepDataL   = [];
StimDataL   = [];
ActionDataL = [];
RewardDataL = [];

BeepDataR   = [];
StimDataR   = [];
ActionDataR = [];
RewardDataR = [];


sessionz = 1:length(BehPhotoM(animal_ID).Session);


if strcmp(Implant,'Un')

    iter = 1;
for iSession = sessionz
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStim;
    
    TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronAction;

    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;
    
    BeepData = [BeepData;TempBeepData];
    
    StimData = [StimData;TempStimData];
    
    ActionData = [ActionData;TempActionData];

    RewardData = [RewardData;TempRewardData];
    
    BehData = [BehData; TempBehData];
    
    
end

elseif strcmp(Implant,'Bi')
    
    iter = 2;
    
   for iSession = sessionz % left hem
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    BehData = [BehData; TempBehData];
    
    % left
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStim;
    
    TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronAction;

    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;
    
    BeepDataL = [BeepDataL;TempBeepData];
    
    StimDataL = [StimDataL;TempStimData];
    
    ActionDataL = [ActionDataL;TempActionData];

    RewardDataL = [RewardDataL;TempRewardData];
     
    % right
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepR;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
    
    TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronActionR;

    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardR;
    
    BeepDataR = [BeepDataR;TempBeepData];
    
    StimDataR = [StimDataR;TempStimData];
    
    ActionDataR = [ActionDataR;TempActionData];

    RewardDataR = [RewardDataR;TempRewardData];
    
    
end
    
    
end

        
RT = BehData(:,10) - BehData(:,13);

toRemove = find ( RT > RTLimit);
BehData(toRemove,:) = [];

for HemIter = 1:iter

    if iter == 2 && HemIter ==1
      
    BeepData = BeepDataL;  
    StimData = StimDataL;   
    ActionData = ActionDataL;
    RewardData = RewardDataL; 
    
    end
       
     if iter == 2 && HemIter ==2
      
    BeepData = BeepDataR;    
    StimData = StimDataR;   
    ActionData = ActionDataR;
    RewardData = RewardDataR; 
    
     end

BeepData(toRemove,:) = [];
StimData(toRemove,:) = [];
ActionData(toRemove,:) = [];
RewardData(toRemove,:) = [];


ToLargeR = find((BehData(:,3)==-1 & BehData(:,8)==1)  | ...
    (BehData(:,3)==1 & BehData(:,8)==2));

BehData(ToLargeR,16)=1;

ToSmallR = setdiff(1:size(BehData,1),ToLargeR)';

BehData(ToSmallR, 16)=-1;

abzStim = unique(abs(BehData(:,2)))';

% to plots rasters for all stimuli
abzStim = unique(BehData(:,2))';


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

subplot(6,3,1); hold on
xlabel('Contrast')
ylabel('P(R)')
title( 'Psychometric curves')
set(gca,'TickDir','out','Box','off');

plot(unique(BehData(:,2))',performance(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(BehData(:,2))',performance(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
legend('LargeRew@L','LargeRew@R','Location','southeast')


subplot(6,3,2); hold on
xlabel('Contrast')
ylabel('Norm RT')
title( 'Reaction Time')
set(gca,'TickDir','out','Box','off');


plot(unique(BehData(:,2))',RTAv(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)

plot(unique(BehData(:,2))',RTAv(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)



%% Stimulus fig panels

c=1;
for iStim = abzStim
    
    AbsStimRaster(c,:)=nanmean(StimData((BehData(:,2))==iStim, :));
    AbsActionRaster(c,:)=nanmean(ActionData((BehData(:,2))==iStim, :));
    
    subplot(6,3,13); hold on
    plot((AbsStimRaster(c,:)),'color',colorGray(c,:),'LineWidth',2)
    
    subplot(6,3,16); hold on
    plot((AbsActionRaster(c,:)),'color',colorGray(c,:),'LineWidth',2)
    
    
    c=c+1;
end

if length(abzStim)==3
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)))
    
elseif  length(abzStim)==4
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)))
    
elseif length(abzStim)==5
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)),num2str(abzStim(5)))

    elseif length(abzStim)==6
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)),num2str(abzStim(5)),num2str(abzStim(6)))

    elseif length(abzStim)==7
    legend (num2str(abzStim(1)),num2str(abzStim(2)),num2str(abzStim(3)),num2str(abzStim(4)),...
        num2str(abzStim(5)),num2str(abzStim(6)),num2str(abzStim(7)))

    
end
subplot(6,3,13);
title('Stimulus Align')

xlim([3500 4900])


set(gca, 'XTick', [3700, 4300, 4900]);
set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')

subplot(6,3,16);
title('Action Align')

xlim([3000 4400])


set(gca, 'XTick', [3000, 3700, 4400]);
set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')

%%


end



%%
% Grand Summary data of the animal





