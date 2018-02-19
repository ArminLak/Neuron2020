clear all
load('BehPhotoM_Exp23')

animal_ID = 48;

color = [
1 0 0         % red
1 0.5 0       % orange
1 1 0         % yellow
0.5 1 0.5     % light green
0 1 1         % light blue
0  0.5 1      % medium blue 
0 0 1];       % blue
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

BehData(RT > 2,:) = [];

BeepData(RT > 2,:) = [];
StimData(RT > 2,:) = [];
RewardData(RT > 2,:) = [];


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

plot(performance(1,:),'k')

plot(performance(2,:),'r')

subplot(4,2,2); hold on

plot(RT(1,:),'k')

plot(RT(2,:),'r')

%%
subplot(4,2,3); hold on

plot((nanmean((BeepData))))


%%
color1 = [1 0.5 0.5 ; 1 0 0];
color2 = [0.0 0.5 1 ; 0 0 1];


c=1;
for iAction = [-1 1]
    
    Block1(c,:)=nanmean(StimData(BehData(:,8)==1& BehData(:,3)==iAction,:));
    Block2(c,:)=nanmean(StimData(BehData(:,8)==2& BehData(:,3)==iAction,:));
   
    subplot(4,2,4); hold on
plot((Block1(c,:)),'color',color1(c,:)) 
plot((Block2(c,:)),'color',color2(c,:)) 



    c=c+1
end

%%

c=1;
for iStim = unique(BehData(:,2))'
    
    StimCor(c,:)=nanmean(StimData(BehData(:,2)==iStim & BehData(:,9)==1,:));
    StimErr(c,:)=nanmean(StimData(BehData(:,2)==iStim & BehData(:,9)==0,:));
    
    RewCor(c,:)=nanmean(StimData(BehData(:,2)==iStim & BehData(:,9)==1,:));
    RewErr(c,:)=nanmean(StimData(BehData(:,2)==iStim & BehData(:,9)==0,:));
    
    subplot(4,2,5); hold on
plot((StimCor(c,:)),'color',color(c,:)) 

    subplot(4,2,6); hold on

plot((StimErr(c,:)),'color',color(c,:))


    c=c+1;
end

%%
