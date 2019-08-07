clear all
close all

% give animal name and it will plot data averaged across all session -
% separating for left and right stimuli and left and right hemispheres
% figure 1 = left hemisphere =) enjoy 

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)
% adapted from Armin 2018
% Morgane October 2018


% September 2019 Morgane inspected this code, it is redundant to
% VisualiseFullStimRaster_MultMice



%VTA : [48, 50,51]  coresponding to ALK068, 70 and 71
% DMS : [53, 55] coresponding to ALK074(Bi), ALK075(Bi)
% NAc : [56, 57,59] coresponding to  ALK078(Bi), MMM001(Un), MMM002(Un)


% select animal
animal_ID = 56

% select database
% load('BehPhotoM_Exp23')

load('BehPhotoM_Exp23_NAc')

% load('BehPhotoM_Exp23_DMS')

% define implant
Implant = 'Un' 


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

stimcolors = [
    1 0.2 0.6
    1 0.4 0.7
    1 0.6 0.8
    1 0.8 0.9
%     0.9 0.8 1
    0.8 0.6 1
    0.7 0.4 1
    0.6 0.2 1];

errorRed = [1 0.2 0.2];


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

if animal_ID == 55 % excluding one bad session from ALK075
    sessionz = 2:length(BehPhotoM(animal_ID).Session);
elseif  animal_ID == 56 % excluding from ALK078
    sessionz = [1 2 5 6 7 8 9 10];
elseif animal_ID == 57 % excluding from MMM001
    sessionz = [1:7, 9:12];
elseif animal_ID == 59 % excluding from MMM002
    sessionz = [1, 4:length(BehPhotoM(animal_ID).Session)];
end

if strcmp(Implant,'Un')

    iter = 1;
for iSession = sessionz
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepL;
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimL;
    TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronActionL;
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardL;
    
    elseif isfield(BehPhotoM(animal_ID).Session,'NeuronRewardR')
       
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepR;
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
    TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronActionR;
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardR;
    end
    
    
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

abzStim = unique(abs(BehData(:,2)))';
% to plots rasters for all stimuli
abzStim = unique(BehData(:,2))';



BeepDataNorm = BeepData;
StimDataNorm = StimData;
ActionDataNorm = ActionData;
RewardDataNorm = RewardData;

c = 1;

ToLargeR = find((BehData(:,3)==-1 & BehData(:,8)==1)  | ...
    (BehData(:,3)==1 & BehData(:,8)==2));

BehData(ToLargeR,16)=1;

ToSmallR = setdiff(1:size(BehData,1),ToLargeR)';

BehData(ToSmallR, 16)=-1;




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

subplot(4,2,1); hold on
xlabel('Contrast')
ylabel('P(R)')
title( 'Psychometric curves')
set(gca,'TickDir','out','Box','off');

plot(unique(BehData(:,2))',performance(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
plot(unique(BehData(:,2))',performance(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
legend('LargeRew@L','LargeRew@R','Location','southeast')


subplot(4,2,2); hold on
xlabel('Contrast')
ylabel('Norm RT')
title( 'Reaction Time')
set(gca,'TickDir','out','Box','off');


plot(unique(BehData(:,2))',RTAv(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)

plot(unique(BehData(:,2))',RTAv(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)



%% Stimulus fig panels
% figure; hold on 
StimNormRaster = [];
ActionNormRaster = [];
RewardNormRaster = [];
c=1;
for iStim = abzStim
    
    StimNormRaster(c,:)=nanmean(StimDataNorm((BehData(:,2))==iStim, :));
    ActionNormRaster(c,:)=nanmean(ActionDataNorm((BehData(:,2))==iStim, :));
    RewardNormRaster(c,:) = nanmean(RewardDataNorm((BehData(:,2))==iStim, :));
    
    subplot(4,2,3); hold on
    plot((StimNormRaster(c,:)./max(max(StimNormRaster))),'color',stimcolors(c,:),'LineWidth',2)
    
    subplot(4,2,5); hold on
    plot((ActionNormRaster(c,:)./max(max(ActionNormRaster))),'color',stimcolors(c,:),'LineWidth',2)
    
    subplot(4,2,7); hold on
    plot((RewardNormRaster(c,:)./max(max(RewardNormRaster))),'color',stimcolors(c,:),'LineWidth',2)
    
    
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
subplot(4, 2,3);
title('Stimulus Align')

xlim([3500 4500])


set(gca, 'XTick', [3700, 4500]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')

subplot(4, 2, 5);
title('Action Align')

xlim([3000 4400])


set(gca, 'XTick', [3000, 3700, 4400]);
set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')

subplot(4, 2, 7);
title('Outcome Align')

xlim([3500 4500])


set(gca, 'XTick', [3700, 4500]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');
xlabel('Time (s)')
ylabel('Norm response')


%% Error/Correct , small/large , left/right

subplot(4, 2, 4)
icorrect = find(BehData(:,9)==1);
ierror = find(BehData(:,9)==0);
irightchoice = find(BehData(:,3)==1);
ileftchoice = find(BehData(:,3)==-1);
ileftblock = find(BehData(:,8)==1);
irightblock = find(BehData(:,8)==2);
plotnum = 4;


TempRasterNorm = [];

M = {StimDataNorm, ActionDataNorm, RewardDataNorm};

for i = 1:3
    data = M{i};
    TempRasterNorm(1,:) = nanmean(data(mintersect(icorrect, ileftchoice, ileftblock), :)); %correct large L choice
    TempRasterNorm(2,:) = nanmean(data(mintersect(icorrect, irightchoice, irightblock), :)); %correct large R choice
    TempRasterNorm(3,:) = nanmean(data(mintersect(icorrect, ileftchoice, irightblock), :)); %correct small L choice
    TempRasterNorm(4,:) = nanmean(data(mintersect(icorrect, irightchoice, ileftblock), :)); %correct small R choice
    TempRasterNorm(5,:) = nanmean(data(mintersect(ierror, ileftchoice), :)); %error L choice
    TempRasterNorm(6,:) = nanmean(data(mintersect(ierror, irightchoice), :)); %error R choice
    i = i+1;
    
    subplot(4, 2, plotnum); hold on
        plot(TempRasterNorm(1,:)./max(max(TempRasterNorm)), 'LineWidth', 2, 'color', stimcolors(1,:)) %correct large L
        plot(TempRasterNorm(2,:)./max(max(TempRasterNorm)), 'LineWidth', 2, 'color', stimcolors(end,:)) %correct large R
        plot(TempRasterNorm(3,:)./max(max(TempRasterNorm)), '--', 'LineWidth', 2, 'color', stimcolors(1,:)) %correct small L
        plot(TempRasterNorm(4,:)./max(max(TempRasterNorm)), '--', 'LineWidth', 2, 'color', stimcolors(end,:)) %correct small R
        plot(TempRasterNorm(5,:)./max(max(TempRasterNorm)), 'LineWidth', 2, 'color', [0.4 0 0.2]) %error L
        plot(TempRasterNorm(6,:)./max(max(TempRasterNorm)), 'LineWidth', 2, 'color', [0.2 0 0.4]) %error R 

    plotnum = plotnum +2;

end


subplot(4, 2,4);
    title('Stimulus Align')
    xlim([3500 4500])
    set(gca, 'XTick', [3700, 4500]);
    set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')

subplot(4, 2, 6);
    title('Action Align')
    xlim([3000 4400])
    set(gca, 'XTick', [3000, 3700, 4400]);
    set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')

subplot(4, 2, 8);
    title('Outcome Align')
    xlim([3500 4500])
    set(gca, 'XTick', [3700, 4500]);
    set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')

legend('Large + L', 'Large + R', 'Small + L', 'Small + R', 'Error + L', 'Error + R')

end



%%
% Grand Summary data of the animal





