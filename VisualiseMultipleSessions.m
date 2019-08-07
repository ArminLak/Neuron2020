clear all
close all

% give animal name and it will plot data averaged across all session

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)
% Armin Feb 2018
% Armin July 2018 added bilateral recoding
% Armin Dec 2018, changed to better handle bilateral recordings

%VTA : [48, 50,51, 64, 69]  coresponding to ALK068, 70 and 71, ALK084,ALK085

% DMS 
% Full list: [53, 55,62, 63,64, 68, 70, 71, 72] coresponding to ALK074(Bi), ALK075(Bi), MMM003(Un), ALK083(Bi),
% ALK084(Un), MMM006(Un), MMM008(Un), MMM009(Un), MMM010(Un), 
% Usefull list: : [53, 62, 63, 71,72] 
% 

% NAc : [56, 57,59,66] coresponding to  ALK078(Bi), MMM001(Un), MMM002(Un),% MMM005(UN)

% Exp:
% 7: learning from scratch
% 7Adv: single reward after learning is advanced ( this one is in the
% getSessionList_photpM code.

% 23: double reward size 


% select animal

animal_ID = 57
BrainStrucutre = 'NAc'
ExpID = '23'

save2file = 0; % decide if you want to overwrite GrandSummary or not

load(['BehPhotoM_Exp', ExpID, '_', BrainStrucutre]);


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

colorRed = [ 1 0.8 0.8
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

iter = 0;
if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')
    
    for iSession = sessionz % left hem
        
        TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        BehData = [BehData; TempBehData];
        
        % left
        TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepL;
        
        TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimL;
        
        TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronActionL;
        
        TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardL;
        
        BeepDataL = [BeepDataL;TempBeepData];
        
        StimDataL = [StimDataL;TempStimData];
        
        ActionDataL = [ActionDataL;TempActionData];
        
        RewardDataL = [RewardDataL;TempRewardData];
        
    end
    
    iter = iter + 1;
end


if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardR')
    BehData = [];
    
    for iSession = sessionz % left hem
        
        TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        BehData = [BehData; TempBehData];
        
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
    
    iter = iter + 1;
    
end


RT = BehData(:,10) - BehData(:,13);

toRemove = find ( RT > RTLimit);

toRemove2= find(BehData(:,1) < 20);
toRemove = unique([toRemove; toRemove2]);

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
    
    if iter == 1
        
        if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')
            
            BeepData = BeepDataL;
            StimData = StimDataL;
            ActionData = ActionDataL;
            RewardData = RewardDataL;
            
            
        elseif isfield(BehPhotoM(animal_ID).Session,'NeuronRewardR')
            
            BeepData = BeepDataR;
            StimData = StimDataR;
            ActionData = ActionDataR;
            RewardData = RewardDataR;
            
        end
        
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
    
    % normalisation
   StimTimeDenom=max([ mean(StimData(mintersect(ToLargeR,find(BehData(:,9)),find(BehData(:,2)==max(BehData(:,2)))),:)),...
mean(StimData(mintersect(ToLargeR,find(BehData(:,9)),find(BehData(:,2)==min(BehData(:,2)))),:))]);

    StimData = StimData ./ StimTimeDenom;
ActionData = ActionData ./ StimTimeDenom;
RewardData = RewardData ./ StimTimeDenom;

    
    
    figure; hold on
    
    for iBlock = unique(BehData(:,8))'
        
        TempData = BehData(BehData(:,8)==iBlock,:);
        
        c = 1;
        for istim = unique(BehData(:,2))'
            
            if strcmpi (ExpID,'38')
                performance(iBlock,c) = nanmean (TempData(TempData(:,2)==istim,9));
                
            else
                
                performance(iBlock,c) = nanmean (TempData(TempData(:,2)==istim,3));
            end
            RTAv(iBlock,c) = nanmean (TempData(TempData(:,2)==istim,7));
            
            c=c+1;
        end
        
    end
    
    subplot(6,3,1); hold on
    xlabel('Contrast')
    ylabel('P(R)')
    title( 'Psychometric curves')
    set(gca,'TickDir','out','Box','off');
    
    if strcmpi (ExpID,'23')
        plot(unique(BehData(:,2))',performance(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
        plot(unique(BehData(:,2))',performance(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
        legend('LargeRew@L','LargeRew@R','Location','southeast')
        
    else
        plot(unique(BehData(:,2))',performance(4,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
    end
    
    
    
    subplot(6,3,2); hold on
    xlabel('Contrast')
    ylabel('Norm RT')
    title( 'Reaction Time')
    set(gca,'TickDir','out','Box','off');
    
    if strcmpi (ExpID,'23')
        
        plot(unique(BehData(:,2))',RTAv(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
        plot(unique(BehData(:,2))',RTAv(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
    else
        plot(unique(BehData(:,2))',RTAv(4,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
        
    end
    
    %% beep figure
    
    
    subplot(6,3,9); hold on
    
    plot(nanmean(BeepData(BehData(:,4)==1,:)),'k','LineWidth',2)
        plot(nanmean(BeepData(BehData(:,4)==0,:)),'r','LineWidth',2)

    
    title('Beep Align')
    
    xlim([200 4900])
    %ylim([-0.3 2])
    
    
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
    
    %% Stimulus fig panels
    
    c=1;
    for iStim = abzStim
        
        AbsStimRaster(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim, :));
        AbsActionRaster(c,:)=nanmean(ActionData(abs(BehData(:,2))==iStim, :));
        
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
    if animal_ID == 48
        NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2)- mean(StimData(:,3400:3800),2);
        
       % NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2);
        
        
    elseif animal_ID == 50
        NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,3100:3400),2);
        
        
    elseif animal_ID == 51
        
        %NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2)- mean(StimData(:,3400:3800),2);
        NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,3400:3800),2);
        
        
    elseif animal_ID == 56
        
        NormBinStim = mean(StimData(:,4500:5000),2)- mean(StimData(:,3100:3500),2);
        NormBinStim = mean(StimData(:,4600:5300),2);
        
        
    elseif animal_ID == 57
        
        NormBinStim = mean(StimData(:,4300:5000),2)- mean(StimData(:,3400:3800),2);

        NormBinStim = mean(StimData(:,4600:5300),2)- mean(StimData(:,3400:3800),2);

    elseif animal_ID == 59
        
        NormBinStim = mean(StimData(:,4400:4800),2);
        NormBinStim = mean(StimData(:,4900:5400),2);
        
    
        elseif animal_ID == 64 && strcmp(BrainStrucutre,'DMS')
        
        NormBinStim = mean(StimData(:,4500:5000),2);
        
    elseif animal_ID == 66
        
        NormBinStim = mean(StimData(:,5000:6000),2);
        
                
    elseif animal_ID == 69
        
        NormBinStim = mean(StimData(:,3900:5200),2)- mean(StimData(:,3400:3800),2);
             
        
    elseif animal_ID == 70
        
        NormBinStim = mean(StimData(:,4000:5000),2);
        
       
    else
        
        NormBinStim = mean(StimData(:,4200:5000),2)- mean(StimData(:,3400:3800),2);
        
    end
    
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
    
    subplot(6,3,4); hold on
    
    plot(unique(BehData(:,2))',PopNormBinStimBlocksNoFold(1,:),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
    plot(unique(BehData(:,2))',PopNormBinStimBlocksNoFold(2,:),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
    title('Stimulus Align')
    xlabel('Contrast')
    ylabel('Norm response')
    
    % this figure plots rasters at the stimulus time separated based on the
    % pendding outcome
    subplot(6,3,14); hold on
    
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
    
    subplot(6,3,7); hold on
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
    
    subplot(6,3,10); hold on
    plot(unique(abs(BehData(:,2)))',PopNormBinStimCorrect(1,:),'--g','LineWidth',2)
    plot(unique(abs(BehData(:,2)))',PopNormBinStimCorrect(2,:),'g','LineWidth',2)
    
    plot(unique(abs(BehData(:,2)))',PopNormBinStimError(1,:),'--r','LineWidth',2)
    plot(unique(abs(BehData(:,2)))',PopNormBinStimError(2,:),'r','LineWidth',2)
    
    set(gca, 'TickDir','out','Box','off');
    
    title('Stimulus Align')
    
    xlabel('Contrast')
    ylabel('Norm response')
    
    % make raster for large/small or correct/error
    
    c=1;
    for iStim = abzStim
        
        AbsStimRasterCorrect(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,9)==1, :),1);
        AbsStimRasterError(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,9)==0, :),1);
        
        AbsStimRasterLargeCorrect(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,16)==1 & BehData(:,9)==1, :),1);
        AbsStimRasterSmallCorrect(c,:)=nanmean(StimData(abs(BehData(:,2))==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :),1);
        
        
        AbsActionRasterCorrect(c,:)=nanmean(ActionData(abs(BehData(:,2))==iStim & BehData(:,9)==1, :),1);
        AbsActionRasterError(c,:)=nanmean(ActionData(abs(BehData(:,2))==iStim & BehData(:,9)==0, :),1);
        
        AbsActionRasterLargeCorrect(c,:)=nanmean(ActionData(abs(BehData(:,2))==iStim & BehData(:,16)==1 & BehData(:,9)==1, :),1);
        AbsActionRasterSmallCorrect(c,:)=nanmean(ActionData(abs(BehData(:,2))==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :),1);
        
       
        AbsRewardRasterCorrect(c,:) = nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,9)==1, :),1);
        AbsRewardRasterError(c,:) = nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,9)==0, :),1);
        
        
        AbsRewardRasterLargeCorrect(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,16)==1 & BehData(:,9)==1, :),1);
        AbsRewardRasterSmallCorrect(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :),1);
        
        c=c+1;
    end
    
    
    c=1;
    for iStim = unique((BehData(:,2)))'
        
        StimRasterLargeCorrect(c,:)=nanmean(StimData(BehData(:,2)==iStim & BehData(:,16)==1 & BehData(:,9)==1, :),1);
        StimRasterSmallCorrect(c,:)=nanmean(StimData(BehData(:,2)==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :),1);
        
        StimRasterLargeError(c,:)  =nanmean(StimData(BehData(:,2)==iStim & BehData(:,16)==1 & BehData(:,9)==0, :),1);
        
         ActRasterLargeCorrect(c,:)=nanmean(ActionData(BehData(:,2)==iStim & BehData(:,16)==1 & BehData(:,9)==1, :),1);
        ActRasterSmallCorrect(c,:)=nanmean(ActionData(BehData(:,2)==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :),1);
        
        ActRasterLargeError(c,:)  =nanmean(ActionData(BehData(:,2)==iStim & BehData(:,16)==1 & BehData(:,9)==0, :),1);
        
        RewardRasterLargeCorrect(c,:) = nanmean(RewardData(BehData(:,2)==iStim & BehData(:,16)==1 & BehData(:,9)==1, :),1);
        RewardRasterSmallCorrect(c,:)=nanmean(RewardData(BehData(:,2)==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :),1);
        
        RewardRasterLargeError(c,:)  =nanmean(RewardData(BehData(:,2)==iStim & BehData(:,16)==1 & BehData(:,9)==0, :),1);
        
        c=c+1;
    end
    
    
    %plot((AbsStimRasterCorrect(2,:)),'g','LineWidth',2)
    %plot((AbsStimRasterError(2,:)),'r','LineWidth',2)
    for i = 1:length(abzStim)
        
        subplot(6,3,3); hold on
        
        plot((smooth(AbsStimRasterLargeCorrect(i,:),70)),'color',colorGray(i,:),'LineWidth',2)
        
        plot((smooth(AbsStimRasterSmallCorrect(i,:),70)),'color',colorRed(i,:),'LineWidth',2)
        
     
       % plot((smooth(AbsActionRasterLargeCorrect(i,:),70)),'color',colorGray(i,:),'LineWidth',2)
        
       % plot((smooth(AbsActionRasterSmallCorrect(i,:),70)),'color',colorRed(i,:),'LineWidth',2)
        
    end
    
    title('Stim Align')
    
    xlim([3500 4900])
    ylim([-0.3 2])
    
    
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
    subplot(6,3,6); hold on
    plot((AbsActionRasterLargeCorrect(2,:)),'b','LineWidth',2)
    plot((AbsActionRasterSmallCorrect(2,:)),'--b','LineWidth',2)
    
    
    title('Action Align')
    
    xlim([3500 4900])
    %ylim([-0.3 2])
    
    
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    legend('Large reward', 'Small reward')
    
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
    
    
    for iStim = abzStim
        
        AbsStimRasterCorrectREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,9)==1, :),1);
        AbsStimRasterErrorREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,9)==0, :),1);
        
        AbsStimRasterLargeCorrectREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,16)==1 & BehData(:,9)==1, :),1);
        AbsStimRasterSmallCorrectREw(c,:)=nanmean(RewardData(abs(BehData(:,2))==iStim & BehData(:,16)==-1  & BehData(:,9)==1, :),1);
        
        
        c=c+1;
    end
    
    for i = 1:length(abzStim)
        
        subplot(6,3,15); hold on
        
        plot(smooth(AbsStimRasterCorrectREw(i,:),70),'color',colorGray(i,:),'LineWidth',2)
        
        plot(smooth(AbsStimRasterErrorREw(i,:),70),'color',colorRed(i,:),'LineWidth',2)
        
    end
    
    
    xlim([3500 4900])
    ylim([-2 2])
    
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
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
    
    if animal_ID == 56
        NormBinReward = mean(RewardData(:,4000:4800),2) - mean(RewardData(:,3400:3800),2);
    end
    
    if animal_ID == 57
        NormBinReward = mean(RewardData(:,4000:4800),2) - mean(RewardData(:,3200:3700),2);
    end
    
    if animal_ID == 59
        NormBinReward = mean(RewardData(:,4000:4800),2);
    end
    
        if animal_ID == 66 &&  strcmp(ExpID , '38')
 
        NormBinReward = mean(RewardData(:,4100:4600),2) - mean(RewardData(:,3400:3800),2);
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
    
    subplot(6,3,5); hold on
    
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
    
    subplot(6,3,8); hold on
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
    
    
    subplot(6,3,11); hold on
    plot(unique(abs(BehData(:,2)))',PopNormBinRewardCorrect(1,:),'--g','LineWidth',2)
    plot(unique(abs(BehData(:,2)))',PopNormBinRewardCorrect(2,:),'g','LineWidth',2)
    
    plot(unique(abs(BehData(:,2)))',PopNormBinRewardError(1,:),'--r','LineWidth',2)
    plot(unique(abs(BehData(:,2)))',PopNormBinRewardError(2,:),'r','LineWidth',2)
    
    
    set(gca, 'TickDir','out','Box','off');
    
    xlabel('Contrast')
    ylabel('Norm response')
    
    title('Outcome Align')
    legend('CorSmall','CorLarge','ErrSmall','ErrLarge')
    
    % looking at reward responses separated for short and long RT
    
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
    
    
    subplot(6,3,12); hold on
    plot(unique(abs(BehData(:,2)))',PopNormBinRewardshortRT,'--b','LineWidth',2)
    plot(unique(abs(BehData(:,2)))',PopNormBinRewardlongRT,'b','LineWidth',2)
    
    
    set(gca, 'TickDir','out','Box','off');
    
    xlabel('Contrast')
    ylabel('Norm response')
    xlim([0 0.5])
    title('Outcome Align')
    legend('shortRT','longRT')
    %%
    % conditional psych function ( condition to DA response)
    figure
    
    BehData01 = BehData; % changed response to 0 and 1
    BehData01 (BehData01(:,3)==-1,3) = 0;
    cz = 1;
    for thresholdRange = 80 %[60 70 80 90 95]
        
        DA_threshold = prctile(NormBinStim, thresholdRange);
        
        
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
%    %   session by session (conditional to Da responses at
%    %  the reward time) controlled for the reward size
%     firstTrialOfBlock=[1 ; find(diff(BehData01(:,1)) < 0)+1]';  % define block onset
%     firstTrialOfBlock = [firstTrialOfBlock, size(BehData01,1)];
% 
%      for     Block2test = 1:2 %(or 4 for single value task)
% 
%       performance_BL  = nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
%      
%             performance_SL=nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
%                
%             performance_BR=nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
%             performance_SR=nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
%     
%     for iblockZ = 1:length(firstTrialOfBlock)-1
%         
%         BehDataTemp = BehData01(firstTrialOfBlock(iblockZ):firstTrialOfBlock(iblockZ+1)-1,:);
%         NormBinRewardTemp = NormBinReward(firstTrialOfBlock(iblockZ):firstTrialOfBlock(iblockZ+1)-1,:);
%              
%     DA_threshold_Reward = prctile(NormBinRewardTemp, 65);
%     
%     % neurons
%             index_BigRPE_Right=mintersect(find(BehDataTemp(:,3)==1) , find(BehDataTemp(:,9)==1) , find(NormBinRewardTemp > DA_threshold_Reward),find(BehDataTemp(:,8)==Block2test));
%             index_BigRPE_Right = index_BigRPE_Right+ 1; % go to the next trial
%             
%             
%             index_BigRPE_Left=mintersect(find(BehDataTemp(:,3)==0) , find(BehDataTemp(:,9)==1), find(NormBinRewardTemp > DA_threshold_Reward),find(BehDataTemp(:,8)==Block2test));
%             index_BigRPE_Left = index_BigRPE_Left+ 1; % go to the next trial
%             
%             index_SmallRPE_Right=mintersect(find(BehDataTemp(:,3)==1) , find(BehDataTemp(:,9)==1), find(NormBinRewardTemp < DA_threshold_Reward),find(BehDataTemp(:,8)==Block2test));
%             index_SmallRPE_Right = index_SmallRPE_Right+ 1; % go to the next trial
%             
%             
%             index_SmallRPE_Left=mintersect(find(BehDataTemp(:,3)==0) , find(BehDataTemp(:,9)==1), find(NormBinRewardTemp < DA_threshold_Reward),find(BehDataTemp(:,8)==Block2test));
%             index_SmallRPE_Left = index_SmallRPE_Left+ 1; % go to the next trial
%       
%        
%              c = 1;
%         for istim = unique(BehData(:,2))'
%             
%             indexBL = intersect(find(BehDataTemp(:,2)==istim),index_BigRPE_Left);
%             indexSL = intersect(find(BehDataTemp(:,2)==istim),index_SmallRPE_Left);
%             
%             indexBR = intersect(find(BehDataTemp(:,2)==istim),index_BigRPE_Right);
%             indexSR = intersect(find(BehDataTemp(:,2)==istim),index_SmallRPE_Right);
%            
%             
%             performance_BL(iblockZ,c) = nanmean (BehDataTemp(indexBL,3));
%             performance_SL(iblockZ,c) = nanmean (BehDataTemp(indexSL,3));
%                
%             performance_BR(iblockZ,c) = nanmean (BehDataTemp(indexBR,3));
%             performance_SR(iblockZ,c) = nanmean (BehDataTemp(indexSR,3));
%             
%             
%             c=c+1;
%         end
%             
%     end
%        
%         DAConditionedPerf.mice(1).Block(Block2test).perf(1,:)=nanmean(performance_BL);
%         DAConditionedPerf.mice(1).Block(Block2test).perf(2,:)=nanmean(performance_SL);
%         DAConditionedPerf.mice(1).Block(Block2test).perf(3,:)=nanmean(performance_BR);
%         DAConditionedPerf.mice(1).Block(Block2test).perf(4,:)=nanmean(performance_SR);
% 
%          
%     
%        figure
%         
%         plot(nanmean(performance_BL),'k')
%         hold on
%         plot(nanmean(performance_SL),'-.k')
%         plot(nanmean(performance_BR),'b')
%         plot(nanmean(performance_SR),'-.b')
%             
%      end

     %%
      %  session by session (conditional to Da responses at
   %  the reward time) controlled for  stimulus 
    firstTrialOfBlock=[1 ; find(diff(BehData01(:,1)) < 0)+1]';  % define block onset
    firstTrialOfBlock = [firstTrialOfBlock, size(BehData01,1)];

      
    cc=1;
          for istimPast = unique(BehData(:,2))'   
              
                    performance_ConL  = nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
                    performance_ConS  = nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));

    for iblockZ = 1:length(firstTrialOfBlock)-1
        
        BehDataTemp = BehData01(firstTrialOfBlock(iblockZ):firstTrialOfBlock(iblockZ+1)-1,:);
        NormBinRewardTemp = NormBinReward(firstTrialOfBlock(iblockZ):firstTrialOfBlock(iblockZ+1)-1,:);
             
    DA_threshold_Reward = prctile(NormBinRewardTemp, 65);
           
           index_ConL=mintersect(find(BehDataTemp(:,9)==1) ,find(NormBinRewardTemp > DA_threshold_Reward),find(BehDataTemp(:,2)==istimPast));
           
            index_ConL = index_ConL+ 1; % go to the next trial
            
            
          index_ConS=mintersect(find(BehDataTemp(:,9)==1) ,find(NormBinRewardTemp < DA_threshold_Reward),find(BehDataTemp(:,2)==istimPast));
           
            index_ConS = index_ConS+ 1; % go to the next trial
            
            
             
            
             c = 1;
        for istim = unique(BehData(:,2))'
              
            indexL = intersect(find(BehDataTemp(:,2)==istim),index_ConL);
          
            performance_ConL(iblockZ,c) = nanmean (BehDataTemp(indexL,3));
            
            indexS = intersect(find(BehDataTemp(:,2)==istim),index_ConS);
          
            performance_ConS(iblockZ,c) = nanmean (BehDataTemp(indexS,3));
            
            
            c=c+1;
        end
            
       
          end
          DAConditionedPerf.mice(1).Block(cc).perf(1,:)=nanmean(performance_ConL);
          DAConditionedPerf.mice(1).Block(cc).perf(2,:)=nanmean(performance_ConS);
          
        
        cc=cc+1;
          end
%%          
   %%
%       %  session by session (conditional to Da responses at
%    %  the reward time) controlled for  stimulus and reward size
%     firstTrialOfBlock=[1 ; find(diff(BehData01(:,1)) < 0)+1]';  % define block onset
%     firstTrialOfBlock = [firstTrialOfBlock, size(BehData01,1)];
% 
%       
%     cc=1;
%           for istimPast = unique(BehData(:,2))'   
%               
%                     performance_ConLL  = nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
%                     performance_ConLS  = nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
%                     performance_ConRL  = nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
%                     performance_ConRS  = nan(length(firstTrialOfBlock),length(unique(BehData(:,2))'));
% 
%     for iblockZ = 1:length(firstTrialOfBlock)-1
%         
%         BehDataTemp = BehData01(firstTrialOfBlock(iblockZ):firstTrialOfBlock(iblockZ+1)-1,:);
%         NormBinRewardTemp = NormBinReward(firstTrialOfBlock(iblockZ):firstTrialOfBlock(iblockZ+1)-1,:);
%              
%     DA_threshold_Reward = prctile(NormBinRewardTemp, 65);
%            
%     
%            index_ConLL=mintersect(find(BehDataTemp(:,9)==1) ,   find(BehDataTemp(:,3)==-1) ,...
%                find(NormBinRewardTemp > DA_threshold_Reward),find(BehDataTemp(:,2)==istimPast));
%            
%             index_ConLL = index_ConLL+ 1; % go to the next trial
%             
%             index_ConLS=mintersect(find(BehDataTemp(:,9)==1) ,   find(BehDataTemp(:,3)==-1) ,...
%                find(NormBinRewardTemp < DA_threshold_Reward),find(BehDataTemp(:,2)==istimPast));
%            
%             index_ConLS = index_ConLS+ 1; % go to the next trial
%             
%          index_ConRL=mintersect(find(BehDataTemp(:,9)==1) ,   find(BehDataTemp(:,3)==1) ,...
%                find(NormBinRewardTemp > DA_threshold_Reward),find(BehDataTemp(:,2)==istimPast));
%            
%             index_ConRL = index_ConRL+ 1; % go to the next trial
%             
%             index_ConRS=mintersect(find(BehDataTemp(:,9)==1) ,   find(BehDataTemp(:,3)==-1) ,...
%                find(NormBinRewardTemp < DA_threshold_Reward),find(BehDataTemp(:,2)==istimPast));
%            
%             index_ConRS = index_ConRS+ 1; % go to the next trial
%             
%            
%             
%             
%             
%              c = 1;
%         for istim = unique(BehData(:,2))'
%               
%             indexL = intersect(find(BehDataTemp(:,2)==istim),index_ConL);
%           
%             performance_ConL(iblockZ,c) = nanmean (BehDataTemp(indexL,3));
%             
%             indexS = intersect(find(BehDataTemp(:,2)==istim),index_ConS);
%           
%             performance_ConS(iblockZ,c) = nanmean (BehDataTemp(indexS,3));
%             
%             
%             c=c+1;
%         end
%             
%        
%           end
%           DAConditionedPerf.mice(1).Block(cc).perf(1,:)=nanmean(performance_ConL);
%           DAConditionedPerf.mice(1).Block(cc).perf(2,:)=nanmean(performance_ConS);
%           
%         
%         cc=cc+1;
% end
          
    %%        
            
            
     % this part is too look at all data at once (not separating sessions)  , it is again for the performance conditional to DA responses at reward time     
            
%      performance_BL=nan(1,length(unique(BehData(:,2))'));
%      
%             performance_SL=nan(1,length(unique(BehData(:,2))'));
%                
%             performance_BR=nan(1,length(unique(BehData(:,2))'));
%             performance_SR=nan(1,length(unique(BehData(:,2))'));
%              
%     DA_threshold_Reward = prctile(NormBinReward, 65);
%     Block2test = 4;
%     
%     % neurons
%             index_BigRPE_Right=mintersect(find(BehData(:,3)==1) , find(BehData(:,9)==1) , find(NormBinReward > DA_threshold_Reward),find(BehData(:,8)==Block2test));
%             index_BigRPE_Right = index_BigRPE_Right+ 1; % go to the next trial
%             
%             
%             index_BigRPE_Left=mintersect(find(BehData(:,3)==-1) , find(BehData(:,9)==1), find(NormBinReward > DA_threshold_Reward),find(BehData(:,8)==Block2test));
%             index_BigRPE_Left = index_BigRPE_Left+ 1; % go to the next trial
%             
%             index_SmallRPE_Right=mintersect(find(BehData(:,3)==1) , find(BehData(:,9)==1), find(NormBinReward < DA_threshold_Reward),find(BehData(:,8)==Block2test));
%             index_SmallRPE_Right = index_SmallRPE_Right+ 1; % go to the next trial
%             
%             
%             index_SmallRPE_Left=mintersect(find(BehData(:,3)==-1) , find(BehData(:,9)==1), find(NormBinReward < DA_threshold_Reward),find(BehData(:,8)==Block2test));
%             index_SmallRPE_Left = index_SmallRPE_Left+ 1; % go to the next trial
%          
%   
% % beh
% %             index_BigRPE_Right=mintersect(find(BehData(:,3)==1) , find(BehData(:,9)==1) , find(BehData(:,2)==0),find(BehData(:,8)==Block2test));
% %             index_BigRPE_Right = index_BigRPE_Right+ 1; % go to the next trial
% %             
% %             
% %             index_BigRPE_Left=mintersect(find(BehData(:,3)==-1) , find(BehData(:,9)==1),  find(BehData(:,2)==0),find(BehData(:,8)==Block2test));
% %             index_BigRPE_Left = index_BigRPE_Left+ 1; % go to the next trial
% %             
% %             index_SmallRPE_Right=mintersect(find(BehData(:,3)==1) , find(BehData(:,9)==1),  find(BehData(:,2)>=0.5),find(BehData(:,8)==Block2test));
% %             index_SmallRPE_Right = index_SmallRPE_Right+ 1; % go to the next trial
% %             
% %             
% %             index_SmallRPE_Left=mintersect(find(BehData(:,3)==-1) , find(BehData(:,9)==1),  find(BehData(:,2)<=-0.5),find(BehData(:,8)==Block2test));
% %             index_SmallRPE_Left = index_SmallRPE_Left+ 1; % go to the next trial
% 
%             
%             
%      c = 1;
%         for istim = unique(BehData(:,2))'
%             
%             indexBL = intersect(find(BehData01(:,2)==istim),index_BigRPE_Left);
%             indexSL = intersect(find(BehData01(:,2)==istim),index_SmallRPE_Left);
%             
%             indexBR = intersect(find(BehData01(:,2)==istim),index_BigRPE_Right);
%             indexSR = intersect(find(BehData01(:,2)==istim),index_SmallRPE_Right);
%            
%             
%             performance_BL(c) = nanmean (BehData01(indexBL,3));
%             performance_SL(c) = nanmean (BehData01(indexSL,3));
%                
%             performance_BR(c) = nanmean (BehData01(indexBR,3));
%             performance_SR(c) = nanmean (BehData01(indexSR,3));
%             
%             
%             c=c+1;
%         end
%     
%         DAConditionedPerf.mice(3).Block4 (1,:)=performance_BL;
%         DAConditionedPerf.mice(3).Block4 (2,:)=performance_SL;
%         DAConditionedPerf.mice(3).Block4 (3,:)=performance_BR;
%         DAConditionedPerf.mice(3).Block4 (4,:)=performance_SR;
%         
%     
%         figure
%         
%         plot(performance_BL,'k')
%         hold on
%         plot(performance_SL,'-.k')
%         plot(performance_BR,'b')
%         plot(performance_SR,'-.b')
        
     %%
    
    if (HemIter ==1 && isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')) || (iter == 2 && HemIter ==1)
        BehPhotoM(animal_ID).GrandSummaryL.Performance = performance;
        BehPhotoM(animal_ID).GrandSummaryL.RT = RTAv;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRaster = AbsStimRaster;
        BehPhotoM(animal_ID).GrandSummaryL.AbsActionRaster = AbsActionRaster;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinStimCorrectErrorNoFold=PopNormBinStimCorrectErrorNoFold;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinStimNoFold = PopNormBinStimBlocksNoFold;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinStimCorrect=PopNormBinStimCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinStimError=PopNormBinStimError;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinRewardCorrectErrorNoFold=PopNormBinRewardCorrectErrorNoFold;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinRewardNoFold = PopNormBinRewardBlocksNoFold;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinRewardCorrect=PopNormBinRewardCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.PopNormBinRewardError=PopNormBinRewardError;
        BehPhotoM(animal_ID).GrandSummaryL.Beep2DACorr = nanmean(BeepData(BehData(:,9)==1 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryL.BeepAwayDACorr = nanmean(BeepData(BehData(:,9)==1 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryL.Beep2DAErr = nanmean(BeepData(BehData(:,9)==0 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryL.BeepAwayDAErr = nanmean(BeepData(BehData(:,9)==0 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryL.Stim2DACorr = nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryL.StimAwayDACorr = nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryL.Stim2DAErr = nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryL.StimAwayDAErr = nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryL.Rew2DACorr = nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryL.RewAwayDACorr = nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryL.Rew2DAErr = nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryL.RewAwayDAErr = nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterCorrect=AbsStimRasterCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterError=AbsStimRasterError;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterLargeCorrect=AbsStimRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterSmallCorrect = AbsStimRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.AbsActionRasterCorrect=AbsActionRasterCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.AbsActionRasterError=AbsActionRasterError;
        BehPhotoM(animal_ID).GrandSummaryL.AbsActionRasterLargeCorrect=AbsActionRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.AbsActionRasterSmallCorrect = AbsActionRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterCorrectREw=AbsStimRasterCorrectREw;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterErrorREw=AbsStimRasterErrorREw;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterLargeCorrectREw=AbsStimRasterLargeCorrectREw;
        BehPhotoM(animal_ID).GrandSummaryL.AbsStimRasterSmallCorrectREw = AbsStimRasterSmallCorrectREw;
        BehPhotoM(animal_ID).GrandSummaryL.performance_lowDA = performance_lowDA;
        BehPhotoM(animal_ID).GrandSummaryL.performance_highDA = performance_highDA;
        
        
        BehPhotoM(animal_ID).GrandSummaryL.StimRasterLargeCorrect = StimRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.StimRasterSmallCorrect = StimRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.StimRasterLargeError = StimRasterLargeError;
        
        
        BehPhotoM(animal_ID).GrandSummaryL.ActionRasterLargeCorrect = ActRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.ActionRasterSmallCorrect = ActRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.ActionRasterLargeError = ActRasterLargeError;
        
        
        BehPhotoM(animal_ID).GrandSummaryL.RewardRasterLargeCorrect = RewardRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.RewardRasterSmallCorrect = RewardRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryL.RewardRasterLargeError = RewardRasterLargeError;
        
        
    end
    
    if (HemIter ==1 && isfield(BehPhotoM(animal_ID).Session,'NeuronRewardR')) || (iter == 2 && HemIter ==2)
        
        BehPhotoM(animal_ID).GrandSummaryR.Performance = performance;
        BehPhotoM(animal_ID).GrandSummaryR.RT = RTAv;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRaster = AbsStimRaster;
        BehPhotoM(animal_ID).GrandSummaryR.AbsActionRaster = AbsActionRaster;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinStimCorrectErrorNoFold=PopNormBinStimCorrectErrorNoFold;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinStimNoFold = PopNormBinStimBlocksNoFold;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinStimCorrect=PopNormBinStimCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinStimError=PopNormBinStimError;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinRewardCorrectErrorNoFold=PopNormBinRewardCorrectErrorNoFold;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinRewardNoFold = PopNormBinRewardBlocksNoFold;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinRewardCorrect=PopNormBinRewardCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.PopNormBinRewardError=PopNormBinRewardError;
        BehPhotoM(animal_ID).GrandSummaryR.Beep2DACorr = nanmean(BeepData(BehData(:,9)==1 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryR.BeepAwayDACorr = nanmean(BeepData(BehData(:,9)==1 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryR.Beep2DAErr = nanmean(BeepData(BehData(:,9)==0 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryR.BeepAwayDAErr = nanmean(BeepData(BehData(:,9)==0 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryR.Stim2DACorr = nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryR.StimAwayDACorr = nanmean(StimData(BehData(:,9)==1 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryR.Stim2DAErr = nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryR.StimAwayDAErr = nanmean(StimData(BehData(:,9)==0 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryR.Rew2DACorr = nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryR.RewAwayDACorr = nanmean(RewardData(BehData(:,9)==1 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryR.Rew2DAErr = nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==1,:));
        BehPhotoM(animal_ID).GrandSummaryR.RewAwayDAErr = nanmean(RewardData(BehData(:,9)==0 & BehData(:,16)==-1,:));
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterCorrect=AbsStimRasterCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterError=AbsStimRasterError;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterLargeCorrect=AbsStimRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterSmallCorrect = AbsStimRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.AbsActionRasterCorrect=AbsActionRasterCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.AbsActionRasterError=AbsActionRasterError;
        BehPhotoM(animal_ID).GrandSummaryR.AbsActionRasterLargeCorrect=AbsActionRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.AbsActionRasterSmallCorrect = AbsActionRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterCorrectREw=AbsStimRasterCorrectREw;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterErrorREw=AbsStimRasterErrorREw;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterLargeCorrectREw=AbsStimRasterLargeCorrectREw;
        BehPhotoM(animal_ID).GrandSummaryR.AbsStimRasterSmallCorrectREw = AbsStimRasterSmallCorrectREw;
        BehPhotoM(animal_ID).GrandSummaryR.performance_lowDA = performance_lowDA;
        BehPhotoM(animal_ID).GrandSummaryR.performance_highDA = performance_highDA;
        
        BehPhotoM(animal_ID).GrandSummaryR.StimRasterLargeCorrect = StimRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.StimRasterSmallCorrect = StimRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.StimRasterLargeError = StimRasterLargeError;
        
            BehPhotoM(animal_ID).GrandSummaryR.ActionRasterLargeCorrect = ActRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.ActionRasterSmallCorrect = ActRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.ActionRasterLargeError = ActRasterLargeError;
        
        
            BehPhotoM(animal_ID).GrandSummaryR.RewardRasterLargeCorrect = RewardRasterLargeCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.RewardRasterSmallCorrect = RewardRasterSmallCorrect;
        BehPhotoM(animal_ID).GrandSummaryR.RewardRasterLargeError = RewardRasterLargeError;
        
    end
    
    figure %stim, action, and reward responses separated by ipsi/contra and contrast 


for istim = 1:length(unique(BehData(:,2)))

        subplot(3,3,1)
    plot(smooth(StimRasterLargeCorrect(istim,:),70))
    hold on
         title('Large Correct') %Stim-aligned
    xlim([3500 4900])
    ylim([-0.2 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
        subplot(3,3,2)
    plot(smooth(StimRasterSmallCorrect(istim,:),80))
    hold on
        title('Small Correct') %Stim-aligned
    xlim([3500 4900])
    ylim([-0.2 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
    
        subplot(3,3,3)
    plot(smooth(StimRasterLargeError(istim,:),120))
    hold on
        title('Large Error') %Stim-aligned
    xlim([3500 4900])
    ylim([-0.2 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')

    
        subplot(3,3,4)
    plot(smooth(ActRasterLargeCorrect(istim,:),70))
    hold on
         title('Large Correct') %Action-aligned
    xlim([3000 4400])
    ylim([-0.2 1])
    set(gca, 'XTick', [3000, 3700, 4400]);
    set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
           subplot(3,3,5)
    plot(smooth(ActRasterSmallCorrect(istim,:),70))
    hold on
        title('Small Correct') %Action-aligned
    xlim([3000 4400])
    ylim([-0.2 1])
    set(gca, 'XTick', [3000, 3700, 4400]);
    set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
           subplot(3,3,6)
    plot(smooth(ActRasterLargeError(istim,:),70))
    hold on
         title('Large Error') %Action-aligned
    xlim([3000 4400])
    ylim([-0.2 1])
    set(gca, 'XTick', [3000, 3700, 4400]);
    set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
    
            subplot(3,3,7)
    plot(smooth(RewardRasterLargeCorrect(istim,:),70))
    hold on
         title('Large Correct') %Reward-aligned
    xlim([3500 4900])
    ylim([-0.2 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
           subplot(3,3,8)
    plot(smooth(RewardRasterSmallCorrect(istim,:),70))
    hold on
        title('Small Correct') %Reward-aligned
    xlim([3500 4900])
    ylim([-0.2 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
           subplot(3,3,9)
    plot(smooth(RewardRasterLargeError(istim,:),70))
    hold on
         title('Large Error') %Outcome-aligned
    xlim([3500 4900])
    ylim([-0.2 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
  
    
end

   figure %figure 4
    
   subplot(3,3,1)
    plot(smooth(mean(StimRasterLargeCorrect(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(StimRasterLargeCorrect(5:7,:)),70),'b')
    hold on
         title('Large Correct') %Stim-aligned
    xlim([3500 4900])
    ylim([-0.5 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')   
    
      subplot(3,3,2)
    plot(smooth(mean(StimRasterSmallCorrect(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(StimRasterSmallCorrect(5:7,:)),70),'b')
    hold on
         title('Small Correct') %Stim-aligned
    xlim([3500 4900])
    ylim([-0.5 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')   
    
    
      subplot(3,3,3)
    plot(smooth(mean(StimRasterLargeError(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(StimRasterLargeError(5:7,:)),70),'b')
    hold on
         title('Small Error') %Stim-aligned
    xlim([3500 4900])
    ylim([-0.5 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
     subplot(3,3,4) %action aligned
    plot(smooth(mean(ActRasterLargeCorrect(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(ActRasterLargeCorrect(5:7,:)),70),'b')
    hold on
             title('act->') 
        xlim([3000 4400])
    ylim([-0.5 1])
    set(gca, 'XTick', [3000, 3700, 4400]);
    set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
      subplot(3,3,5) %action aligned 
    plot(smooth(mean(ActRasterSmallCorrect(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(ActRasterSmallCorrect(5:7,:)),70),'b')
    hold on
          xlim([3000 4400])
    ylim([-0.5 1])
    set(gca, 'XTick', [3000, 3700, 4400]);
    set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
    
      subplot(3,3,6) %action aligned 
    plot(smooth(mean(ActRasterLargeError(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(ActRasterLargeError(5:7,:)),70),'b')
    hold on
  xlim([3000 4400])
    ylim([-0.5 1])
    set(gca, 'XTick', [3000, 3700, 4400]);
    set(gca, 'XTickLabel', {'-.7','0','0.7'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
    
    
       subplot(3,3,7) % large correct outcome
    plot(smooth(mean(RewardRasterLargeCorrect(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(RewardRasterLargeCorrect(5:7,:)),70),'b')
    hold on
         title('rwd->') 
    xlim([3500 4900])
    ylim([-0.5 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')   
    
      subplot(3,3,8) %small correct outcome
    plot(smooth(mean(RewardRasterSmallCorrect(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(RewardRasterSmallCorrect(5:7,:)),70),'b')
    hold on
    xlim([3500 4900])
    ylim([-0.5 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')   
    
    
      subplot(3,3,9) %large error outcome
    plot(smooth(mean(RewardRasterLargeError(1:3,:)),70),'k')
    hold on
    plot(smooth(mean(RewardRasterLargeError(5:7,:)),70),'b')
    hold on
    xlim([3500 4900])
    ylim([-0.5 1])
    set(gca, 'XTick', [3700, 4300, 4900]);
    set(gca, 'XTickLabel', {'0','0.6','1.2'},'TickDir','out','Box','off');
    xlabel('Time (s)')
    ylabel('Norm response')
    
end
%%



figure

scatter(NormBinStim(intersect(ToLargeR,find(BehData(:,9)==1))),NormBinReward(intersect(ToLargeR,find(BehData(:,9)==1))))

hold on

scatter(NormBinStim(intersect(ToSmallR,find(BehData(:,9)==1))),NormBinReward(intersect(ToSmallR,find(BehData(:,9)==1))))

scatter(NormBinStim(find(BehData(:,9)==0)),NormBinReward(find(BehData(:,9)==0)))



%%
if save2file

if strcmpi(getComputerName,'zopamine2')
    cd ('C:\Users\Armin\Dropbox\Work\UCL\Science\Analysis Code\PhotoM')
elseif strcmpi(getComputerName, 'zebrafish')
    cd ('C:\Users\morga\Documents\MATLAB')
end

save(['BehPhotoM_Exp', ExpID, '_', BrainStrucutre], 'BehPhotoM', '-v7.3');
end
