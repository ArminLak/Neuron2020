% creates 2 rasters for contra and ipsi (per animal or per session) ordered by
% Error/Small/Large and by RTs. Works for unihem and bihem animals
% Morgane April 2019

% can you see this? yes

% for information about normalisation, check end of this code. 

% close all

clear all


% ----- enter reqs --------------------------------------------------------

%DMS 
 if exist('brain_region', 'var') && strcmp(brain_region, 'DMS')
     clearvars -except BehPhotoM
 else clear all
 end
 Animals = [53, 62, 63, 71,72]
 brain_region = 'DMS'

% NAC
%  if exist('brain_region', 'var') && strcmp(brain_region, 'NAc')
%      clearvars -except BehPhotoM
%   else clear all
%  end
% Animals = [56 57 59 66]
% brain_region = 'NAc'

% VTA: 
%  if exist('brain_region', 'var') && strcmp(brain_region, 'VTA')
%      clearvars -except BehPhotoM
%  else clear all
%  end
% Animals = [48 50 51 64]
% brain_region = 'VTA'



exp_ID = '23'
plotRasters = 0;

concatenate = 1; %show all trials across all sessions for an animal
IpsiContra = 1; %separate rasters based on ipsi / contra? 

stim_2_plot = 0.5; %should be positive
smooth_factor = 100;

LargeErrorSmallColor = [0.15 0.4 0
    1 0 0
    0.15 0.75 0
    ];

RT_min = 0.2;
RT_max = 3; % reaction time range to include

OT_min = 0.2;
OT_max = 3; % outcome time range to include 

% ------------------ start stop times for task events in second ------------------

sStart = -0.1; %stimulus
sStop = 3;

aStart = -0.6; %action
aStop = 0.2;

rStart = -0.2; %reward
rStop = 0.6;

downSample = 1200; % sampling rate AFTER DOWNSAMPLING in makeDataSummaryMatrix
eventOnset = 3700;  % in the saved matrix, the event occurs in this column

% -------------------------------------------------------------------------

if ~exist('BehPhotoM', 'var')
    load(['BehPhotoM_Exp',exp_ID,'_',brain_region,'.mat'])
end

fullStim = [-stim_2_plot stim_2_plot];

animal_count = 0;

for animal_ID = Animals
    animal_count = animal_count +1;
nSessions = 1:length(BehPhotoM(animal_ID).Session);

for iSession = nSessions
    
    if animal_count ==1 && iSession ==1
        leftStimTrials = [];
        rightStimTrials = [];
        
        LargeStimDataIpsi        = [];
        LargeStimDataContra      = [];
        SmallStimDataIpsi        = [];
        SmallStimDataContra      = [];
        ErrorStimDataContra      = [];
        ErrorStimDataIpsi        = [];
        
        LargeBehDataIpsi         = [];
        LargeBehDataContra       = [];
        SmallBehDataIpsi         = [];
        SmallBehDataContra       = [];
        ErrorBehDataIpsi         = [];
        ErrorBehDataContra       = [];
        
    end
    

    % this is interval between outcome and stimulus 
    RTs = BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,10) - BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,13);
    OTs = BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,14) - BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,13);
    
    RTExcludeTrials = unique([(find(RTs < RT_min)) ; (find(RTs >RT_max)) ; (find(OTs < OT_min)) ; (find(OTs > OT_max))]);

    leftStimTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2)== -(stim_2_plot)),RTExcludeTrials);
    rightStimTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2)==stim_2_plot), RTExcludeTrials);
   
    correctTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1), RTExcludeTrials);
    
    largeRewardTrials = setdiff(unique([mintersect(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,8)==1), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,3))==-1); mintersect(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,8)==2), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,3)==1))]), RTExcludeTrials);
    
    smallRewardTrials = setdiff(unique([mintersect(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,8)==1), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,3))==1); mintersect(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,8)==2), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1), ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,3)==-1))]), RTExcludeTrials);

    errorTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==0), RTExcludeTrials);
    
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimL')%load hemispheric deltaFoverF. add Action/Reward to the following two 'if's if desired. 
        LargeStimDataIpsi        = [LargeStimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimL(intersect(leftStimTrials,largeRewardTrials),:)];
        LargeStimDataContra      = [LargeStimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimL(intersect(rightStimTrials,largeRewardTrials),:)];
        SmallStimDataIpsi        = [SmallStimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimL(intersect(leftStimTrials,smallRewardTrials),:)];
        SmallStimDataContra      = [SmallStimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimL(intersect(rightStimTrials,smallRewardTrials),:)];
        ErrorStimDataContra      = [ErrorStimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimL(intersect(rightStimTrials,errorTrials),:)];
        ErrorStimDataIpsi        = [ErrorStimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimL(intersect(leftStimTrials,errorTrials),:)];
        
        LargeBehDataIpsi         = [LargeBehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(leftStimTrials,largeRewardTrials),:)];
        LargeBehDataContra       = [LargeBehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(rightStimTrials,largeRewardTrials),:)];
        SmallBehDataIpsi         = [SmallBehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(leftStimTrials,smallRewardTrials),:)];
        SmallBehDataContra       = [SmallBehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(rightStimTrials,smallRewardTrials),:)];
        ErrorBehDataIpsi         = [ErrorBehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(leftStimTrials,errorTrials),:)];
        ErrorBehDataContra       = [ErrorBehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(rightStimTrials,errorTrials),:)];
    end
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimR') 
        LargeStimDataIpsi        = [LargeStimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimR(intersect(rightStimTrials,largeRewardTrials),:)];
        LargeStimDataContra      = [LargeStimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimR(intersect(leftStimTrials,largeRewardTrials),:)];
        SmallStimDataIpsi        = [SmallStimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimR(intersect(rightStimTrials,smallRewardTrials),:)];
        SmallStimDataContra      = [SmallStimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimR(intersect(leftStimTrials,smallRewardTrials),:)];
        ErrorStimDataContra      = [ErrorStimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimR(intersect(leftStimTrials,errorTrials),:)];
        ErrorStimDataIpsi        = [ErrorStimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimR(intersect(rightStimTrials,errorTrials),:)];
        
        LargeBehDataIpsi         = [LargeBehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(rightStimTrials,largeRewardTrials),:)];
        LargeBehDataContra       = [LargeBehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(leftStimTrials,largeRewardTrials),:)];
        SmallBehDataIpsi         = [SmallBehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(rightStimTrials,smallRewardTrials),:)];
        SmallBehDataContra       = [SmallBehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(leftStimTrials,smallRewardTrials),:)];
        ErrorBehDataIpsi         = [ErrorBehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(rightStimTrials,errorTrials),:)];
        ErrorBehDataContra       = [ErrorBehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(intersect(leftStimTrials,errorTrials),:)];
    end
    
end
end

        actionTimesLargeIpsi = downSample*(abs(sStart) + LargeBehDataIpsi(:,10)-LargeBehDataIpsi(:,13)); 
        actionTimesLargeContra = downSample*(abs(sStart) + LargeBehDataContra(:,10)-LargeBehDataContra(:,13));
        actionTimesSmallIpsi = downSample*(abs(sStart) + SmallBehDataIpsi(:,10)-SmallBehDataIpsi(:,13));
        actionTimesSmallContra = downSample*(abs(sStart) + SmallBehDataContra(:,10)-SmallBehDataContra(:,13));
        actionTimesErrorIpsi = downSample*(abs(sStart) + ErrorBehDataIpsi(:,10)-ErrorBehDataIpsi(:,13));
        actionTimesErrorContra = downSample*(abs(sStart) + ErrorBehDataContra(:,10)-ErrorBehDataContra(:,13));
        
        outcomeTimesLargeIpsi = downSample*(abs(sStart) + LargeBehDataIpsi(:,14)-LargeBehDataIpsi(:,13)); 
        outcomeTimesLargeContra = downSample*(abs(sStart) + LargeBehDataContra(:,14)-LargeBehDataContra(:,13));
        outcomeTimesSmallIpsi = downSample*(abs(sStart) + SmallBehDataIpsi(:,14)-SmallBehDataIpsi(:,13));
        outcomeTimesSmallContra = downSample*(abs(sStart) + SmallBehDataContra(:,14)-SmallBehDataContra(:,13));
        outcomeTimesErrorIpsi = downSample*(abs(sStart) + ErrorBehDataIpsi(:,14)-ErrorBehDataIpsi(:,13));
        outcomeTimesErrorContra = downSample*(abs(sStart) + ErrorBehDataContra(:,14)-ErrorBehDataContra(:,13));
        
    % ------------------- plot  ----------------------------------
%     if concatenate == 0 || iSession == max(nSessions)
        
        M = [];
        M = {LargeStimDataContra, LargeStimDataIpsi, SmallStimDataContra, SmallStimDataIpsi, ErrorStimDataContra, ErrorStimDataIpsi};
        A = [];
        A = {actionTimesLargeContra, actionTimesLargeIpsi, actionTimesSmallContra, actionTimesSmallIpsi, actionTimesErrorContra, actionTimesErrorIpsi};
        Ar = [];
        Ar = {outcomeTimesLargeContra, outcomeTimesLargeIpsi, outcomeTimesSmallContra, outcomeTimesSmallIpsi, outcomeTimesErrorContra, outcomeTimesErrorIpsi};


%% average signal through trial aligned  at stimulus onset 

        figure; 

subplot(3,1,1) % contra: large reward
data = smooth2a(M{1},0,smooth_factor);
yyaxis left
tempdata = mean(data(:,(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))));
largeCorrectmax = max(tempdata);
tempdata = tempdata./largeCorrectmax;
plot(tempdata, 'LineWidth', 1, 'Color', LargeErrorSmallColor(1,:));
xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) % 
xticklabels([sStart 0 1 2 sStop])
ylim([-0.2 1.1])
ax = gca; 
ax.TickDir = 'out';
hold on; % action time distribution
yyaxis right
histogram(actionTimesLargeContra, 'BinWidth', 90, 'FaceColor', [0 0.58 0.77] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
histogram(outcomeTimesLargeContra, 'BinWidth', 90, 'FaceColor', [0.27 0.74 0], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
ylim([0 400])
hold on; 
yyaxis left
txt = '\nabla';
text(nanmedian(actionTimesLargeContra), tempdata(floor(nanmedian(actionTimesLargeContra)))+0.2, txt)
txt2 = '\nabla';
text(nanmedian(outcomeTimesLargeContra), tempdata(floor(nanmedian(outcomeTimesLargeContra)))+0.2, txt2)


subplot(3, 1, 2) % contra: large reward and error
lines = [1 5];
for i = 1:2
    data = [];
    data = smooth2a(M{lines(i)},0,smooth_factor);
    yyaxis left
    tempdata = mean(data(:,(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))));
    tempdata = tempdata./largeCorrectmax;
    hold on;
    plot(tempdata, 'LineWidth', 1, 'Color', LargeErrorSmallColor(i,:), 'LineStyle', '-');
    
end
xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) %
xticklabels([sStart 0 1 2 sStop])
ylim([-0.2 1.1])
ax = gca;
ax.TickDir = 'out';

subplot(3,1,3) % contra: large reward, small reward, and error
lines = [1 5 3];
for i = 1:3
    data = [];
    data = smooth2a(M{lines(i)},0,smooth_factor);
    yyaxis left
     tempdata = mean(data(:,(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))));
     tempdata = tempdata./largeCorrectmax;
    hold on;
    plot(tempdata, 'LineWidth', 1, 'Color', LargeErrorSmallColor(i,:), 'LineStyle', '-');
    
end
xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) %
xticklabels([sStart 0 1 2 sStop])
ylim([-0.2 1.1])
ax = gca;
ax.TickDir = 'out';

 linkaxes(get(gcf,'children'),'x')
 
 % alt action color = [0.07 0.62 1.00]

