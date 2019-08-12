
close all
clear all



% ----- enter reqs --------------------------------------------------------
animal_ID = 70
brain_region = 'DMS'
exp_ID = '7'


stim_2_plot = 0.25; %should be positive
smooth_factor = 150;
colorRange= getColorRange(animal_ID, exp_ID); % color range for plotting imagesc data

RT_min = 0.2;
RT_max = 2; % reaction time range to include

sStart = -0.3; %stimulus
sStop = 3;
downSample = 1200; % sampling rate AFTER DOWNSAMPLING in makeDataSummaryMatrix
eventOnset = 3700;  % in the saved matrix, the event occurs in this column

load(['BehPhotoM_Exp',exp_ID,'_',brain_region,'.mat'])

fullStim = [-stim_2_plot stim_2_plot];

errorData = [];
largeRewardData = [];
smallRewardData = [];

errorOutcomeTimes = [];
smallRewardOutcomeTimes = [];
largeRewardOutcomeTimes = [];

errorActionTimes = [];
smallRewardActionTimes = [];
largeRewardActionTimes = [];

for iSession = 1:length(BehPhotoM(animal_ID).Session)
    
    RTs = [];
    OutcomeTs = [];
    TempBehData = [];
    iErrors = [];
    iLargeReward = [];
    iSmallReward = [];
    
    RTs = BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,10) - BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,13);
    OutcomeTs = BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,14) - BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,13);
    RTExcludeTrials = [(find(RTs < RT_min)) ; (find(RTs >RT_max))];
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    iErrors = setdiff(intersect(find(TempBehData(:,9)==0), find(abs(TempBehData(:,2))==stim_2_plot)), RTExcludeTrials); 
    iLargeReward = [setdiff(mintersect(find(TempBehData(:,8)==1), find(TempBehData(:,9)==1), find(TempBehData(:,3)==-1)), RTExcludeTrials);...
        setdiff(mintersect(find(TempBehData(:,8)==2), find(TempBehData(:,9)==1), find(TempBehData(:,3)==1), find(abs(TempBehData(:,2))==stim_2_plot)), RTExcludeTrials)];
    iSmallReward = [setdiff(mintersect(find(TempBehData(:,8)==1), find(TempBehData(:,9)==1), find(TempBehData(:,3)==1)), RTExcludeTrials); ...
        setdiff(mintersect(find(TempBehData(:,8)==2), find(TempBehData(:,9)==1), find(TempBehData(:,3)==-1), find(abs(TempBehData(:,2))==stim_2_plot)), RTExcludeTrials)];
    
    TempStimData = [];
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimL')
        TempStimData = BehPhotoM(animal_ID).Session(iSession).NeuronStimL;
        errorData = [errorData; TempStimData(iErrors,:)];
        smallRewardData = [smallRewardData; TempStimData(iSmallReward,:)];
        largeRewardData = [largeRewardData; TempStimData(iLargeReward,:)];
        
        errorOutcomeTimes = [errorOutcomeTimes; OutcomeTs(iErrors,:)];
        smallRewardOutcomeTimes = [smallRewardOutcomeTimes; OutcomeTs(iSmallReward,:)];
        largeRewardOutcomeTimes = [largeRewardOutcomeTimes; OutcomeTs(iLargeReward,:)];
        
        errorActionTimes = [errorActionTimes; RTs(iErrors)];
        smallRewardActionTimes = [smallRewardActionTimes; RTs(iSmallReward)];
        largeRewardActionTimes = [largeRewardActionTimes; RTs(iLargeReward)];        
    end
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimR')
        TempStimData = BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
        errorData = [errorData; TempStimData(iErrors,:)];
        smallRewardData = [smallRewardData; TempStimData(iSmallReward,:)];
        largeRewardData = [largeRewardData; TempStimData(iLargeReward,:)];
        
        errorOutcomeTimes = [errorOutcomeTimes; OutcomeTs(iErrors,:)];
        smallRewardOutcomeTimes = [smallRewardOutcomeTimes; OutcomeTs(iSmallReward,:)];
        largeRewardOutcomeTimes = [largeRewardOutcomeTimes; OutcomeTs(iLargeReward,:)];

        errorActionTimes = [errorActionTimes; RTs(iErrors)];
        smallRewardActionTimes = [smallRewardActionTimes; RTs(iSmallReward)];
        largeRewardActionTimes = [largeRewardActionTimes; RTs(iLargeReward)];
        
    end        
        
end

%sorting by RT
[errorActionTimes, errorActionSort] = sort(errorActionTimes);
[smallRewardActionTimes, smallRewardActionSort] = sort(smallRewardActionTimes);
[largeRewardActionTimes, largeRewardActionSort] = sort(largeRewardActionTimes);

errorOutcomeTimes = errorOutcomeTimes(errorActionSort);
smallRewardOutcomeTimes = smallRewardOutcomeTimes(smallRewardActionSort);
largeRewardOutcomeTimes = largeRewardOutcomeTimes(largeRewardActionSort);

plotIndex = [errorActionSort; smallRewardActionSort+length(errorActionSort); ...
    largeRewardActionSort+length(errorActionSort)+length(smallRewardActionSort)];

figure; 
%%
masterData = [errorData; smallRewardData; largeRewardData];
actionTimes = [errorActionTimes; smallRewardActionTimes; largeRewardActionTimes]*downSample + abs(sStart*downSample);
outcomeTimes = [errorOutcomeTimes; smallRewardOutcomeTimes; largeRewardOutcomeTimes]*downSample + abs(sStart*downSample);

imagesc(smooth2a(masterData(plotIndex, (eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))), 0, smooth_factor), colorRange)
colormap('bluewhitered')
hold on
for ievent = 1:length(actionTimes)
    if ievent <= size(errorData,1)
        plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0.7 0 0.2 ], 'MarkerFaceColor', [0.7 0 0.2 ], 'MarkerSize', 3)
    elseif size(errorData,1) < ievent && ievent <= size(errorData,1)+size(smallRewardData,1)
        plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0 1 0 ], 'MarkerFaceColor', [0 1 0 ], 'MarkerSize', 3)
    elseif ievent > size(errorData,1)+size(smallRewardData,1)
        plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 3)
    end
end

for ievent = 1:length(outcomeTimes)
    if ievent <= size(errorData,1)
        plot(outcomeTimes(ievent),ievent, 'x', 'MarkerEdgeColor', [0.7 0 0.2 ], 'MarkerFaceColor', [0.7 0 0.2 ], 'MarkerSize', 5)
    elseif size(errorData,1) < ievent && ievent <= size(errorData,1)+size(smallRewardData,1)
        plot(outcomeTimes(ievent),ievent, 'x', 'MarkerEdgeColor', [0 1 0 ], 'MarkerFaceColor', [0 1 0 ], 'MarkerSize', 5)
    elseif ievent > size(errorData,1)+size(smallRewardData,1)
        plot(outcomeTimes(ievent),ievent, 'x', 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 5)
    end
end
xticks([1 abs(sStart*downSample) (sStop-sStart)*downSample-1])
xticklabels([sStart 0 sStop])
ylabel('Trials')
xlabel('Time from stimulus (s)')
line([abs(sStart*downSample) abs(sStart*downSample)], [0 size(masterData,1)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 1.5); % stim onset line




%%

function [colorRange] = getColorRange(animal_ID, expID)

if animal_ID == 48
    if strcmp(expID, '23')
        [colorRange] = [-3.5 5]
    end
elseif animal_ID == 51
    if strcmp(expID, '23')
        [colorRange] = [-10 25]
    end    
elseif animal_ID == 72
    if strcmp(expID, '23')
        [colorRange] = [-3 4]
    end
elseif animal_ID == 71
    if strcmp(expID, '23')
        [colorRange] = [-3 8]
    end
elseif animal_ID == 64
    if strcmp(expID, '23')
        [colorRange] = [-8 20]
    end
elseif animal_ID == 63
    if strcmp(expID, '23')
        [colorRange] = [-3 3]
    end
elseif animal_ID == 57
    if strcmp(expID, '23') || strcmp(expID, '7')
        [colorRange] = [-2 2]
    end
else [colorRange] = [-3 6]
    
end

end