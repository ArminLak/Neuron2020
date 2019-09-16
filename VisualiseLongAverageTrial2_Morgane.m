% creates 2 rasters for contra and ipsi (per animal or per session) ordered by
% Error/Small/Large and by RTs. Works for unihem and bihem animals
% Morgane April 2019

% can you see this? yes

% for information about normalisation, check end of this code. 

% close all

clear all


% ----- enter reqs --------------------------------------------------------

%DMS 
%  if exist('brain_region', 'var') && strcmp(brain_region, 'DMS')
%      clearvars -except BehPhotoM
%  else clear all
%  end
%  Animals = [53, 62, 63, 71,72]
%  brain_region = 'DMS'

% NAC
%  if exist('brain_region', 'var') && strcmp(brain_region, 'NAc')
%      clearvars -except BehPhotoM
%   else clear all
%  end
% Animals = [56 57 59 66]
% brain_region = 'NAc'

% VTA: 
 if exist('brain_region', 'var') && strcmp(brain_region, 'VTA')
     clearvars -except BehPhotoM
 else clear all
 end
Animals = [48 50 51 64]
brain_region = 'VTA'




exp_ID = '23'
plotRasters = 0;

concatenate = 1; %show all trials across all sessions for an animal
IpsiContra = 1; %separate rasters based on ipsi / contra? 

stim_2_plot = 0.5; %should be positive
smooth_factor = 100;

[ActionColor, IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors();
% colorRange= getColorRange(animal_ID, exp_ID); % color range for plotting imagesc data

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


load(['BehPhotoM_Exp',exp_ID,'_',brain_region,'.mat'])

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
    
    RTExcludeTrials = [(find(RTs < RT_min)) ; (find(RTs >RT_max)) ; (find(OTs < OT_min)) ; (find(OTs > OT_max))];

    leftStimTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2)== -(stim_2_plot)),RTExcludeTrials);
    rightStimTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2)==stim_2_plot), RTExcludeTrials);
   
    correctTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1), RTExcludeTrials);
    largeRewardTrials = setdiff(unique([find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,8)==1); ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1); ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,3)==-1)]), RTExcludeTrials);
    smallRewardTrials = setdiff(unique([find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,8)==2); ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,9)==1); ...
        find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,3)==1)]), RTExcludeTrials);
    
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

        actionTimesLargeIpsi = downsample*(abs(sStart) + LargeBehDataIpsi(:,10)-LargeBehDataIpsi(:,13)); 
        actionTimesLargeContra = downsample*(abs(sStart) + LargeBehDataContra(:,10)-LargeBehDataContra(:,13));
        actionTimesSmallIpsi = downsample*(abs(sStart) +SmallBehDataIpsi(:,10)-SmallBehDataIpsi(:,13));
        actionTimesSmallContra = downsample*(abs(sStart) + SmallBehDataContra(:,10)-SmallBehDataContra(:,13));
        actionTimesErrorIpsi = downsample*(abs(sStart) + ErrorBehDataIpsi(:,10)-ErrorBehDataIpsi(:,13));
        actionTimesErrorContra = downsample*(abs(sStart) + ErrorBehDataContra(:,10)-ErrorBehDataContra(:,13));
        
        outcomeTimesLargeIpsi = downsample*(abs(sStart) + LargeBehDataIpsi(:,14)-LargeBehDataIpsi(:,13)); 
        outcomeTimesLargeContra = downsample*(abs(sStart) + LargeBehDataContra(:,14)-LargeBehDataContra(:,13));
        outcomeTimesSmallIpsi = downsample*(abs(sStart) +SmallBehDataIpsi(:,14)-SmallBehDataIpsi(:,13));
        outcomeTimesSmallContra = downsample*(abs(sStart) + SmallBehDataContra(:,14)-SmallBehDataContra(:,13));
        outcomeTimesErrorIpsi = downsample*(abs(sStart) + ErrorBehDataIpsi(:,14)-ErrorBehDataIpsi(:,13));
        outcomeTimesErrorContra = downsample*(abs(sStart) + ErrorBehDataContra(:,14)-ErrorBehDataContra(:,13));
        
    % ------------------- plot  ----------------------------------
%     if concatenate == 0 || iSession == max(nSessions)
        
        M = [];
        M = {LargeStimDataContra, LargeStimDataIpsi, SmallStimDataContra, SmallStimDataIpsi, ErrorStimDataContra, ErrorStimDataIpsi};
        A = [];
        A = {actionTimesLargeContra, actionTimesLargeIpsi, actionTimesSmallContra, actionTimesSmallIpsi, actionTimesErrorContra, actionTimesErrorIpsi};
        Ar = [];
        Ar = {outcomeTimesLargeContra, outcomeTimesLargeIpsi, outcomeTimesSmallContra, outcomeTimesSmallIpsi, outcomeTimesErrorContra, outcomeTimesErrorIpsi};


%% average signal through trial aligned  at stimulus onset 

        actionTimesContra = [];
        actionTimesContra = abs(sStart)*downSample + (largeRewTrialsContra(:,1))*downSample;
        outcomeTimesContra = [];
        outcomeTimesContra = abs(sStart)*downSample + (largeRewTrialsContraOutcome(:,1))*downSample;
        
        figure; 

subplot(3,1,1) % contra: large reward
data = smooth2a(M{1},0,smooth_factor);
yyaxis left
tempdata = mean(data(length(errorTrialsContra)+length(smallRewTrialsContra)+1:end,(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))));
tempdata = tempdata./max(tempdata);
plot(tempdata, 'LineWidth', 1, 'Color', 'k');
xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) % 
xticklabels([sStart 0 1 2 sStop])
ylim([-0.2 1.1])
ax = gca; 
ax.TickDir = 'out';


hold on; % action time distribution
yyaxis right
histogram(actionTimesContra, 'BinWidth', 90, 'FaceColor', [0 0.58 0.77] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
txt = '\nabla';
text(nanmedian(actionTimesContra), tempdata(floor(nanmedian(actionTimesContra)))+0.2, txt)
ylim([0 400])

hold on; % reward time distribution
yyaxis right 
histogram(outcomeTimesContra, 'BinWidth', 90, 'FaceColor', [0.27 0.74 0], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
txt2 = '\nabla';
text(nanmedian(outcomeTimesContra), tempdata(floor(nanmedian(outcomeTimesContra)))+0.2, txt2)

subplot(3, 1, 2) % contra: large reward and error
for i = [1 5]
data = smooth2a(M{1},0,smooth_factor);
yyaxis left
tempdata = mean(data(length(errorTrialsContra)+length(smallRewTrialsContra)+1:end,(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))));
tempdata = tempdata./max(tempdata);
plot(tempdata, 'LineWidth', 1, 'Color', 'k');
xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) % 
xticklabels([sStart 0 1 2 sStop])
ylim([-0.2 1.1])
ax = gca; 
ax.TickDir = 'out';

subplot(3,1,2) % contra: large reward, small reward, and error


 linkaxes(get(gcf,'children'),'x')
 
 % alt action color = [0.07 0.62 1.00]
            
%% script-specficic functions

function [colorRange] = getColorRange(animal_ID, expID)

if animal_ID == 48
    if strcmp(expID, '23')
        [colorRange] = [-3 4]
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
    if strcmp(expID, '23')
        [colorRange] = [-2 2]
    end
else [colorRange] = [-3 6]
    
end

end

function [ActionColor, IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors()
IpsiContraColor = [ 0 0 0
                    153/255 51/255 1
                    ];
                
ErrorCorrectColor = [   'r'
                        'g'
                        ];
                    
SmallLargeColor = [ 0.15 0.75 0
                    0.15 0.4 0
                    ];
                
ActionColor = [0.35 0.65 0.86];

end
