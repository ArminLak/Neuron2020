% creates 2 rasters for contra and ipsi (per animal or per session) ordered by
% Error/Small/Large and by RTs. Works for unihem and bihem animals
% Morgane April 2019



close all
clear all

% DMS animals: [ALK074 is 53; ALK075 is 55; ALK083 is 63; MMM003 is 62; MMM008 is 70; MMM009/10 are 71/72]
% NAc animals: [56, 57,59] coresponding to  ALK078(Bi), MMM001(Un), MMM002(Un)

% ----- enter reqs --------------------------------------------------------
animal_ID = 57
brain_region = 'NAC'
exp_ID = '23'

concatenate = 1; %show all trials across all sessions for an animal

stim_2_plot = 0.25; %should be positive
smooth_factor = 100;
colorRange=[-1 3]; % color range for plotting imagesc data

% ------------------ start stop times for task events in second ------------------

sStart = -0.3; %stimulus
sStop = 2;

aStart = -0.6; %action
aStop = 0.2;

rStart = -0.2; %reward
rStop = 0.6;

downSample = 1200; % sampling rate AFTER DOWNSAMPLING in makeDataSummaryMatrix
eventOnset = 3700;  % in the saved matrix, the event occurs in this column

% -------------------------------------------------------------------------


load(['BehPhotoM_Exp',exp_ID,'_',brain_region,'.mat'])

fullStim = [-stim_2_plot stim_2_plot];

nSessions = 1:length(BehPhotoM(animal_ID).Session);

for iSession = nSessions
    
    if iSession ==1 ||concatenate == 0
        StimDataR       = [];
        ActionDataR     = [];
        RewardDataR     = [];
        
        StimDataL       = [];
        ActionDataL     = [];
        RewardDataL     = [];
        
        leftStimTrials = [];
        rightStimTrials = [];
        errorTrials = [];
        smallRewTrials = [];
        
        BehData = [];
        
        StimDataContra      = [];
        ActionDataContra    = [];
        RewardDataContra    = [];
        
        StimDataIpsi      = [];
        ActionDataIpsi    = [];
        RewardDataIpsi    = [];
        
    end
    
    
    BehData = [BehData; BehPhotoM(animal_ID).Session(iSession).TrialTimingData(ismember(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2), fullStim),:)];
    
    % trial type indexes
    leftStimTrials  = [leftStimTrials; find(ismember(BehData(:,2), (-(stim_2_plot))))];
    rightStimTrials = [rightStimTrials; find(ismember(BehData(:,2), (stim_2_plot)))];
    
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimR') %load hemispheric deltaFoverF
        StimDataR        = [StimDataR;      BehPhotoM(animal_ID).Session(iSession).NeuronStimR(ismember(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2), fullStim),:)];
%         ActionDataR      = [ActionDataR;    BehPhotoM(animal_ID).Session(iSession).NeuronActionR(ismember(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2), fullStim),:)];
%         RewardDataR      = [RewardDataR;    BehPhotoM(animal_ID).Session(iSession).NeuronRewardR(ismember(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2), fullStim),:)];
    end
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimL')
        StimDataL        = [StimDataL;      BehPhotoM(animal_ID).Session(iSession).NeuronStimL(ismember(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2), fullStim),:)];
%         ActionDataL      = [ActionDataL;    BehPhotoM(animal_ID).Session(iSession).NeuronActionL(ismember(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2), fullStim),:)];
%         RewardDataL      = [RewardDataL;    BehPhotoM(animal_ID).Session(iSession).NeuronRewardL(ismember(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2), fullStim),:)];
    end
    
    if concatenate == 0 || iSession == max(nSessions) %sort once for each session if not concatenating or at the very end if concatenating
        
        errorTrials(:,2) = [errorTrials; find(BehData(:,9)==0)];
        smallRewTrials(:,2) = [intersect(find(BehData(:,9)==1), intersect(find(BehData(:,8)==1), find(BehData(:,2)<0))); ...
            intersect(find(BehData(:,9)==1), intersect(find(BehData(:,8)==2), find(BehData(:,2)>0)))]; % WIP indexing; to be ordered by RT
        largeRewTrials(:,2) = [intersect(find(BehData(:,9)==1), intersect(find(BehData(:,8)==2), find(BehData(:,2)<0))); ...
            intersect(find(BehData(:,9)==1), intersect(find(BehData(:,8)==1), find(BehData(:,2)>0)))]; % to be ordered by RT
        
        %get RTs for different trial types
        errorTrials(:,1) = BehData(errorTrials(:,2), 10) - BehData(errorTrials(:,2), 13);
        smallRewTrials(:,1) = BehData(smallRewTrials(:,2), 10) - BehData(smallRewTrials(:,2), 13);
        largeRewTrials(:,1) = BehData(largeRewTrials(:,2), 10) - BehData(largeRewTrials(:,2), 13);
        
        %now sort by RT
        errorTrials     = sortrows(errorTrials);
        smallRewTrials  = sortrows(smallRewTrials);
        largeRewTrials  = sortrows(largeRewTrials);
        
        %plotting order
        plotIndex = [errorTrials(:,2); smallRewTrials(:,2); largeRewTrials(:,2)];
        RTs = [errorTrials(:,1); smallRewTrials(:,1); largeRewTrials(:,1)];
        actionTimes = abs(sStart*downSample) + (RTs*downSample);
    end
    
    % ------------------- contra matrices -------------------------------
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimR')
        StimDataContra      = [StimDataContra; StimDataR(leftStimTrials,:)];
%         ActionDataContra    = [ActionDataContra; ActionDataR(leftStimTrials,:)];
%         RewardDataContra    = [RewardDataContra; RewardDataR(leftStimTrials,:)];
    end
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimL')
        StimDataContra      = [StimDataContra; StimDataL(rightStimTrials,:)];
%         ActionDataContra    = [ActionDataContra; ActionDataL(rightStimTrials,:)];
%         RewardDataContra    = [RewardDataContra; RewardDataL(rightStimTrials,:)];
    end
    
    % ------------ ipsi matrices ---------------------------------
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimR')
        StimDataIpsi      = [StimDataIpsi; StimDataR(rightStimTrials,:)];
%         ActionDataIpsi    = [ActionDataIpsi; ActionDataR(rightStimTrials,:)];
%         RewardDataIpsi    = [RewardDataIpsi; RewardDataR(rightStimTrials,:)];
    end
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimL')
        StimDataIpsi      = [StimDataIpsi; StimDataL(leftStimTrials,:)];
%         ActionDataIpsi    = [ActionDataIpsi; ActionDataL(leftStimTrials,:)];
%         RewardDataIpsi    = [RewardDataIpsi; RewardDataL(leftStimTrials,:)];
    end
    
    %%
    % ------------------- plot  ----------------------------------
    if concatenate == 0 || iSession == max(nSessions)
        
        M = [];
        M = {StimDataContra, StimDataIpsi};
        
        figure;
        
        for iSubplot = 1:2
            data = M{iSubplot}; % looping through StimDataContra / Ipsi
            
            subplot(1,2,iSubplot)
            imagesc(smooth2a(data(plotIndex, (eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))), 0, smooth_factor), colorRange)
            colormap('bluewhitered')
            hold on
            for ievent = 1:length(actionTimes)
                if ievent <= length(errorTrials)
                    plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0.7 0 0.2 ], 'MarkerFaceColor', [0.7 0 0.2 ], 'MarkerSize', 3)
                elseif length(errorTrials) < ievent && ievent <= length(errorTrials)+length(smallRewTrials)
                    plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0 1 0 ], 'MarkerFaceColor', [0 1 0 ], 'MarkerSize', 3)
                elseif ievent > length(errorTrials)+length(smallRewTrials)
                    plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 3)
                end
            end
            
            xticks([1 abs(sStart*downSample) (sStop-sStart)*downSample-1])
            xticklabels([sStart 0 sStop])
            ylabel('Trials')
            xlabel('Time from stimulus (s)')
            line([abs(sStart*downSample) abs(sStart*downSample)], [0 length(plotIndex)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 1.5); % stim onset line
            
        end
    end
end

