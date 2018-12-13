
function [TrialTimingData] = Salvatore_PrepareNeuronBehData(animal_name,exp_date,exp_series)


% This function is preparing datamatrix from saved data files.

% animals_name      : animal name as stored in the MiceBehInfo database
% exp_id            : as stored in the MiceBehInfo database
% good_sessions     : 1 only good session, 0 only bad sessions, [] all sessions
% bait              : 1 only sessions with bait, 0 only sessions without bait
%                     []all sessions
% single_session    : if not empty then the function is only looking at that
%                     particular session. If empty then it looks at all sessions which satisfy
%                     the requirements set by other parameters.
% ignorefirsttrials : if trial number then in each session trials nefore
%                     the number will be removed. if empty then all trials included


% Note : be carefull about asinging laser to reward time. changed and old
% sessions not checked. 17 March 2015

% Note : block types..only valid for reward time

% block_type = 1;  % laser at rewards on left
% block_type = 2;  % laser at rewards on right
% block_type = 3;  % laser at rewards on both sides
% block_type = 4;  % no laser at reward time

% AL 05/03/2015 London
% AL modified to include hide cue delay Dec 2015
% AL modified to include situation when laser@reward on both sides have
% different number of pulses   22 Dec 2015

%addpath '\\ZSERVER\Code\PsychoKiller';
%addpath '\\ZSERVER\Code\Matteobox\';


% AL 02/08/2015 London
% AL added Action_Onset_Time 2016-01-13
% AL separately added Go_Cue_Time 2016-06-28
% AL added feedbackdeliveryDelay to correctly comute reward time



%load('MiceExpInfo')


path2Beh= ['\\zserver\Data\expInfo\',animal_name,'\',exp_date,'\',exp_series];

path2Neuron= ['\\ZSERVER\Data\multichanspikes',animal_name,'\',exp_date];

addpath (genpath(path2Beh))
addpath (genpath(path2Neuron))



filename_block = [exp_date,'_',num2str(exp_series),'_',animal_name,'_Block','.mat']

filename_parameter = [exp_date,'_',num2str(exp_series),'_',animal_name,'_parameters','.mat']


DataMatrixs=[];



load(filename_block);
load(filename_parameter);


%% pre-allocation
Sign_Contrast    = nan(1,size(block.trial,2)-1);
ResponseSide     = nan(1,size(block.trial,2)-1);
liquid           = nan(1,size(block.trial,2)-1);
laser            = nan(1,size(block.trial,2)-1);
Correction_Trial = nan(1,size(block.trial,2)-1);
TrialCount       = nan(1,size(block.trial,2)-1);
laser_Stimtime   = nan(1,size(block.trial,2)-1);
response_time    = nan(1,size(block.trial,2)-1);
block_type       = nan(1,size(block.trial,2)-1);
session_ID_Num   = nan(1,size(block.trial,2)-1);
correct          = nan(1,size(block.trial,2)-1);
InivsFinalMove   = nan(1,size(block.trial,2)-1);
reaction_time    = nan(1,size(block.trial,2)-1);
trajectory_ratio = nan(1,size(block.trial,2)-1);
Sound_Onset_Time = nan(1,size(block.trial,2)-1);
Stim_Onset_Time  = nan(1,size(block.trial,2)-1);
Reward_Onset_Time = nan(1,size(block.trial,2)-1);
Action_Onset_Time = nan(1,size(block.trial,2)-1);
Go_Cue_Time       = nan(1,size(block.trial,2)-1);
Liquidsize       = nan(1,size(block.trial,2)-1);




for trial = 1 : size(block.trial,2)-1
    
    if block.trial(trial).condition.visCueContrast(1) > 0
        
        Sign_Contrast(trial) = -block.trial(trial).condition.visCueContrast(1);
    elseif block.trial(trial).condition.visCueContrast(2) > 0
        
        Sign_Contrast(trial) = block.trial(trial).condition.visCueContrast(2);
    else
        
        Sign_Contrast(trial) = 0;
        
    end
    
    ResponseSide(trial) = block.trial(trial).responseMadeID;
    
    liquid(trial) = block.trial(trial).feedbackType;
    liquid(liquid==-1) = 0;
    
    correct(trial) = block.trial(trial).feedbackType;
    correct(correct==-1) = 0;
    
    
    Correction_Trial(trial) = block.trial(trial).condition.repeatNum;
    
    Sound_Onset_Time (trial) = block.trial(trial).onsetToneSoundPlayedTime(1);  % go cue time is the second component
    if length(block.trial(trial).onsetToneSoundPlayedTime) > 1  % both initial beep and go cue
        
        Go_Cue_Time (trial) = block.trial(trial).onsetToneSoundPlayedTime(2);
        
    end
    
    Stim_Onset_Time  (trial) = block.trial(trial).stimulusCueStartedTime;
    
    Reward_Onset_Time (trial) = block.trial(trial).feedbackStartedTime + block.parameters.feedbackDeliveryDelay;
    
    
    
    if isfield (block.trial(1).condition,'rewardVolume')   %  we dont have this field in some of the sessions (globalised param)
                                                Liquidsize(trial)=block.trial(trial).condition.rewardVolume(1);
                                       

        
        if length(find(block.trial(trial).condition.rewardVolume))==1
            
            laser(trial) = 0;
            
        elseif length(find(block.trial(trial).condition.rewardVolume)) > 1 && liquid(trial)==1
            
            laser(trial) = block.trial(trial).condition.rewardVolume(2); % 1;
            
        else
            laser(trial) = 0;
        end
        
        %% putting solenoid valve click -small liquid is set to make the click-
        if block.trial(trial).condition.rewardVolume(1) < 0.6
            liquid(trial) = 0;
        end
        
    else
        
                Liquidsize(trial)=parameters.rewardVolume(1);

        if size(parameters.rewardVolume,1)==1
            
            laser(trial) = 0;
            
        elseif size(parameters.rewardVolume,1)>1
            
            %                     if size(parameters.rewardVolume,2)==1 && parameters.rewardVolume(2)==0
            %
            %                         laser(trial) = 0;
            if size(parameters.rewardVolume,2)==1 &&  liquid(trial)==1 % be careful here..might not work for old sessions
                laser(trial) = parameters.rewardVolume(2); % 1;
            elseif size(parameters.rewardVolume,2)==1  && liquid(trial)==0 % be careful here..might not work for old sessions
                laser(trial) = 0;
                
            elseif size(parameters.rewardVolume,2)>1
                
                if parameters.rewardVolume(2) && ResponseSide(trial)==1 && liquid(trial)==1
                    laser(trial) = parameters.rewardVolume(ResponseSide,find(abs(Sign_Contrast)==parameters.visCueContrast(ResponseSide,:)));  %1;
                    
                    %                         elseif parameters.rewardVolume(2)==0 && ResponseSide(trial)==2 && liquid(trial)==1
                    %                             laser(trial) = 1;
                    
                else
                    laser(trial) = 0;
                end
                
            end
            
        end
        if parameters.rewardVolume(1) < 0.6
            liquid(trial) = 0;
        end
        
    end
    
    
    
    %% laser on stimulus time
    if isfield (block.trial(1).condition,'rewardOnStimulus')  %  we dont have this field in some of the sessions (globalised param)
        if any(block.trial(trial).condition.rewardOnStimulus(2))
            laser_Stimtime(trial)=1;
            
        else
            laser_Stimtime(trial) = 0;
            
        end
    elseif isfield (parameters,'rewardOnStimulus')
        
        if size(parameters.rewardOnStimulus,1) > 1 && parameters.rewardOnStimulus(2)>0
            laser_Stimtime(trial)=1;
            
        else
            laser_Stimtime(trial) = 0;
            
        end
    end
    %% response time
    response_time(trial) = block.trial(trial).responseMadeTime - block.trial(trial).interactiveStartedTime;
    
    
    %% reaction time (29.06.2016 changed stimonset_index (for now for animals with go cue.))
    % stimonset_index=find(floor(100*(block.inputSensorPositionTimes))==floor(100*(block.trial(trial).interactiveStartedTime)),1);
    stimonset_index=find(floor(100*(block.inputSensorPositionTimes))==floor(100*(block.trial(trial).stimulusCueStartedTime)),1); 
    
    
    
    if ~isempty(stimonset_index)
        
        actiononset_index= find(abs(diff(block.inputSensorPositions(stimonset_index:end)))>1,1);  % number defines a arbit threshold
        
        
        if  ~isempty(actiononset_index)
            %  if length(block.trial(trial).onsetToneSoundPlayedTime) < 2  % only initial beep but no go time ..therefore we are interested in action onset
            
            Action_Onset_Time (trial) = block.inputSensorPositionTimes(actiononset_index+stimonset_index);
            
            %  end
            
            reaction_time(trial) =block.inputSensorPositionTimes(actiononset_index+stimonset_index) -  block.inputSensorPositionTimes(stimonset_index);
            
            
            if isempty(block.trial(trial).inputThresholdCrossedTime)
                
                threshold_index =[];
            else
                threshold_index=find(floor(100*(block.inputSensorPositionTimes))==floor(100*(block.trial(trial).inputThresholdCrossedTime)),1);
                
            end
            
            
            if ~isempty(threshold_index)
                
                movement_trace=block.inputSensorPositions(stimonset_index:threshold_index)-block.inputSensorPositions(stimonset_index);
                movement_time= block.inputSensorPositionTimes(stimonset_index:threshold_index)-block.inputSensorPositionTimes(stimonset_index);
                
                if ResponseSide(trial) == 1 && nansum(movement_trace(find(movement_trace < 0))) ~= 0  % left choice
                    
                    trajectory_ratio(trial) = abs(nansum(movement_trace(find(movement_trace > 0))) ./ nansum(movement_trace(find(movement_trace < 0))));
                elseif ResponseSide(trial) == 2 && nansum(movement_trace(find(movement_trace > 0))) ~= 0  % right choice
                    
                    trajectory_ratio(trial) = abs(nansum(movement_trace(find(movement_trace < 0))) ./ nansum(movement_trace(find(movement_trace > 0))));
                    
                end
                     
                % comparing initial movement direction versus final movement direction
                InivsFinalMove(trial)=sign((block.inputSensorPositions(actiononset_index+stimonset_index)-block.inputSensorPositions(stimonset_index)) ...
                    * (block.inputSensorPositions(threshold_index)-block.inputSensorPositions(stimonset_index)));
            end
        end

    else
        reaction_time(trial)=nan;
    
end
    
    TrialCount(trial)=trial;
end

response_time= zscore(response_time);

reaction_time = (reaction_time - nanmean(reaction_time)) ./ nanstd(reaction_time);

ResponseSide(ResponseSide==1) = -1;
ResponseSide(ResponseSide==2) = 1 ;

%% define the block type
if sum(laser)==0
    laserTimesResponseSide=laser .* ResponseSide;
else
    laser=laser== max(laser);
    laserTimesResponseSide=laser .* ResponseSide;
end

%laserTimesResponseSide=laser .* ResponseSide;

if sum(ismember(unique(laserTimesResponseSide(isfinite(laserTimesResponseSide))),[-1 0 1]))  ==3
    block_type(1:trial) = 3;  % laser at rewards on both sides
elseif sum(ismember(unique(laserTimesResponseSide(isfinite(laserTimesResponseSide))),[-1 0])) ==2
    block_type(1:trial) = 1;  % laser at rewards on left
elseif sum(ismember(unique(laserTimesResponseSide(isfinite(laserTimesResponseSide))),[0 1]))   ==2
    block_type(1:trial) = 2;  % laser at rewards on right
    
elseif sum(ismember(unique(laserTimesResponseSide(isfinite(laserTimesResponseSide))),[0]))     ==1
    block_type(1:trial) = 4;  % no laser at reward time
    
end

% define blocks for asymmetric reward size
        if length(unique(Liquidsize)) > 1   % 2 different liquid size
            
            AchievedLiquid = Liquidsize .* correct;
            
            if nonzeros(unique(AchievedLiquid(ResponseSide>0))) > nonzeros(unique(AchievedLiquid(ResponseSide<0)))
                block_type(1:trial) = 2;  % large liquid on right
            elseif nonzeros(unique(AchievedLiquid(ResponseSide>0))) < nonzeros(unique(AchievedLiquid(ResponseSide<0)))
                block_type(1:trial) = 1;  % large liquid on left
            end
            
        end
        
if strcmp(animal_name,'ALK011') && ~ ismember(filename_block,{...
        
    '2016-01-08_3_ALK011_Block.mat';...
    '2016-01-11_3_ALK011_Block.mat';...
    '2016-01-13_3_ALK011_Block.mat';...
    '2016-01-13_4_ALK011_Block.mat';...
    '2016-01-14_3_ALK011_Block.mat';...
    '2016-01-14_4_ALK011_Block.mat';...
    '2016-01-14_6_ALK011_Block.mat';...
    '2016-01-14_7_ALK011_Block.mat';...
    '2016-01-16_3_ALK011_Block.mat';...
    '2016-01-16_4_ALK011_Block.mat'})

liquid(1:trial) = zeros(1,trial);

end

%   Sign_Contrast(Sign_Contrast==0) = 0.01;

%   correct = sign(Sign_Contrast .* ResponseSide);
%   correct =
%   correct(correct==-1)=0;

%  Sign_Contrast(Sign_Contrast==0.01) = 0;


%DataMatrix=[TrialCount',Sign_Contrast',ResponseSide',liquid', laser',laser_Stimtime', response_time', block_type', correct',trajectory_ratio',...
%    reaction_time', Sound_Onset_Time',Stim_Onset_Time',Reward_Onset_Time' Correction_Trial'];

DataMatrix=[TrialCount',Sign_Contrast',ResponseSide',liquid', laser',laser_Stimtime', response_time', block_type', correct',Action_Onset_Time',...
    reaction_time', Sound_Onset_Time',Stim_Onset_Time',Reward_Onset_Time' Correction_Trial' Go_Cue_Time']; % Go_Cue_Time will eventually stay 15. see Salvatore_PrepareTrialTimingData



%%
%  session = session + 1;

%   DataMatrix_AllSessions = [DataMatrix; DataMatrix_AllSessions];


if ~strcmp(animal_name,'ALK068') && ~strcmp(animal_name,'ALK070') && ~strcmp(animal_name,'ALK083')...
    && ~strcmp(animal_name,'ALK084')
    
DataMatrix(DataMatrix(:,3)==3,:) = [];
end


TrialTimingData = DataMatrix;



