% this code does all the t-tests and linear regressions for morgane's masters dissertation 
% inspired by 'VisualiseMultipleSessions'

% to use this code run section 1 and then run whichever section you want according to which tests you want to carry out


%% section 1: load and organise data 

clear all
close all

%VTA : [48, 50,51]  coresponding to ALK068, 70 and 71
% NAc : [57] coresponding to MMM001
% DMS : [53] coresponding to ALK074


% select animal
animal_ID = 50

% select database
load('BehPhotoM_Exp23')
%load('BehPhotoM_Exp23_NAc')
%load('BehPhotoM_Exp23_DMS')

sample_rate = 12000;                                        % photoM recording sampling rate
downSample = 1200;                                       % factor downsampling the Ca responses
eventOnset = 3700;  % in the saved matrix, the ev

% define implant
Implant = 'Un' 

RTLimit = 10; % in s, excluding trials with RT longer than this

BehData = [];
BeepData = [];
StimData = [];
ActionData = [];
RewardData = [];

sessionz = 1:length(BehPhotoM(animal_ID).Session);



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

%exclude trials where RT is too long (defined above)
% 
% RT = BehData(:,10) - BehData(:,13);
% 
% toRemove = find ( RT > RTLimit);
% BehData(toRemove,:) = [];


correct = find(TempBehData(:,9)==1);
error = find(TempBehData(:,9)==0);

%% section 2.1 : non parametric for each animal pre- and post- event


%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0;
postStop = 0.05;


%if looking only at e.g. correct trials of 0.5 contrast use this as 1st index of line 86:
tempIndex = intersect(find(abs(TempBehData(:,2))==0.5),correct);


% pre-post Stimulus significance -----------------------------------------
preStimAverages = nanmean(StimData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postStimAverages = nanmean(StimData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostStimpvalue, prepostStimNullReject] = ranksum(preStimAverages, postStimAverages)


%pre-post Outcome significance (only correct trials) ---------------------
preOutcomeAverages = nanmean(RewardData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postOutcomeAverages = nanmean(RewardData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostOutcomepvalue, prepostOutcometNullReject] = ranksum(preOutcomeAverages, postOutcomeAverages)


%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.2;
postStop = 0.8;

%pre-post Reward significance (only correct trials) ---------------------
preCorrectAverages = nanmean(RewardData(correct, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postCorrectAverages = nanmean(RewardData(correct, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostCorrectpvalue, prepostCorrectNullReject] = ranksum(preCorrectAverages, postCorrectAverages)


%pre-post Error significance (only error trials) ------------------------
preErrorAverages = nanmean(RewardData(error, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postErrorAverages = nanmean(RewardData(error, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostErrorpvalue, prepostErrorNullReject] = ranksum(preErrorAverages, postErrorAverages)


%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0;
postStop = 0.05;

%pre-post Action significance ------------------------
preActionAverages = nanmean(ActionData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postActionAverages = nanmean(ActionData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostActionpvalue, prepostActionNullReject] = ranksum(preActionAverages, postActionAverages)


%% section 2.2 : non parametric test for response large vs small reward same contrast 

%pre and post event time windows (seconds):
Start = 0.2;
Stop = 0.8;

largeRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
smallRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);
    
% large vs small side Stimulus significance -----------------------------------------
largeRewStimAverages = nanmean(StimData(largeRew, (eventOnset+(Start*downSample)):(eventOnset+(Stop*downSample))), 2);
smallRewStimAverages = nanmean(StimData(smallRew, (eventOnset+(Start*downSample)):(eventOnset+(Stop*downSample))), 2);

[largesmallStimpvalue, largesmallStimNullReject] = ranksum(largeRewStimAverages, smallRewStimAverages)



%% section 3 : linear regressions for contrast sensitivity of response


% create list of responses to stim 
% create list of contrasts in each response to stim 

% ---- stim contrast v stim response regression analysis -----------
%pre and post event time windows (seconds):
start = 0.2;
stop = 0.4;

stimResponses = nanmean(StimData(:, eventOnset + (start*downSample): eventOnset + (stop*downSample)), 2);
stimContrasts = abs(BehData(:,2));

stimRegStats = regstats( stimContrasts, stimResponses);
StimContrastFstat = stimRegStats.fstat.f
StimContrastFstatPVal = stimRegStats.fstat.pval


% ---- stim contrast v stim response regression analysis -----------
%pre and post event time windows (seconds):
start = 0.2;
stop = 0.4;

stimResponses = nanmean(StimData(:, eventOnset + (start*downSample): eventOnset + (stop*downSample)), 2);
stimContrasts = abs(BehData(:,2));

stimRegStats = regstats( stimContrasts, stimResponses);
StimContrastFstat = stimRegStats.fstat.f
StimContrastFstatPVal = stimRegStats.fstat.pval




%% section 4: ANOVA for response to large rwd, small rwd, or no rwd










