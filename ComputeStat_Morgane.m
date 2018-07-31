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
animal_ID = 51

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

%% section 2 : mann whitney U for each animal pre- and post- event


%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.2;
postStop = 0.4;


%if looking only at e.g. correct trials of 0.5 contrast use this as 1st index of line 86:
tempIndex = intersect(find(abs(TempBehData(:,2))==0.5),correct);


% pre-post Stimulus significance -----------------------------------------
preStimAverages = nanmean(StimData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postStimAverages = nanmean(StimData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostStimpvalue, prepostStimNullReject] = ranksum(preStimAverages, postStimAverages)


%pre-post Reward significance (only correct trials) ---------------------
preCorrectAverages = nanmean(RewardData(correct, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postCorrectAverages = nanmean(RewardData(correct, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostCorrectpvalue, prepostCorrectNullReject] = ranksum(preCorrectAverages, postCorrectAverages)


%pre-post Reward significance (only error trials) ------------------------
preErrorAverages = nanmean(RewardData(error, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postErrorAverages = nanmean(RewardData(error, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostErrorpvalue, prepostErrorNullReject] = ranksum(preErrorAverages, postErrorAverages)




%% section 3 : linear regressions for contrast sensitivity of response





