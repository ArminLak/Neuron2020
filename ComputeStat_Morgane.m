% this code does all the t-tests and linear regressions for morgane's masters dissertation 
% inspired by 'VisualiseMultipleSessions'

% Morgane Moss Jul 2018

% to use this code run section 1 and then run whichever section you want according to which tests you want to carry out
% unless running section 6 in which case no need to run section 1 first :) 

% CONTENTS:
% Section 2: Non-param pre- vs post-event comparisons (including comparing large vs small reward etc)
% Section 3: Linear regressions for stim response sensitivity to contrast
% Section 4: ANOVA for large/small/no reward comparisons for stim and outcome
% Section 5: Non-param point by point comparisons for summary figs (by contrast)
% Section 6: Comparison of linear regression slopes for correct/error and large/small reward for stim and outcome


%% section 1: load and organise data 

clear all
close all

%VTA : [48, 50,51]  coresponding to ALK068, 70 and 71
% NAc : [57] coresponding to MMM001
% DMS : [53] coresponding to ALK074


% select animal
animal_ID = 48

% select database
load('BehPhotoM_Exp23')
%load('BehPhotoM_Exp23_NAc')
%load('BehPhotoM_Exp23_DMS')

sample_rate = 12000;                                        % photoM recording sampling rate
downSample = 1200;                                       % factor downsampling the Ca responses
eventOnset = 3700;  % in the saved matrix, the ev

% define implant
Implant = 'Un' 

RTLimit = 6; % in s, excluding trials with RT longer than this

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

% exclude trials where RT is too long (defined above)

RT = BehData(:,10) - BehData(:,13);

toRemove = find ( RT > RTLimit);
BehData(toRemove,:) = [];
StimData(toRemove,:) = [];
ActionData(toRemove,:) = [];
RewardData(toRemove,:) = [];
BeepData(toRemove,:) = [];



%% section 2.1 : pre vs post STIMULUS non parametric for one animal


%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.2;
postStop = 0.8;


%if looking only at e.g. correct trials of 0.5 contrast use this as 1st index of line 86:
tempIndex = intersect(find(abs(TempBehData(:,2))==0.5),correct);


% pre-post Stimulus significance -----------------------------------------
preStimAverages = nanmean(StimData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postStimAverages = nanmean(StimData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostStimpvalue, prepostStimNullReject] = ranksum(preStimAverages, postStimAverages)


%% section 2.2 : pre vs post OUTCOME non parametric for one animal

%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.2;
postStop = 0.4;

correct = find(BehData(:,9)==1);
error = find(BehData(:,9)==0);

%pre-post Outcome significance (allt trials) ---------------------
preOutcomeAverages = nanmean(RewardData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postOutcomeAverages = nanmean(RewardData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostOutcomepvalue, prepostOutcometNullReject] = signrank(preOutcomeAverages, postOutcomeAverages)


%pre-post Reward significance (only correct trials) ---------------------
preCorrectAverages = nanmean(RewardData(correct, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postCorrectAverages = nanmean(RewardData(correct, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostCorrectpvalue, prepostCorrectNullReject] = signrank(preCorrectAverages, postCorrectAverages)


%pre-post Error significance (only error trials) ------------------------
preErrorAverages = nanmean(RewardData(error, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postErrorAverages = nanmean(RewardData(error, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostErrorpvalue, prepostErrorNullReject] = signrank(preErrorAverages, postErrorAverages)

%% section 2.3 : pre vs post ACTION non parametric for one animal

% pre and post event time windows (seconds):
preStart = -0.5;
preStop = 0;
postStart = 0;
postStop = 0.05;

% pre-post Action significance ------------------------
preActionAverages = nanmean(ActionData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postActionAverages = nanmean(ActionData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

[prepostActionpvalue, prepostActionNullReject] = signrank(preActionAverages, postActionAverages)


%% section 2.4 : LARGE VS SMALL REWARD non parametric for one animal

%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.4;
postStop = 0.8;

largeRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
smallRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);

largeRewpreStimAverages = nanmean(StimData(largeRew, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
largeRewpostStimAverages = nanmean(StimData(largeRew, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
smallRewpreStimAverages = nanmean(StimData(smallRew, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
smallRewpostStimAverages = nanmean(StimData(smallRew, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
    
% large vs small side Stimulus significance -----------------------------------------
largeRewStimAverages = largeRewpostStimAverages - largeRewpreStimAverages;
smallRewStimAverages = smallRewpostStimAverages - smallRewpreStimAverages;

[largesmallStimpvalue, largesmallStimNullReject] = ranksum(largeRewStimAverages, smallRewStimAverages)

%% section 2.5: REWARD VS ERROR non parametric for one animals (stimulus)

%pre and post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.4;
postStop = 0.8;

correct = find(BehData(:,9)==1);
error = find(BehData(:,9)==0);


correctPreStimAverages = nanmean(StimData(correct, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
correctPostStimAverages = nanmean(StimData(correct, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
errorPreStimAverages = nanmean(StimData(error, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
errorPostStimAverages = nanmean(StimData(error, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);


% correct vs error Stimulus response significance -----------------------------------------
correctStimAverages = correctPostStimAverages - correctPreStimAverages;
errorStimAverages = errorPostStimAverages - errorPreStimAverages;

[correrrStimpvalue, correrrStimNullReject] = ranksum(correctStimAverages, errorStimAverages)


%% section 3.1 : linear regressions for contrast sensitivity of response to stimulus

% post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0;
postStop = 0.8;


preStimAverages = nanmean(StimData(:, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
postStimAverages = nanmean(StimData(:, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
stimResponses = postStimAverages - preStimAverages;
stimContrasts = abs(BehData(:,2));


% ---- stim contrast v stim response regression analysis -----------
stimRegStats = regstats( stimContrasts, stimResponses);
StimContrastFstat = stimRegStats.fstat.f
StimContrastFstatPVal = stimRegStats.fstat.pval

%% section 4.1: ANOVA for response to STIM broken by REWARD (large/small/none)

% post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.4;
postStop = 0.8;


%data for large and small reward trials:
largeRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
smallRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);
largeRewpreStimAverages = nanmean(StimData(largeRew, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
largeRewpostStimAverages = nanmean(StimData(largeRew, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
smallRewpreStimAverages = nanmean(StimData(smallRew, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
smallRewpostStimAverages = nanmean(StimData(smallRew, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

% data for error trials:
error = find(BehData(:,9)==0);
errorPreStimAverages = nanmean(StimData(error, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
errorPostStimAverages = nanmean(StimData(error, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
 


largeRewStimAverages = largeRewpostStimAverages - largeRewpreStimAverages;
smallRewStimAverages = smallRewpostStimAverages - smallRewpreStimAverages;
errorStimAverages = errorPostStimAverages - errorPreStimAverages;

[outcome] = padcat(largeRewStimAverages, smallRewStimAverages, errorStimAverages);
[p] = anova1(outcome)


%% section 4.2: ANOVA for response to OUTCOME broken by REWARD (large/small/none)

% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.3;
postStop = 0.7;


%data for large and small reward trials:
largeRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
smallRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);
largePreRewAverages = nanmean(RewardData(largeRew, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
largePostRewAverages = nanmean(RewardData(largeRew, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
smallPreRewAverages = nanmean(RewardData(smallRew, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
smallPostRewAverages = nanmean(RewardData(smallRew, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);

% data for error trials: s0 
error = find(BehData(:,9)==0);
errorPreRewAverages = nanmean(RewardData(error, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
errorPostRewAverages = nanmean(RewardData(error, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);
 


largeRewAverages = largePostRewAverages - largePreRewAverages;
smallRewAverages = smallPostRewAverages -  smallPreRewAverages;
errorRewAverages = errorPostRewAverages - errorPreRewAverages;

[outcome] = padcat(largeRewAverages, smallRewAverages, errorRewAverages);
[p] = anova1(outcome)



%% section 5.1: point-by-point non parametric comparison of STIMULUS response for each contrast, large vs small reward 

% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.3;
postStop = 0.7;

stimz = unique (abs(BehData(:,2)));

%find index of large reward trials with -1 contrast 
largeRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
smallRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);

%pre-define responses. these will be organised with 0 contrast in first column and largest contrast in 4th column
largeRewResponses = [];
smallRewResponses = [];


for istim = stimz'

    tempIndex = intersect(find(abs(BehData(:,2))==istim),largeRew);
    addThis = nanmean(StimData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(StimData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [largeRewResponses] = padconcatenation(largeRewResponses, addThis, 2);
    
    tempIndex = intersect(find(abs(BehData(:,2))==istim),smallRew);
    addThis = nanmean(StimData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(StimData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [smallRewResponses] = padconcatenation(smallRewResponses, addThis, 2);
    
end


stimContPVals = nan(2, length(stimz));
stimContPVals(1,:) = stimz;



for istim = 1:length(stimz)
    
    [p] = ranksum(smallRewResponses(:,istim), largeRewResponses(:,istim));
    stimContPVals(2,istim) = p;
    
    
end


%% section 5.2: point-by-point non parametric comparison of REWARD response for each contrast, large vs small reward 

% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.3;
postStop = 0.7;

stimz = unique (abs(BehData(:,2)));

%find index of large reward / small reward trials
largeRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
smallRew = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);

%pre-define responses. these will be organised with 0 contrast in first column and largest contrast in 4th column
largeRewResponses = [];
smallRewResponses = [];


for istim = stimz'

    tempIndex = intersect(find(BehData(:,2)==istim),largeRew);
    addThis = nanmean(RewardData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(RewardData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [largeRewResponses] = padconcatenation(largeRewResponses, addThis, 2);
    
    tempIndex = intersect(find(BehData(:,2)==istim),smallRew);
    addThis = nanmean(RewardData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(RewardData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [smallRewResponses] = padconcatenation(smallRewResponses, addThis, 2);
    
end


stimContPVals = nan(2, length(stimz));
stimContPVals(1,:) = stimz;



for istim = 1:length(stimz)
    
    [p] = ranksum(smallRewResponses(:,istim), largeRewResponses(:,istim));
    stimContPVals(2,istim) = p;
    
    
end



%% section 5.3: point-by-point non parametric comparison of STIMULUS response for correct vs error

% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.3;
postStop = 0.7;

stimz = unique (abs(BehData(:,2)));

%find index of correct and error trials
correct = find(BehData(:,9)==1);
error = find(BehData(:,9)==0);

%pre-define responses. these will be organised with 0 contrast in first column and largest contrast in 4th column
correctResponses = [];
errorResponses = [];


for istim = stimz'

    tempIndex = intersect(find(BehData(:,2)==istim),correct);
    addThis = nanmean(StimData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(StimData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [correctResponses] = padconcatenation(correctResponses, addThis, 2);
    
    tempIndex = intersect(find(BehData(:,2)==istim),error);
    addThis = nanmean(StimData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(StimData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [errorResponses] = padconcatenation(errorResponses, addThis, 2);
    
end


stimOutcomePVals = nan(2, length(stimz));
stimOutcomePVals(1,:) = stimz;



for istim = 1:length(stimz)
    
    [p] = ranksum(errorResponses(:,istim), correctResponses(:,istim));
    stimOutcomePVals(2,istim) = p;
    
    
end

%% section 5.4: point-by-point non parametric comparison of OUTCOME response for correct vs error

% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.3;
postStop = 0.7;

stimz = unique (abs(BehData(:,2)));

%find index of correct and error trials
correct = find(BehData(:,9)==1);
error = find(BehData(:,9)==0);

%pre-define responses. these will be organised with 0 contrast in first column and largest contrast in 4th column
correctResponses = [];
errorResponses = [];


for istim = stimz'

    tempIndex = intersect(find(BehData(:,2)==istim),correct);
    addThis = nanmean(RewardData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(RewardData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [correctResponses] = padconcatenation(correctResponses, addThis, 2);
    
    tempIndex = intersect(find(BehData(:,2)==istim),error);
    addThis = nanmean(RewardData(tempIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(RewardData(tempIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
    [errorResponses] = padconcatenation(errorResponses, addThis, 2);
    
end


rewOutcomePVals = nan(2, length(stimz));
rewOutcomePVals(1,:) = stimz;



for istim = 1:length(stimz)
    
    [p] = ranksum(errorResponses(:,istim), correctResponses(:,istim));
    rewOutcomePVals(2,istim) = p;
    
    
end

%% Section 6.1: non-parametric comparison of regression slopes for contrast - dependent responses to stimulus for large vs small reward trials

clear all
close all

% list of animals
Animals = [48 50 51]

% Animals = 48

sample_rate = 12000;                                        % photoM recording sampling rate
downSample = 1200;                                       % factor downsampling the Ca responses
eventOnset = 3700;  % in the saved matrix, the ev

RTLimit = 6; % in s, excluding trials with RT longer than this


% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.7;
postStop = 1.9;


load('BehPhotoM_Exp23')

BehData=[];
StimData=[];
RewardData=[];

SRewardSlopes = [];
LRewardSlopes = [];

LResponses = [];
SResponses = [];

% get all sessions for each animal

for animal_ID = Animals
    sessionz = 1:length(BehPhotoM(animal_ID).Session);

    for iSession = sessionz

        % for each of these, the event is at column 3700
        % getting all the data from this session
        BehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        StimData   = BehPhotoM(animal_ID).Session(iSession).NeuronStim;
        
        RT = BehData(:,10) - BehData(:,13);

        % remove trials where RT > RTLimit (defined above)
        toRemove = find ( RT > RTLimit);
        BehData(toRemove,:) = [];
        StimData(toRemove,:) = [];

        
        % index for large and small reward trials
        largeRewIndex = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
        smallRewIndex = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);
        
        % stim contrast vectors for large and small reward trials
        LStimz = abs(BehData(largeRewIndex, 2));
        SStimz = abs(BehData(smallRewIndex, 2));

        % response vectors for large reward and small reward trials
        LResponses = nanmean(StimData(largeRewIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2);% - nanmean(StimData(largeRewIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
        SResponses = nanmean(StimData(smallRewIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2); %- nanmean(StimData(smallRewIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);

        % find regression slope for large and small reward trials for this session
        stats = regstats( LResponses, LStimz);
        LSlope = stats.beta(2);
        stats = regstats( SResponses, SStimz);
        SSlope = stats.beta(2);
        
        % and finally add these to the vector of all slopes for L/S trials
        [LRewardSlopes] = [LRewardSlopes; LSlope];
        [SRewardSlopes] = [SRewardSlopes; SSlope];
    end


    
    end



% now we have all the large trial regression slops and all the small trial
% regression slopes. So now we just do a rank sum to find whether they are
% the same. 


[largesmallStimPVal, largesmallStimNullReject] = ranksum(LRewardSlopes, SRewardSlopes)


%% Section 6.2: non-parametric comparison of regression slopes for contrast - dependent responses to outcome for large vs small reward trials

clear all
close all

% list of animals
Animals = [48 50 51]

% Animals = 48

sample_rate = 12000;                                        % photoM recording sampling rate
downSample = 1200;                                       % factor downsampling the Ca responses
eventOnset = 3700;  % in the saved matrix, the ev

RTLimit = 6; % in s, excluding trials with RT longer than this

% post event time windows (seconds):
preStart = -0.4;
preStop = 0;
postStart = 0.3;
postStop = 1.3;


load('BehPhotoM_Exp23')

BehData=[];
StimData=[];
RewardData=[];

SRewardSlopes = [];
LRewardSlopes = [];

    LResponses = [];
    SResponses = [];

% get all sessions for each animal

for animal_ID = Animals
    sessionz = 1:length(BehPhotoM(animal_ID).Session);

    for iSession = sessionz

        % for each of these, the event is at column 3700
        % getting all the data from this session
        BehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        RewardData = BehPhotoM(animal_ID).Session(iSession).NeuronReward;
        
                
        RT = BehData(:,10) - BehData(:,13);

        % remove trials where RT > RTLimit (defined above)
        toRemove = find ( RT > RTLimit);
        BehData(toRemove,:) = [];
        RewardData(toRemove,:) = [];

         
        % index for large and small reward trials
        largeRewIndex = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==2)))]);
        smallRewIndex = sort([(intersect(find(BehData(:,9)==1 & BehData(:,3)==1), find(BehData(:,8)==1))); (intersect(find(BehData(:,9)==1 & BehData(:,3)==-1), find(BehData(:,8)==2)))]);
        
        % stim contrast vectors for large and small reward trials
        LStimz = abs(BehData(largeRewIndex, 2));
        SStimz = abs(BehData(smallRewIndex, 2));

        % response vectors for large reward and small reward trials
        LResponses = nanmean(RewardData(largeRewIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(RewardData(largeRewIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
        SResponses = nanmean(RewardData(smallRewIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2) - nanmean(RewardData(smallRewIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);

        % find regression slope for large and small reward trials for this session
        stats = regstats( LResponses, LStimz);
        LSlope = stats.beta(2);
        stats = regstats( SResponses, SStimz);
        SSlope = stats.beta(2);
        
        % and finally add these to the vector of all slopes for L/S trials
        [LRewardSlopes] = [LRewardSlopes; LSlope];
        [SRewardSlopes] = [SRewardSlopes; SSlope];
    end


    
    end
    


% now we have all the large trial regression slops and all the small trial
% regression slopes. So now we just do a rank sum to find whether they are
% the same. 


[largesmallOutcomePVal, largesmallOutcomeNullReject] = ranksum(LRewardSlopes, SRewardSlopes)



%% Section 6.3: non-parametric comparison of regression slopes for contrast - dependent responses to stimulus for correct vs error trials

clear all
close all

% list of animals
Animals = [48 50 51]

% Animals = 48

sample_rate = 12000;                                        % photoM recording sampling rate
downSample = 1200;                                       % factor downsampling the Ca responses
eventOnset = 3700;  % in the saved matrix, the ev

RTLimit = 6; % in s, excluding trials with RT longer than this


% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.7;
postStop = 1.9;


load('BehPhotoM_Exp23')

BehData=[];
StimData=[];

CorrectSlopes = [];
ErrorSlopes = [];

CorrectResponses = [];
ErrorResponses = [];
CorrectPVals = [];
ErrorPVals = [];

% get all sessions for each animal

for animal_ID = Animals
    sessionz = 1:length(BehPhotoM(animal_ID).Session);

    for iSession = sessionz

        % for each of these, the event is at column 3700
        % getting all the data from this session
        BehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        StimData   = BehPhotoM(animal_ID).Session(iSession).NeuronStim;
        
        RT = BehData(:,10) - BehData(:,13);

        % remove trials where RT > RTLimit (defined above)
        toRemove = find ( RT > RTLimit);
        BehData(toRemove,:) = [];
        StimData(toRemove,:) = [];

        
        % index for correct and error trials
        correctIndex = find(BehData(:,9)==1);
        errorIndex = find(BehData(:,9)==0);
        
        % stim contrast vectors for large and small reward trials
        CorrectStimz = abs(BehData(correctIndex, 2));
        ErrorStimz = abs(BehData(errorIndex, 2));

        % response vectors for large reward and small reward trials
        CorrectResponses = nanmean(StimData(correctIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2); - nanmean(StimData(correctIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
        ErrorResponses = nanmean(StimData(errorIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2); - nanmean(StimData(errorIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);

        % find regression slope for large and small reward trials for this session
        stats = regstats( CorrectResponses, CorrectStimz);
        CSlope = stats.beta(2);
        CPVal = stats.fstat.pval;
        stats = regstats( ErrorResponses, ErrorStimz);
        ESlope = stats.beta(2);
        EPVal = stats.fstat.pval;
        
        % and finally add these to the vector of all slopes for L/S trials
        [CorrectSlopes] = [CorrectSlopes; CSlope];
        [ErrorSlopes] = [ErrorSlopes; ESlope];
        [CorrectPVals] = [CorrectPVals; CPVal];
        [ErrorPVals] = [ErrorPVals; EPVal];
    end


    
    end



% now we have all the large trial regression slops and all the small trial
% regression slopes. So now we just do a rank sum to find whether they are
% the same. 


[correrrStimPVal, correrrStimNullReject] = ranksum(CorrectSlopes, ErrorSlopes)


%% Section 6.4: non-parametric comparison of regression slopes for contrast - dependent responses to outcome for correct vs error trials

clear all
close all

% list of animals
Animals = [48 50 51]

% Animals = 48

sample_rate = 12000;                                        % photoM recording sampling rate
downSample = 1200;                                       % factor downsampling the Ca responses
eventOnset = 3700;  % in the saved matrix, the ev

RTLimit = 6; % in s, excluding trials with RT longer than this


% post event time windows (seconds):
preStart = -0.2;
preStop = 0;
postStart = 0.3;
postStop = 0.9;


load('BehPhotoM_Exp23')

BehData=[];
RewardData=[];

CorrectSlopes = [];
ErrorSlopes = [];

CorrectResponses = [];
ErrorResponses = [];
CorrectPVals = [];
ErrorPVals = [];

% get all sessions for each animal

for animal_ID = Animals
    sessionz = 1:length(BehPhotoM(animal_ID).Session);

    for iSession = sessionz

        % for each of these, the event is at column 3700
        % getting all the data from this session
        BehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        RewardData = BehPhotoM(animal_ID).Session(iSession).NeuronReward;
        
        RT = BehData(:,10) - BehData(:,13);

        % remove trials where RT > RTLimit (defined above)
        toRemove = find ( RT > RTLimit);
        BehData(toRemove,:) = [];
        RewardData(toRemove,:) = [];

        
        % index for correct and error trials
        correctIndex = find(BehData(:,9)==1);
        errorIndex = find(BehData(:,9)==0);
        
        % stim contrast vectors for large and small reward trials
        CorrectStimz = abs(BehData(correctIndex, 2));
        ErrorStimz = abs(BehData(errorIndex, 2));

        % response vectors for large reward and small reward trials
        CorrectResponses = nanmean(RewardData(correctIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2); - nanmean(RewardData(correctIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);
        ErrorResponses = nanmean(RewardData(errorIndex, (eventOnset+(postStart*downSample)):(eventOnset+(postStop*downSample))), 2); - nanmean(RewardData(errorIndex, (eventOnset+(preStart*downSample)):(eventOnset+(preStop*downSample))), 2);

        % find regression slope for large and small reward trials for this session
        stats = regstats( CorrectResponses, CorrectStimz);
        CSlope = stats.beta(2);
        CPVal = stats.fstat.pval;
        stats = regstats( ErrorResponses, ErrorStimz);
        ESlope = stats.beta(2);
        EPVal = stats.fstat.pval;
        
        % and finally add these to the vector of all slopes for L/S trials
        [CorrectSlopes] = [CorrectSlopes; CSlope];
        [ErrorSlopes] = [ErrorSlopes; ESlope];
        [CorrectPVals] = [CorrectPVals; CPVal];
        [ErrorPVals] = [ErrorPVals; EPVal];
    end


    
    end



% now we have all the large trial regression slops and all the small trial
% regression slopes. So now we just do a rank sum to find whether they are
% the same. 


[correrrStimPVal, correrrStimNullReject] = ranksum(CorrectSlopes, ErrorSlopes)

