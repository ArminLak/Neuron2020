clear all
close all

% give animal name and it will plot data averaged across all session

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)


%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 50
load('BehPhotoM_Exp23')

RTLimit = 7; % in s, Dont change. excluding trials with RT longer than this
Timereso = 40; %this is in unit of 1200 in s. 
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
RewardData = [];

Kernel_beep=nan(1,30);
Kernel_stim=nan(1,30);
Kernel_action=nan(1,30);
Kernel_outcome=nan(1,30);


sessionz = 1:length(BehPhotoM(animal_ID).Session);
%if animal_ID == 51
%sessionz = [1 2 3 4 5 8]; end % sessions used for ALK071 Exp 23 [1 2 3 4 5 8]

for iSession = sessionz
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStim;
    
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;
    
    BeepData = [BeepData;TempBeepData];
    
    StimData = [StimData;TempStimData];
    
    RewardData = [RewardData;TempRewardData];
    
    BehData = [BehData; TempBehData];
    
   
end


% trial by trial reaction time (in fact reward stim interval)
RT = BehData(:,14) - BehData(:,13);


ToLargeR = find((BehData(:,3)==-1 & BehData(:,8)==1)  | ...
    (BehData(:,3)==1 & BehData(:,8)==2));

BehData(ToLargeR,16)=1;

ToSmallR = setdiff(1:size(BehData,1),ToLargeR)';

BehData(ToSmallR, 16)=-1;

abzStim = unique(abs(BehData(:,2)))';

totalSamples = size(StimData,1) * size(StimData,2); 


toExclude = unique( [find(BehData(:,1)==1) ; find(RT > RTLimit)]); % exlcuding trials with very long RT
%and also the first trials of each sessnio since the beep timing in the
%very first trial is sometimes wrong

% make continous data
 inverted_Matrix=StimData;
 inverted_Matrix(toExclude,:) = [];
 inverted_Matrix =inverted_Matrix';
 invertedStimData_to_row = inverted_Matrix(:)';
 
% make continous time 
t = 1:1:totalSamples; % time
t=reshape(t,13100,size(BehData,1))';
t(toExclude,:) = [];
t=reshape(t',length(invertedStimData_to_row),1);


% this is according to my saved data (13100 samples convering -3 to 8 s
% after the event)
stim_time =  3700 : 13100 : totalSamples;
beep_time        = stim_time - 1200.* (BehData (:,13) - BehData (:,12))';
action_onsetTime = stim_time + 1200.* (BehData (:,10) - BehData (:,13))';
outcome_time     = stim_time + 1200.* (BehData (:,14) - BehData (:,13))';

beep_time(toExclude) = [];
stim_time (toExclude) = [];
action_onsetTime(toExclude) = [];
outcome_time(toExclude) = [];
action_onsetTime (isnan(action_onsetTime))=outcome_time(isnan(action_onsetTime))-200 ; % very rare trials we could not measure action_onset time and we use this one


BehData(toExclude,:) = [];
BeepData(toExclude,:) = [];
StimData(toExclude,:) = [];
RewardData(toExclude,:) = [];

% moving average and downsampling
invertedStimData_to_row = movingmean(invertedStimData_to_row',100)';
inverted_Matrix_to_row = downsample(invertedStimData_to_row',Timereso);
t=downsample(t,Timereso);

% kernle regression
eventTimes = {beep_time;stim_time; action_onsetTime; outcome_time};
eventValues = {[];[];[];[]};
windows = {[-400 800];[-300 900];[-800 400];[-200 1000]}; % 

[fitKernels, predictedSignals,cvErr] = kernelRegression(inverted_Matrix_to_row', t', eventTimes, eventValues, windows , [5], [0 0]);

%%
% find trial by trail coefficents (regression)

StimDataSmooth=movingmean(StimData,100,2);
StimDataSmoothDownSample=downsample(StimDataSmooth',Timereso)';

    Coefiz = nan(size(StimData,1),4);
 
for itrial = 1: size(StimData,1)
    
    zeroPadBeep = zeros(1,328);
     
     zeroPadAction = zeros(1,328);
        zeroPadStim = zeros(1,328);
        zeroPadOutcome = zeros(1,328);
        
        % 90 is the stim onset 
        zeroPadStim (90 -8:90+29-8 ) = fitKernels{2}';
        
        
        t_beep = floor((BehData(itrial, 13 ) - BehData(itrial, 12))*30); % 11 s is 328 samples, thus 1 s is 30
        
        zeroPadBeep(90-t_beep -10 : 90 -t_beep+ 29-10) = fitKernels{1}';
        
        t_act = floor((BehData(itrial, 10 ) - BehData(itrial, 13))*30); % 11 s is 328 samples, thus 1 s is 30
        
        t_act (isnan(t_act))=floor((BehData(itrial, 14 ) - BehData(itrial, 13))*30) - 5; % very rare trials we could not measure action_onset time and we use this one

        
        zeroPadAction(90+t_act-10 : 90+t_act-10 +29) = fitKernels{3}';
        
        t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*30); % 11 s is 328 samples, thus 1 s is 30
        
        zeroPadOutcome(90+t_out-5 : 90+t_out-5 +29) = fitKernels{4}';
        
        
           Coefiz(itrial,:)= regress(StimDataSmoothDownSample(itrial,:)',[zeroPadBeep',zeroPadStim',zeroPadAction',zeroPadOutcome'])';

           
           
end
%%
Coefiz(Coefiz > 10) = nan;
Coefiz(Coefiz < -10) = nan;


StimsAllowed = unique(BehData(:,2))';
Stimz = BehData(:,2);
block = BehData(:,8);
correct = BehData(:,9);

for i=unique(BehData(:,2))'
        
        tempIndex = find(Stimz==i);
        
        [i j]=ismember(i, StimsAllowed);
        
        
        
        resp_Block1StimCor(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),1));
        resp_Block2StimCor(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),1));
        
        resp_Block1StimErr(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==0)),1));
        resp_Block2StimErr(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==0)),1));
        
        resp_Block1ActCor(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),2));
        resp_Block2ActCor(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),2));
        
        resp_Block1ActErr(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==0)),2));
        resp_Block2ActErr(j) = nanmean(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==0)),2));
        
        
end
    


