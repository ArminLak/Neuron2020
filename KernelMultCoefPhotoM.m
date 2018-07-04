clear all
close all

% give animal name and it will plot data averaged across all session

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)


%[48, 50,51]  coresponding to ALK068, 70 and 71

load('BehPhotoM_Exp23')
%
animal_ID = 50

ModelArrangment = 11

RTLimit = 5.9; % in s, Dont change. excluding trials with RT longer than this
Timereso = 50; %this is in unit of 1200 in s. dont change
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


BehData = [];
BeepData = [];
StimData = [];
RewardData = [];


sessionz = 1:length(BehPhotoM(animal_ID).Session);


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

%StimData(:,5400:end) = StimData(:,5400:end) - repmat(mean(StimData (:,5400:5500),2), 1, 13100-5399);
%StimData(:,5500:end) = StimData(:,5500:end) - repmat(mean(StimData (:,5200:5500),2), 1, 13100-5499);

%StimData = StimData - repmat(mean(StimData,2), 1, 13100);

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

inverted_Matrix_to_row = invertedStimData_to_row;
inverted_Matrix_to_row = downsample(inverted_Matrix_to_row,Timereso);
StimDataSmooth2Visualise = reshape(inverted_Matrix_to_row', 262,size(BehData,1))';


[inverted_Matrix_to_row, s] = deconvolveCa(inverted_Matrix_to_row,'ar2'); % smooo

StimDataSmoothDownSample = reshape(inverted_Matrix_to_row', 262,size(BehData,1))';

t=downsample(t,Timereso);


% kernle regression

% full
eventTimes = {stim_time; action_onsetTime; outcome_time};
eventValues = {[];[];[]};

windows = {[-400 800];[-200 1600];[-600 200];[200 2000]}; %

windows = {[-400 2600];[-1000 200];[-400 2000]}; %


% stim,outcome
eventTimes = {stim_time; outcome_time};

windows = {[-400 2600];[-400 2000]}; %


for fititer=1:1:4
    
    
    [fitKernels, predictedSignals,cvErr] = kernelRegression(inverted_Matrix_to_row', t', eventTimes, eventValues, windows , [0], [0 0]);
 
    
    Coefiz=ones(length(beep_time),3);
    
    
    % recover observed and estimated PSTHs
    for itrial = 1 : length(beep_time)-1
        
 
        zeroPadAction = zeros(1,262);
        zeroPadStim = zeros(1,262);
        zeroPadOutcome = zeros(1,262);
        
       
        if ModelArrangment ==8
            
        zeroPadStim (71 -8:71+59-8) =  (fitKernels{1}');
        
        
        t_act = floor((BehData(itrial, 10 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
        zeroPadAction(71+t_act-20 : 71+t_act-20 +23) = (fitKernels{2}');
        
       
        t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
        zeroPadOutcome(71+t_out-8: 71+t_out-8 +47) = (fitKernels{3}');
        
        [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadStim',zeroPadAction',zeroPadOutcome']);
        
        Coefiz(itrial,L)= B';
        end
        
        
                if ModelArrangment ==11

        zeroPadStim (71 -8:71+59-8) =  (fitKernels{1}');
        
       
        t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
        zeroPadOutcome(71+t_out-8: 71+t_out-8 +47) = (fitKernels{2}');
        
        [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadStim',zeroPadOutcome']);
        
          tofill = [1,3];
        Coefiz(itrial,tofill)= B';
                end
        
        
        for i=1:length(eventTimes)
            
            % this is for shuffling coef across trials
            
            % Coefiz(:,i) = Coefiz(randperm(length(Coefiz(:,i))),1);
            
            eventValues{i} = Coefiz(:,tofill(i))';
            
        end
        
        
    end
    
    % compute EV
    EV(fititer) = 1-mean(mean((predictedSignals-inverted_Matrix_to_row').^2))/mean(mean(inverted_Matrix_to_row.^2));
    
end

Coefiz = Coefiz  ./ mean (Coefiz);

SignalStim = [];
SignalAction = [];
SignalOutcome = [];

EstimatedSignalStim = [];
EstimatedSignalAction = [];
EstimatedSignalOutcome = [];


for iter = 1 : length(stim_time) -1  % we ingore the last trial
    
    [diff_recover i_recover] = min (abs(t - stim_time(iter)));
    
    
    EstimatedSignalStim (iter,:) =predictedSignals(i_recover-10 : i_recover+200);
    SignalStim (iter,:) =inverted_Matrix_to_row(i_recover-10 : i_recover+200)';
    
    
    [diff_recover i_recover] = min (abs(t - outcome_time(iter)));
    
    
    EstimatedSignalOutcome (iter,:) =predictedSignals(i_recover-10 : i_recover+200);
    SignalOutcome (iter,:) =inverted_Matrix_to_row(i_recover-10 : i_recover+200)';
    
    
end



%%
BehData(end,:)=[];
StimsAllowed = unique(BehData(:,2))';
Stimz = BehData(:,2);
block = BehData(:,8);
correct = BehData(:,9);
BehDataAbs = BehData;
BehDataAbs(:,2) = abs(BehDataAbs(:,2));

[i j] = sortrows(BehDataAbs,[9 2]);
figure

subplot(1,2,1)
imagesc(EstimatedSignalStim(j,3:30),[-1 15])
title('estimated stim')
colormap('Copper')

subplot(1,2,2)
imagesc(SignalStim(j,3:30),[-1 15])
colormap('Copper')
title('real stim')

figure
subplot(1,2,1)
imagesc(EstimatedSignalOutcome(j,2:50),[0 5])
title('estimated outcome')
colormap('Copper')

subplot(1,2,2)
imagesc(SignalOutcome(j,2:50),[0 5])
colormap('Copper')
title('real outcome')


%


for i=abzStim
    
    tempIndex = find(abs(Stimz)==i);
    
    [i j]=ismember(i, abzStim);
    
 
    
    % correct large
    resp_CorrStimLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==1)),1));
    resp_ErrStimLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0),find(BehData(:,16)==1)),1));
    
    resp_CorrStimSmall(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==-1)),1));

     resp_CorrActionLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==1)),2));
    resp_ErrActionLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0),find(BehData(:,16)==1)),2));
    
    resp_CorrActionSmall(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==-1)),2));

    resp_CorrOutcomeLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==1)),3));
    resp_ErrOutcomeLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0),find(BehData(:,16)==1)),3));
    
    resp_CorrOutcomeSmall(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==-1)),3));

    
    resp_Block1StimCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),1));
    resp_Block2StimCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),1));
    
    resp_Block1ActCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==0)),2));
    resp_Block2ActCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==0)),2));
    
    resp_Block1RewCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),3));
    resp_Block2RewCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),3));
    
   
end
%
figure

subplot(2,3,1)
plot(resp_CorrStimLarge,'g')
hold on
plot(resp_ErrStimLarge,'r')

subplot(2,3,2)

plot(resp_CorrStimLarge,'g')
hold on
plot(resp_CorrStimSmall,'--g')


% plot(resp_CorrAct,'g')
% hold on
% plot(resp_ErrAct,'r')

subplot(2,3,3)


plot(resp_CorrOutcomeLarge,'g')
hold on
plot(resp_ErrOutcomeLarge,'r')

subplot(2,3,4)
plot(resp_CorrOutcomeLarge,'g')
hold on
plot(resp_CorrOutcomeSmall,'--g')

subplot(2,3,5)
% plot(resp_Block1ActCor,'k')
% hold on
% plot(resp_Block2ActCor,'-.k')



% subplot(2,3,6)
% 
% plot(resp_Block1RewCor,'k')
% hold on
% plot(resp_Block2RewCor,'-.k')


%save
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).Kernels = fitKernels;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrStimLarge = resp_CorrStimLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrStimSmall = resp_CorrStimSmall;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrStimLarge=resp_CorrStimLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrOutcomeLarge=resp_CorrOutcomeLarge;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_ErrOutcomeLarge=resp_ErrOutcomeLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrOutcomeSmall=resp_CorrOutcomeSmall;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_Block1StimCor=resp_Block1StimCor;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_Block2StimCor=resp_Block2StimCor;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_Block1ActCor=resp_Block1ActCor;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_Block2ActCor=resp_Block2ActCor;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_Block1RewCor=resp_Block1RewCor;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_Block2RewCor=resp_Block2RewCor;

%
StimzAbs = abs(BehData(1:end-1,2));

StimsAllowed = unique(StimzAbs');
for i=unique(StimzAbs')
    
    tempIndex = find(StimzAbs==i);
    
    [i j]=ismember(i, StimsAllowed);
    
    
    PredictedPopStimAlign(j,:)=nanmean(EstimatedSignalStim(tempIndex,:));
    
end
%
figure;

coefStim = [0.6 1 2 3]

coefStim = [1 1 1 1]


for i=1:4
    
    PredictedPopStimAlign(i,15:end)=  coefStim(i) .* PredictedPopStimAlign(i,15:end);
    hold on
    plot(PredictedPopStimAlign(i,:),'color',colorGray(i,:))
    
end

BehData(end,:)=[];
IndexLargeRew = [find(BehData(:,9)==1 & BehData(:,3)==-1 & BehData(:,8)==1 ) ; find(BehData(:,9)==1 & BehData(:,3)==1 & BehData(:,8)==2) ];
IndexSmallRew = [find(BehData(:,9)==1 & BehData(:,3)==-1 & BehData(:,8)==2 ) ; find(BehData(:,9)==1 & BehData(:,3)==1 & BehData(:,8)==1) ];
IndexNoRew = find(BehData(:,9)==1);

figure; hold on
plot(nanmean(EstimatedSignalOutcome(IndexLargeRew,:)));

plot(nanmean(EstimatedSignalOutcome(IndexSmallRew,:)));

plot(nanmean(EstimatedSignalOutcome(IndexNoRew,:)));




BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopStimAlign=PredictedPopStimAlign;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopRewAlignLarge=nanmean(EstimatedSignalOutcome(IndexLargeRew,:));
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopRewAlignSmall=nanmean(EstimatedSignalOutcome(IndexSmallRew,:));
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopRewAlignNoRew=nanmean(EstimatedSignalOutcome(IndexNoRew,:));

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).EV=EV;


