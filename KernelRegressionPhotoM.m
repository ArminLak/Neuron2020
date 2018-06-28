clear all
close all

% give animal name and it will plot data averaged across all session

% it also save average data of the animal into 'GrandSummary' (you will need to save mannulay)


%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 51
load('BehPhotoM_Exp23')

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

% moving average and downsampling
%invertedStimData_to_row = movingmean(invertedStimData_to_row',100)';
%inverted_Matrix_to_row = (inverted_Matrix_to_row -
%movmean(inverted_Matrix_to_row,3000)); this was a bad idea


%  bpFilt =   designfilt('highpassfir','FilterOrder',20,'CutoffFrequency',50, 'SampleRate', 1200);




% inverted_Matrix_to_row = smooth(invertedStimData_to_row,300);

% StimDataSmooth = reshape(inverted_Matrix_to_row', 13100,size(BehData,1))';
% StimDataSmoothDownSample=downsample(StimDataSmooth',Timereso)';
%
% [inverted_Matrix_to_row, s] = deconvolveCa(inverted_Matrix_to_row,'ar2');
%
%
% inverted_Matrix_to_row = downsample(inverted_Matrix_to_row,Timereso);



inverted_Matrix_to_row = invertedStimData_to_row;
inverted_Matrix_to_row = downsample(inverted_Matrix_to_row,Timereso);
StimDataSmooth2Visualise = reshape(inverted_Matrix_to_row', 262,size(BehData,1))';


[inverted_Matrix_to_row, s] = deconvolveCa(inverted_Matrix_to_row,'ar2');

StimDataSmoothDownSample = reshape(inverted_Matrix_to_row', 262,size(BehData,1))';

t=downsample(t,Timereso);


% kernle regression
%%
% full
eventTimes = {beep_time;stim_time; action_onsetTime; outcome_time};
eventValues = {[];[];[];[]};

windows = {[-400 800];[-200 1600];[-600 200];[200 2000]}; %
windows = {[-400 800];[-400 2600];[-1000 200];[-400 2000]}; %

%windows = {[-1000 2000];[-1000 3000];[-2000 600];[-800 2000]}; %

%    eventTimes = {beep_time; action_onsetTime; outcome_time};
%    eventValues = {[];[];[]};
% % % % windows = {[-400 800];[-600 200];[200 2000]}; %
%   windows = {[-400 800];[-1000 200];[-400 2000]}; %
% % %
% % % %
%      eventTimes = {beep_time;stim_time;  outcome_time};
%      eventValues = {[];[];[]};
% % % %   windows = {[-400 800];[-200 1600];[200 2000]}; %
%   windows = {[-400 800];[-400 2600];[-400 2000]}; %
% % %
% % %
%   eventTimes = {beep_time;stim_time; action_onsetTime};
%   eventValues = {[];[];[]};
% % %
%   windows = {[-400 800];[-400 2600];[-1000 200]}; %
% %
% %
%  eventTimes = {stim_time; action_onsetTime; outcome_time};
% % eventValues = {[];[];[]};
%  %windows = {[-200 1600];[-600 200];[200 2000]};
%  windows = {[-400 2600];[-1000 200];[-400 2000]}; %


[fitKernels, predictedSignals,cvErr] = kernelRegression(inverted_Matrix_to_row', t', eventTimes, eventValues, windows , [0], [0 0]);

% explained Variance
EV = 1-mean(mean((predictedSignals-inverted_Matrix_to_row').^2))/mean(mean(inverted_Matrix_to_row.^2));

for iter = 1 : length(stim_time) -1  % we ingore the last trial
    
    [diff_recover i_recover] = min (abs(t - stim_time(iter)));
    
    
    Sliced_PredictedSignalStim (iter,:) =predictedSignals(i_recover-10 : i_recover+200);
    
    [diff_recover i_recover] = min (abs(t - outcome_time(iter)));
    
    
    Sliced_PredictedSignalRew (iter,:) =predictedSignals(i_recover-10 : i_recover+200);
    
    
end




figure

subplot(1,4,1)
plot(fitKernels{1})

subplot(1,4,2)
plot(fitKernels{2})

subplot(1,4,3)
plot(fitKernels{3})

subplot(1,4,4)
plot(fitKernels{4})



%%
% find trial by trail coefficents (regression)



%StimDataSmoothDownSample = smooth2a(StimDataSmoothDownSample,10,0);

scaling =1;

Coefiz = nan(size(StimData,1),4); % full
%Coefiz = nan(size(StimData,1),3); % reduced, no action


for itrial = 1: size(StimData,1)
    
    % 72 is the stim onset (for 262) - full model
    
    
    zeroPadBeep = zeros(1,262);
    zeroPadAction = zeros(1,262);
    zeroPadStim = zeros(1,262);
    zeroPadOutcome = zeros(1,262);
    
    %zeroPadStim (71 -4:71+35-4) = scaling .* (fitKernels{2}');
    zeroPadStim (71 -8:71+59-8) = scaling .* (fitKernels{2}');
    
    t_beep = floor((BehData(itrial, 13 ) - BehData(itrial, 12))*24); % 11 s is 262 samples, thus 1 s is 24
    zeroPadBeep(71-t_beep -8 : 71 -t_beep+ 23-8) = scaling.*( fitKernels{1}');
    
    
    t_act = floor((BehData(itrial, 10 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
    t_act (isnan(t_act))=floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24) - 4; % very rare trials we could not measure action_onset time and we use this one
    %zeroPadAction(71+t_act-12 : 71+t_act-12 +15) = scaling .*(fitKernels{3}');
    zeroPadAction(71+t_act-20 : 71+t_act-20 +23) = scaling .*(fitKernels{3}');
    
    t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
    % zeroPadOutcome(71+t_out+4 : 71+t_out+4 +35) = scaling .*(fitKernels{4}');
    zeroPadOutcome(71+t_out-8: 71+t_out-8 +47) = scaling .*(fitKernels{4}');
    
    [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadBeep',zeroPadStim',zeroPadAction',zeroPadOutcome']);
    
    
    % reduced model with good EV
    %      zeroPadBeep = zeros(1,262);
    %     zeroPadStim = zeros(1,262);
    %     zeroPadOutcome = zeros(1,262);
    %
    %     zeroPadStim (72 -4:72+35-4) = scaling .* (fitKernels{2}');
    %
    %     t_beep = floor((BehData(itrial, 13 ) - BehData(itrial, 12))*24); % 11 s is 262 samples, thus 1 s is 24
    %     zeroPadBeep(72-t_beep -8 : 72 -t_beep+ 23-8) = scaling.*( fitKernels{1}');
    %
    %     t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
    %     zeroPadOutcome(72+t_out+4 : 72+t_out+4 +35) = scaling .*(fitKernels{3}');
    %
    %
    %    [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadBeep',zeroPadStim',zeroPadOutcome']);
    
    
    Coefiz(itrial,:)= B';
    
    
    
    % plotting example trials:
    %     if ismember(itrial,[400:430])
    %
    %         figure
    %         plot(smooth(sum([B(1) .* zeroPadBeep',B(2) .*zeroPadStim',B(3) .*zeroPadAction',B(4) .*zeroPadOutcome']'),20) )
    %         hold on
    %         plot(smooth(StimDataSmooth2Visualise(itrial,:)',20))
    %
    %
    %
    %         ax = gca;
    %
    %         InterStimAction = 71 +  t_act;
    %         h=rectangle(ax, 'Position',[InterStimAction 0  0.1 5],'EdgeColor',[1 0 0]);
    %
    %
    %         InterBeepStim = 71 -  t_beep ;
    %         h=rectangle(ax, 'Position',[InterBeepStim 0  0.1 5],'EdgeColor',[0 1 0]);
    %
    %
    %         InterStimRew = 71 +  t_out ;
    %         h=rectangle(ax, 'Position',[InterStimRew 0  0.1 5],'EdgeColor',[0 0 1]);
    %
    %     end
    
    
end
%%
%Coefiz(Coefiz > 3) = nan;
%Coefiz(Coefiz < -2) = nan;

%Coefiz = zscore(Coefiz')';


StimsAllowed = unique(BehData(:,2))';
Stimz = BehData(:,2);
block = BehData(:,8);
correct = BehData(:,9);
for i=unique(BehData(:,2))'
    
    tempIndex = find(Stimz==i);
    
    [i j]=ismember(i, StimsAllowed);
    
    
    % full
    resp_CorrStim(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1)),2));
    resp_ErrStim(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0)),2));
    
    
    resp_CorrAct(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1)),3));
    resp_ErrAct(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0)),3));
    
    resp_CorrRew(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1)),4));
    resp_ErrRew(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0)),4));
    
    resp_Block1StimCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),2));
    resp_Block2StimCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),2));
    
    resp_Block1ActCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==0)),3));
    resp_Block2ActCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==0)),3));
    
    resp_Block1RewCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),4));
    resp_Block2RewCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),4));
    
    
    % reduced
    %     resp_CorrStim(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1)),2));
    %     resp_ErrStim(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0)),2));
    %
    %
    %     resp_CorrRew(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1)),3));
    %     resp_ErrRew(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0)),3));
    %
    %     resp_Block1StimCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),2));
    %     resp_Block2StimCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),2));
    %
    %     resp_Block1RewCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==1),find(    correct==1)),3));
    %     resp_Block2RewCor(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    block==2),find(    correct==1)),3));
    %
    
end

figure

subplot(2,3,1)
plot(resp_CorrStim,'g')
hold on
plot(resp_ErrStim,'r')

subplot(2,3,2)

% plot(resp_CorrAct,'g')
% hold on
% plot(resp_ErrAct,'r')

subplot(2,3,3)


plot(resp_CorrRew,'g')
hold on
plot(resp_ErrRew,'r')

subplot(2,3,4)
plot(resp_Block1StimCor,'k')
hold on
plot(resp_Block2StimCor,'-.k')

subplot(2,3,5)
% plot(resp_Block1ActCor,'k')
% hold on
% plot(resp_Block2ActCor,'-.k')



subplot(2,3,6)

plot(resp_Block1RewCor,'k')
hold on
plot(resp_Block2RewCor,'-.k')


%save
BehPhotoM(animal_ID).KernelSummary.Kernels = fitKernels;
BehPhotoM(animal_ID).KernelSummary.resp_CorrStim = resp_CorrStim;
BehPhotoM(animal_ID).KernelSummary.resp_ErrStim = resp_ErrStim;

BehPhotoM(animal_ID).KernelSummary.resp_CorrAct=resp_CorrAct;
BehPhotoM(animal_ID).KernelSummary.resp_ErrAct=resp_ErrAct;

BehPhotoM(animal_ID).KernelSummary.resp_CorrRew=resp_CorrRew;
BehPhotoM(animal_ID).KernelSummary.resp_ErrRew=resp_ErrRew;

BehPhotoM(animal_ID).KernelSummary.resp_Block1StimCor=resp_Block1StimCor;
BehPhotoM(animal_ID).KernelSummary.resp_Block2StimCor=resp_Block2StimCor;

BehPhotoM(animal_ID).KernelSummary.resp_Block1ActCor=resp_Block1ActCor;
BehPhotoM(animal_ID).KernelSummary.resp_Block2ActCor=resp_Block2ActCor;

BehPhotoM(animal_ID).KernelSummary.resp_Block1RewCor=resp_Block1RewCor;
BehPhotoM(animal_ID).KernelSummary.resp_Block2RewCor=resp_Block2RewCor;

%%
StimzAbs = abs(BehData(1:end-1,2));

StimsAllowed = unique(StimzAbs');
for i=unique(StimzAbs')
    
    tempIndex = find(StimzAbs==i);
    
    [i j]=ismember(i, StimsAllowed);
    
    
    PredictedPopStimAlign(j,:)=nanmean(Sliced_PredictedSignalStim(tempIndex,:));
    
end

figure;

coefStim = [0.6 1 2 3]


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
plot(nanmean(Sliced_PredictedSignalRew(IndexLargeRew,:)));

plot(nanmean(0.5 * Sliced_PredictedSignalRew(IndexSmallRew,:)));

plot(nanmean(-0.2 * Sliced_PredictedSignalRew(IndexNoRew,:)));




BehPhotoM(animal_ID).KernelSummary.PredictedPopStimAlign=PredictedPopStimAlign;
BehPhotoM(animal_ID).KernelSummary.PredictedPopRewAlignLarge=nanmean(Sliced_PredictedSignalRew(IndexLargeRew,:));
BehPhotoM(animal_ID).KernelSummary.PredictedPopRewAlignSmall=nanmean(0.5 * Sliced_PredictedSignalRew(IndexSmallRew,:));
BehPhotoM(animal_ID).KernelSummary.PredictedPopRewAlignNoRew=nanmean(-0.2 *Sliced_PredictedSignalRew(IndexNoRew,:));



