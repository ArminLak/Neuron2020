clear all
close all

% last checked Jan 2019. it includes figures and stat for the paper.

% give animal name and it estimate kernels and their coefficents. it  plots mega raster for the paper. For this change the following:
% put RTLimit = 3
% comment around line 85
% set model to 11
% in the selectkernel code, set windows for the plotting
% under section dealing with model 11 (in this code), select the line which
% has correct windows duration

% model 11: for real analysis 
% model 11 with modified RT, and kernel duration for mega raster in the
% paper
% model 15 for pop PSTH of model prediction

%% includes regression stats for regressing coefficents on chosenValue and
%% RPE

%[48, 50,51]  coresponding to ALK068, 70 and 71

load('BehPhotoM_Exp23_VTA')
%load('BehPhotoM_Exp23_NAc.mat')

%
animal_ID = 51

if animal_ID ==48
animal_name = 'ALK068'
elseif animal_ID ==50
animal_name = 'ALK070'
elseif animal_ID ==51
animal_name = 'ALK071'
elseif animal_ID ==56
animal_name = 'ALK078'

elseif animal_ID ==57
animal_name = 'MMM001'

elseif animal_ID ==56
animal_name = 'MMM002'

end
    

ModelArrangment = 11 %  ModelArrangement=15 (it is just similar to model 11) for visualing pop PSTH of the model predictions 


RTLimit = 5.9; % in s, Dont change. excluding trials with RT longer than this

% for visualsing mega rasters only 
RTLimit = 3; % for visualisng in s, Dont change. excluding trials with RT longer than this


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
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepR;
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardR;
    
    if animal_ID ==56 % for ALK078 we use right hem data 
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepR;
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardR;
    end
    
    
    BeepData = [BeepData;TempBeepData];
    StimData = [StimData;TempStimData];
    RewardData = [RewardData;TempRewardData];
    BehData = [BehData; TempBehData];
    
end

figure
plot(nanmean(StimData))



% for mega raster of paper, comment this section

% if animal_ID ==48
%     StimData(:,5500:end) = StimData(:,5500:end) - repmat(mean(StimData (:,5200:5500),2), 1, 13100-5499);
% elseif animal_ID ==50
% 
% elseif animal_ID ==51
% 
%     StimData(:,5100:end) = StimData(:,5100:end) - repmat(mean(StimData (:,5100:5500),2), 1, 13100-5099);
% end



hold on
plot(nanmean(StimData))

%StimData = StimData - repmat(mean(StimData,2), 1, 13100);

% trial by trial reaction time (in fact reward stim interval)
RT = BehData(:,14) - BehData(:,13);
%RT = BehData(:,10) - BehData (:,13); % compute choice reaction time from action-stim interval
%RT(RT < 0) = nan;   %few trials with error negative RTs

BehData(:,7)= BehData(:,10) - BehData (:,13);  % put these RTs in the matrix

ToLargeR = find((BehData(:,3)==-1 & BehData(:,8)==1)  | ...
    (BehData(:,3)==1 & BehData(:,8)==2));


BehData(ToLargeR,16)=1;

ToSmallR = setdiff(1:size(BehData,1),ToLargeR)';

BehData(ToSmallR, 16)=-1;

BehData (:,17) = nan;  % adding aquired reward
BehData (intersect(ToLargeR,find(BehData(:,9)==1)),17) = 2;
BehData (intersect(ToSmallR,find(BehData(:,9)==1)),17) = 1;
BehData (isnan(BehData(:,17)),17)= 0;


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

inverted_Matrix_to_row = smooth(inverted_Matrix_to_row,5);

inverted_Matrix_to_row = inverted_Matrix_to_row';

% this doconvolve was working before but it is too slow recently:
%[inverted_Matrix_to_row, s] = deconvolveCa(inverted_Matrix_to_row,'ar2'); % smoooth

inverted_Matrix_to_row = inverted_Matrix_to_row - (mean(inverted_Matrix_to_row)/1.5);


% inverted_Matrix_to_row = zscore(inverted_Matrix_to_row); %not that good

StimDataSmoothDownSample = reshape(inverted_Matrix_to_row', 262,size(BehData,1))';

t=downsample(t,Timereso);


% kernel regression

[eventTimes, eventValues,windows]=Salvatore_SelectKernelsPhotoM(ModelArrangment,beep_time,stim_time, action_onsetTime, outcome_time);


for fititer=1:1:4
    
    
    [fitKernels, predictedSignals,cvErr] = kernelRegression(inverted_Matrix_to_row, t', eventTimes, eventValues, windows , [0], [0 0]);
    
    Coefiz=ones(length(beep_time),4);
    
    
    % recover observed and estimated PSTHs
    for itrial = 1 : length(beep_time)-1
        
        zeroPadBeep = zeros(1,262);
        
        zeroPadAction = zeros(1,262);
        zeroPadStim = zeros(1,262);
        zeroPadOutcome = zeros(1,262);
        
        
        if ModelArrangment ==8
            
            zeroPadStim (71+4 :71+47+4) =  (fitKernels{2}');
            
            t_beep = floor((BehData(itrial, 13 ) - BehData(itrial, 12))*24); % 11 s is 262 samples, thus 1 s is 24
            zeroPadBeep(71-t_beep -8 : 71 -t_beep+ 31-8) = ( fitKernels{1}');
            
            t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
           zeroPadOutcome(71+t_out: 71+t_out +39) = (fitKernels{4}');
            
            t_act = floor((BehData(itrial, 10 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            t_act(isnan(t_act) )=t_out  -4;
            
            zeroPadAction(71+t_act-20 : 71+t_act-20 +23) = (fitKernels{3}');
           
            [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadBeep',zeroPadStim',zeroPadAction',zeroPadOutcome']);
            
            
            tofill = [1,2,3,4];
            Coefiz(itrial,tofill)= B';
        end
        
        if ModelArrangment ==9
            
            zeroPadStim (71+4 :71+47+4) =  (fitKernels{1}');
            
            t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            
            t_act = floor((BehData(itrial, 10 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            t_act(isnan(t_act) )=t_out  -4;
            
            zeroPadAction(71+t_act-20 : 71+t_act-20 +23) = (fitKernels{2}');
           
            [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadStim',zeroPadAction']);
            
            
            tofill = [2,3];
            Coefiz(itrial,tofill)= B';
        end
        
        if ModelArrangment ==10
            
            
            t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            
            t_act = floor((BehData(itrial, 10 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            t_act(isnan(t_act) )=t_out  -4;
            
            zeroPadAction(71+t_act-20 : 71+t_act-20 +23) = (fitKernels{1}');
            
            [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadAction']);
            
            
            tofill = [3];
            Coefiz(itrial,tofill)= B';
        end
        
        
        if ModelArrangment ==11
            
            
%            zeroPadStim (71+4 :71+47+4) =  (fitKernels{2}');
            zeroPadStim (71+4 :71+55+4) =  (fitKernels{2}'); %for plotting mega raster
            
            t_beep = floor((BehData(itrial, 13 ) - BehData(itrial, 12))*24); % 13.1 s is 262 samples, thus 1 s is 24
            zeroPadBeep(71-t_beep -8 : 71 -t_beep+ 31-8) = ( fitKernels{1}');
            
            t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
             
%            zeroPadOutcome(71+t_out: 71+t_out +39) = (fitKernels{3}');
            
            zeroPadOutcome(71+t_out-8: 71+t_out-8 +87) = (fitKernels{3}'); % for plotting mega raster
            
           [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadBeep',zeroPadStim',zeroPadOutcome']);

            
            tofill = [1,2,4];

            Coefiz(itrial,tofill)= B';
            
            
        end
        
         if ModelArrangment ==15
            
            
            zeroPadStim (71+4 :71+47+4) =  (fitKernels{1}');
             
            t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
             
            zeroPadOutcome(71+t_out: 71+t_out +39) = (fitKernels{2}');
            
          
            [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadStim',zeroPadOutcome']);

            tofill = [2,4];

            Coefiz(itrial,tofill)= B';
            
            
         end
        
        if ModelArrangment ==12
            
            
            t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            zeroPadOutcome(71+t_out: 71+t_out +39) = (fitKernels{2}');
            
           
            t_act = floor((BehData(itrial, 10 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            t_act(isnan(t_act) )=t_out  -4;
            
            zeroPadAction(71+t_act-20 : 71+t_act-20 +23) = (fitKernels{1}');
            
            [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadAction',zeroPadOutcome']);
            
            
            tofill = [3,4];
            Coefiz(itrial,tofill)= B';
        end
        
        
        if ModelArrangment ==13
            
            zeroPadStim (71+4 :71+47+4) =  (fitKernels{1}');
            
            [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadStim']);
            
            
            tofill = [2];
            Coefiz(itrial,tofill)= B';
        end
        
        
        if ModelArrangment == 14
            
            
            t_out = floor((BehData(itrial, 14 ) - BehData(itrial, 13))*24); % 11 s is 328 samples, thus 1 s is 30
            zeroPadOutcome(71+t_out: 71+t_out +39) = (fitKernels{1}');
            
            
            [B] = regress(StimDataSmoothDownSample(itrial,:)',[zeroPadOutcome']);
            
            
            tofill = [4];
            Coefiz(itrial,tofill)= B';
        end
         
        
        for i=1:length(eventTimes)
            
            % this is for shuffling coef across trials
            
            % Coefiz(:,i) = Coefiz(randperm(length(Coefiz(:,i))),1);
            
            eventValues{i} = Coefiz(:,tofill(i))';
            
        end
        
        
    end
    
    % compute EV
    EV(fititer) = 1-mean(mean((predictedSignals-inverted_Matrix_to_row).^2))/mean(mean(inverted_Matrix_to_row.^2));
    
end

%%
          % this part is fitting the correct RL model on choices, using saved
          %  parameters.
            ExpID='23';
            
            CrossVal= LoadSavedRLFits (animal_name, ExpID);
            
            
            FullModelParamZ(1,:)=CrossVal.Results(1,1).TopModelParams;
            FullModelParamZ(2,:)=CrossVal.Results(1,2).TopModelParams;
            
            FullModelParam = nanmean(FullModelParamZ);
            iter=10;
            
            TempBehData=BehData;
            
            TempBehData(:,10)=BehData(:,9);
            
            BehFitData.data = TempBehData(:,1:15);
            
            [data_models, action, ~, ~, ~ ] = RunPOMDP_GS_NLL_withLapse(BehFitData,FullModelParam,iter);
            
            
            ChosenVal = [];
            
            ChosenVal(data_models(:,3)==-1) =    data_models(data_models(:,3)==-1,17) ;
            ChosenVal(data_models(:,3)==1) =     data_models(data_models(:,3)==1,18) ;
            
      
            RPE = data_models(:,21)' ;
            
            
            % regress coeffients on chosen value and RPEs
              iChosenValRange = linspace (min(ChosenVal)+0.1, max(ChosenVal)-0.1,9);
     
            iChosenValInt = iChosenValRange(2) - iChosenValRange(1);
             
            iChosenValCount = 1;
            for iChosenVal = iChosenValRange (1:end-1)
                
                indi = intersect(find(ChosenVal > iChosenVal ),find(ChosenVal < (iChosenVal + iChosenValInt))) ;
                
                DAvsChosenValBIN (iChosenValCount)    = nanmean(Coefiz(indi,2));
  
                iChosenValCount = iChosenValCount + 1;
            end
            
             stat_ChosenVal=regstats(DAvsChosenValBIN,iChosenValRange (1:end-1))
            
              iRPERange = linspace (min(RPE)+0.1, max(RPE)-0.1,9);
           

           iRPERangeInt = iRPERange(2) - iRPERange(1);
             
            iRPECount = 1;
            for iRPE = iRPERange (1:end-1)
                
                indi = intersect(find(RPE > iRPE ),find(RPE < (iRPE + iRPERangeInt))) ;
                
                DAvsRPEBIN (iRPECount)     = nanmean(Coefiz(indi,4));
  
                iRPECount = iRPECount + 1;
            end
            
       stat_RPE=regstats(DAvsRPEBIN,iRPERange (1:end-1))
%%
          
            
figure

plot(t,smooth(inverted_Matrix_to_row + randn(1,length(inverted_Matrix_to_row)),20));
hold on

plot(t,smooth(predictedSignals + randn(1,length(predictedSignals)),20));

correct = BehData(:,9);


for i = outcome_time(correct==1)
    plot(i,1,'g*')
    
end

for i = outcome_time(correct==0)
    plot(i,1,'r*')
    
end

for i = stim_time
    plot(i,1,'k*')
    
end


SignalStim = [];
SignalAction = [];
SignalOutcome = [];

EstimatedSignalStim = [];
EstimatedSignalAction = [];
EstimatedSignalOutcome = [];

%predictedSignals = zscore(predictedSignals);

for iter = 1 : length(stim_time) -1  % we ingore the last trial
    
    
    [diff_recover i_recover] = min (abs(t - stim_time(iter)));
    
    EstimatedSignalStim (iter,:) =predictedSignals(i_recover-10 : i_recover+200);
    SignalStim (iter,:) =inverted_Matrix_to_row(i_recover-10 : i_recover+200)';
    
    [diff_recover i_recover] = min (abs(t - action_onsetTime(iter)));
    
    EstimatedSignalAction (iter,:) =predictedSignals(i_recover-20 : i_recover+190);
    SignalAction (iter,:) =inverted_Matrix_to_row(i_recover-20 : i_recover+190)';
    
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

%tempIndex = find(BehDataAbs(:,2)==0.25);  % abs stimulus to visualise

%BehDataAbs = BehDataAbs(tempIndex,:);

Stim2Visul = 0.25;

% to align and plot RT on the plot:
t_act = randi([-1 2],length(Stimz),1) + 12 + floor((BehData(:, 10 ) - BehData(:, 13))*24); % 11 s is 328 samples, thus 1 s is 30

t_out = 10 + floor((BehData(:, 14 ) - BehData(:, 13))*24);
% t_act(isnan(t_act) )=t_out  -4;

BehDataAbs (:,18) = t_out;
BehData2sort = BehDataAbs;

BehData2sort(BehData2sort(:,2)==0,2) =0.01; % to separate zero contraat rewarded and unrewarded

%[BehDatasorted j]= sortrows(BehData2sort,[17,2,7]); % final sorting (Correct/Error, abs contrast and RTs)
[BehDatasorted j]= sortrows(BehData2sort,[17,2,18]); % final sorting (Correct/Error, abs contrast and reward time)



figure

subplot(1,2,1)
imagesc(smooth2a(EstimatedSignalStim(j,5:100),0,5),[-1 5])
title('estimated stim')

subplot(1,2,2)
imagesc(smooth2a(SignalStim(j,8:100),0,5),[-1 5])
title('real stim')
colormap('bluewhitered')

hold on
trace = 1;

for ievent=t_act(j)'
    
    plot(ievent,trace,'k.')
    
    trace  = trace + 1;
    
end

hold on
trace = 1;

for ievent=t_out(j)'
    
    plot(ievent,trace,'k.')
    
    trace  = trace + 1;
    
end


figure

subplot(1,2,1)
imagesc(smooth2a(EstimatedSignalStim(j(BehDatasorted(:,2)==Stim2Visul),2:80),0,5),[-1.1 5.5])
title('estimated stim')
hold on

subplot(1,2,2)
imagesc(smooth2a(SignalStim(j(BehDatasorted(:,2)==Stim2Visul),5:83),0,5),[-1 5])
title('real stim')
colormap('bluewhitered')

hold on

for subplotz = [1 2]
   subplot(1,2,subplotz) 
trace = 1;

for ievent=t_act(j(BehDatasorted(:,2)==Stim2Visul))'
    
    plot(ievent,trace,'k.','markersize',3)
    
    trace  = trace + 1;
    
end

trace = 1;

for ievent=t_out(j(BehDatasorted(:,2)==Stim2Visul))'
    azs= j(BehDatasorted(:,2)==Stim2Visul);
    if BehData(azs(trace),17) == 0
        col = [0.7 0 0.2 ];
    elseif BehData(azs(trace),17) == 1
        col = [0 1 0 ];
    elseif BehData(azs(trace),17) == 2
        col = [0 0.7 0];
    end
    plot(ievent,trace,'.','color',col,'markersize',2)
    
    trace  = trace + 1;
    
end
end

figure
subplot(1,2,1)
imagesc(EstimatedSignalOutcome(j,2:50),[0 5])
title('estimated outcome')

subplot(1,2,2)
imagesc(SignalOutcome(j,2:50),[0 5])
title('real outcome')
colormap('bluewhitered')


%%

for i=abzStim
    
    tempIndex = find(abs(Stimz)==i);
    
    [i j]=ismember(i, abzStim);
    
    % correct large
    resp_CorrStimLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==1)),2));
    resp_ErrStimLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0),find(BehData(:,16)==1)),2));
    
    resp_CorrStimSmall(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==-1)),2));
    
    resp_CorrActionLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==1)),3));
    resp_ErrActionLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0),find(BehData(:,16)==1)),3));
    
    resp_CorrActionSmall(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==-1)),3));
    
    resp_CorrOutcomeLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==1)),4));
    resp_ErrOutcomeLarge(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==0),find(BehData(:,16)==1)),4));
    
    resp_CorrOutcomeSmall(j) = nanmedian(Coefiz(mintersect(tempIndex,find(    correct==1),find(BehData(:,16)==-1)),4));
    
    
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

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).Kernels = fitKernels;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrStimLarge = resp_CorrStimLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_ErrStimLarge = resp_ErrStimLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrStimSmall=resp_CorrStimSmall;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrActionLarge=resp_CorrActionLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_ErrActionLarge=resp_ErrActionLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrActionSmall=resp_CorrActionSmall;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrOutcomeLarge=resp_CorrOutcomeLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_ErrOutcomeLarge=resp_ErrOutcomeLarge;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).resp_CorrOutcomeSmall=resp_CorrOutcomeSmall;


%%
StimzAbs = abs(BehData(1:end-1,2));

StimsAllowed = unique(StimzAbs');

%
figure; hold on
for i=unique(StimzAbs')
    
    tempIndex = find(StimzAbs==i);
    
    [i j]=ismember(i, StimsAllowed);
    
    
    PredictedPopStimAlign(j,:)=nanmean(EstimatedSignalStim(tempIndex,:));
    PredictedPopActionAlign(j,:)=nanmean(EstimatedSignalAction(tempIndex,:));
    
    hold on
    
    plot(PredictedPopStimAlign(j,:),'color',colorGray(j,:))
    
    
end
%
figure;


BehData(end,:)=[];
IndexLargeRew = mintersect(find( BehData(:,9)==1),find(BehData(:,16)==1));

IndexSmallRew = mintersect(find(   BehData(:,9)==1),find(BehData(:,16)==-1));

IndexNoRew = find(BehData(:,9)==0);



figure; hold on
plot(nanmean(EstimatedSignalOutcome(IndexLargeRew,:)));

plot(nanmean(EstimatedSignalOutcome(IndexSmallRew,:)));

plot(nanmean(EstimatedSignalOutcome(IndexNoRew,:)));


BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopStimAlign=PredictedPopStimAlign;
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopActionAlign=PredictedPopActionAlign;

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopRewAlignLarge=nanmean(EstimatedSignalOutcome(IndexLargeRew,:));
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopRewAlignSmall=nanmean(EstimatedSignalOutcome(IndexSmallRew,:));
BehPhotoM(animal_ID).KernelSummary(ModelArrangment).PredictedPopRewAlignNoRew=nanmean(EstimatedSignalOutcome(IndexNoRew,:));

BehPhotoM(animal_ID).KernelSummary(ModelArrangment).EV=EV;


