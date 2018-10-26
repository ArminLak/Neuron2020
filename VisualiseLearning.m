% Visualise avg responses to stimulus and reward in specified sessions

% Armin Oct 2018


% to do instead of z-scoring normalise to max response

clear all
close all




% animal_name= ['ALK068'];
animal_list= [{'ALK068', 'ALK070', 'ALK071'}];


animal_list= [{'ALK078', 'MMM001', 'MMM002'}];


%animal_list= [{'ALK074', 'ALK075'}];

%animal_list= [{'ALK074'}];

animal_ID_list_VTA = [48  50  51];

animal_ID_list_NAc =[56, 57,59];

animal_ID_list_DMS =[53, 55];



% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = -1 % s this should be -1 or less
stop = 2    % s
event_time = 3; % this is the time in the summary matrix where the event took place

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

start2stop = (event_time+start)*sample_rate/downsampleScale:(event_time+stop)*sample_rate/downsampleScale; %window of interest

% event is at 3 seconds i.e. point 3600

preAlignStim = (event_time-0.4)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignStim = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;
preAlignRwd = (event_time-0.2)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignRwd = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;


% plot colours
colorGray4 = [0.8 0.8 0.8 %lightest
    0.6 0.6 0.6
    0.4 0.4 0.4
    0.2 0.2 0.2
    0 0 0]; % black
colorGray3 = [0.8 0.8 0.8
    0.4 0.4 0.4
    0 0 0];
colorGray2 = [0.7 0.7 0.7
    0 0 0];
colorGreen = [0 1 0
    0 0.8 0
    0 0.6 0
    0  0.3 0];



NeuronRewardErr=nan(length(animal_list),20,5,13100);
NeuronRewardCor=nan(length(animal_list),20,5,13100);

NeuronStimCor =nan(length(animal_list),20,5,13100);
NeuronStimErr =nan(length(animal_list),20,5,13100);

NeuronActionCor =nan(length(animal_list),20,5,13100);
NeuronActionErr =nan(length(animal_list),20,5,13100);


for animalcount = 1:length(animal_list)
    animal_name = animal_list(animalcount);
    [animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);
    
    if ismember(animal_ID,animal_ID_list_VTA)
        load('BehPhotoM_Exp7_VTA')                                   % load beh data databse
        
    elseif ismember(animal_ID,animal_ID_list_NAc)
        load('BehPhotoM_Exp7_NAc')                                   % load beh data databse
        
    elseif ismember(animal_ID,animal_ID_list_DMS)
        load('BehPhotoM_Exp7_DMS')                                   % load beh data databse
        
    end
    
    
    SessionList = 1:length(BehPhotoM(animal_ID).Session);
    
    for iSession = SessionList
        
        
        TrialTimingData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        
        
        NeuronStim = BehPhotoM(animal_ID).Session(iSession).NeuronStim;
        
        
        NeuronStimTurned=NeuronStim';
        NeuronStimVector=NeuronStimTurned(:);
        NeuronStimVector = zscore(NeuronStimVector);
        NeuronStimNorm = reshape(NeuronStimVector,size(NeuronStim,2),size(NeuronStim,1));
        NeuronStimNorm = NeuronStimNorm';
        
        NeuronStim = NeuronStimNorm;
        
        
        NeuronAction = BehPhotoM(animal_ID).Session(iSession).NeuronAction;
        
        
        NeuronActionTurned=NeuronAction';
        NeuronActionVector=NeuronActionTurned(:);
        NeuronActionVector = zscore(NeuronActionVector);
        NeuronActionNorm = reshape(NeuronActionVector,size(NeuronAction,2),size(NeuronAction,1));
        NeuronActionNorm = NeuronActionNorm';
        
        
        NeuronAction = NeuronActionNorm;
        
        NeuronReward = BehPhotoM(animal_ID).Session(iSession).NeuronReward;
        
        NeuronRewardTurned=NeuronReward';
        NeuronRewardVector=NeuronRewardTurned(:);
        NeuronRewardVector = zscore(NeuronRewardVector);
        NeuronRewardNorm = reshape(NeuronRewardVector,size(NeuronReward,2),size(NeuronReward,1));
        NeuronRewardNorm = NeuronRewardNorm';
        
        NeuronReward = NeuronRewardNorm;
        
        
        
        c = 5;
        for istim = fliplr(unique(abs(TrialTimingData(:,2))'))
            
            
            
            NeuronStimCor(animalcount,iSession,c,:) = mean(NeuronStim(TrialTimingData(:,9)==1 & abs(TrialTimingData(:,2))==istim, :));
            NeuronActionCor(animalcount,iSession,c,:) = mean(NeuronAction(TrialTimingData(:,9)==1 & abs(TrialTimingData(:,2))==istim, :));
            
            NeuronRewardCor(animalcount,iSession,c,:)  = mean(NeuronReward(TrialTimingData(:,9)==1 & abs(TrialTimingData(:,2))==istim, :));
            
            NeuronStimErr(animalcount,iSession,c,:) = mean(NeuronStim(TrialTimingData(:,9)==0 & abs(TrialTimingData(:,2))==istim, :));
            NeuronActionErr(animalcount,iSession,c,:) = mean(NeuronAction(TrialTimingData(:,9)==0 & abs(TrialTimingData(:,2))==istim, :));
            
            
            NeuronRewardErr(animalcount,iSession,c,:)  = mean(NeuronReward(TrialTimingData(:,9)==0 & abs(TrialTimingData(:,2))==istim, :));
            
            
            
            c=c-1;
            
            
        end
        
    end
    
end


%%
figure

if length(animal_list)==1
    for i=1:9
        
        subplot(9,3,3*i-2); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronStimCor(:,i,istim,:)),2),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        xlim([3700 4500])
        ylim([-1 2])
        
        subplot(9,3,3*i-1); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronActionCor(:,i,istim,:)),2),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        xlim([3500 4500])
        ylim([-1 2])
        
        subplot(9,3,3*i); hold on
        
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronRewardCor(:,i,istim,:)),2),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        ylim([-1 2])
        
        xlim([3500 4500])
    end
    
else
    
    for i=1:9
        
        subplot(9,3,3*i-2); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronStimCor(:,i,istim,:))),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        xlim([3700 4500])
        ylim([-1 2])
        
        subplot(9,3,3*i-1); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronActionCor(:,i,istim,:))),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        xlim([3500 4500])
        ylim([-1 2])
        
        subplot(9,3,3*i); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronRewardCor(:,i,istim,:))),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        ylim([-1 2])
        
        xlim([3500 4500])
    end
end


