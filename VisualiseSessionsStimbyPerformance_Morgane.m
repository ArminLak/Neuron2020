% See how performance and response size changes with stimulus contrast (and later: across sessions)
% Morgane Moss Oct 2018

clear all
% close all

animal_name = 'ALK070'

load('BehPhotoM_Exp7_VTA')  % load beh data database 


[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);  % read animals' ID 

%%

SessionList = [1:9]; %for now keep this one number 


StimPerformance = nan(max(SessionList), 5); % 1 to 5 are 0, 0.12, 0.25, 0.5, 1.0 contrast levels. 
StimResponse = nan(max(SessionList), 5); % same as above but this will store avg response as a Z-score rel. to that session 

StimzAbs = [1 0.5 0.25 0.12 0];

event_time = 3; % this is the time in the summary matrix where the event took place

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

preAlignStim = (event_time-0.4)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignStim = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;

% ------------- create performance matrix --------------------------------
% columns are stim 0 , 0.12 , 0.25, 0.5, 1 (blanks are nan)
% rows are sessions in order 

for isession = SessionList 
    
    TrialTimingData = BehPhotoM(animal_ID).Session(isession).TrialTimingData;  % load beh data 
    TrialTimingDataCor = TrialTimingData(TrialTimingData(:,9)==1, :);
    NeuronStim = BehPhotoM(animal_ID).Session(isession).NeuronStim; % load photom data 

    for stimcount = 1:length(StimzAbs)
        istim = StimzAbs(stimcount)

        StimPerformance(isession, stimcount) = length(find(abs(TrialTimingDataCor(:,2))==istim))...
            / length(find(abs(TrialTimingData(:,2))==istim)); % get P(correct)
        
        StimResponse(isession, stimcount) = nanmean(nanmean(NeuronStim(abs(TrialTimingDataCor(:,2))==istim, postAlignStim)))...
            - nanmean(nanmean(NeuronStim(abs(TrialTimingDataCor(:,2))==istim, preAlignStim))); %difference between before and after; indicator of relative signal change 

    end

end

StimPerformance(StimPerformance == 0) = nan;



%%
% ---------- scatter --------- 

figure; hold on 
c = 1;

for isession = SessionList
    
    for stimcount = 1:length(StimzAbs)
        
        if StimResponse(isession, stimcount) > 0
            plot(StimPerformance(isession, stimcount), StimzAbs(stimcount), 'o', 'color', [1-StimResponse(isession, stimcount) 1-StimResponse(isession, stimcount) 1-StimResponse(isession, stimcount)],...
            'markerFaceColor', [1-StimResponse(isession, stimcount) 1-StimResponse(isession, stimcount) 1-StimResponse(isession, stimcount)])
        elseif StimResponse(isession, stimcount) < 0
            plot(StimPerformance(isession, stimcount), StimzAbs(stimcount), 'o', 'color', [0 0 1+StimResponse(isession, stimcount)], ...
                'markerFacecolor', [0 0 1+StimResponse(isession, stimcount)])
        end
        
    end
    
    
end

xlabel('P(correct)')
yticks([0 0.12 0.25 0.5 1.0])
ylabel('Stim contrast')


%%
% ------ heat map --------





