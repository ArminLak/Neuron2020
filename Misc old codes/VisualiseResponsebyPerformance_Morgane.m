% See how response to stimulus changes with overall performance across sessions
% Morgane Moss Oct 2018

clear all
% close all

animal_name = 'ALK068'

load('BehPhotoM_Exp7_VTA')  % load beh data database 


[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);  % read animals' ID 

%%

SessionList = 1:length(BehPhotoM(animal_ID).Session); %for now keep this one number 


StimPerformance = nan(max(SessionList), 5); % 1 to 5 are 0, 0.12, 0.25, 0.5, 1.0 contrast levels. 
StimResponse = nan(max(SessionList), 5); % same as above but this will store avg response as a Z-score rel. to that session 

StimzAbs = [1 0.5 0.25 0.12 0];

event_time = 3; % this is the time in the summary matrix where the event took place

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

preAlignStim = (event_time-0.4)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignStim = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;

%%
% ------------- create performance matrix --------------------------------
% columns are stim 0 , 0.12 , 0.25, 0.5, 1 (blanks are nan)
% rows are sessions in order 

for isession = SessionList 
    
    TrialTimingData = BehPhotoM(animal_ID).Session(isession).TrialTimingData;  % load beh data 
    TrialTimingDataCor = TrialTimingData(TrialTimingData(:,9)==1, :);
    NeuronStim = BehPhotoM(animal_ID).Session(isession).NeuronStimR; % load photom data 

    for stimcount = 1:length(StimzAbs)
        istim = StimzAbs(stimcount);

        StimPerformance(isession, stimcount) = length(find(abs(TrialTimingDataCor(:,2))==istim))...
            / length(find(abs(TrialTimingData(:,2))==istim)); % get P(correct)
        
        StimResponse(isession, stimcount) = nanmean(nanmean(NeuronStim(intersect(find(abs(TrialTimingData(:,2))==istim), find(TrialTimingData(:,9)==1)), postAlignStim)))...
            - nanmean(nanmean(NeuronStim(intersect(find(abs(TrialTimingData(:,2))==istim), find(TrialTimingData(:,9)==1)), preAlignStim))); %difference between before and after; indicator of relative signal change 

    end

end

StimPerformance(StimPerformance == 0) = nan;


% ---------- scatter --------- 

figure; hold on 

for isession = SessionList
    
    for stimcount = 1:length(StimzAbs)
        
        c = (isession-1)/max(SessionList);
        
        % first session = black, last session = lightest 
        
        plot(StimResponse(isession, stimcount), StimPerformance(isession, stimcount), 'o', 'color', [c c c], 'MarkerFaceColor', [c c c], 'MarkerSize', 12)
        
    end
    
    
end

ylabel('P(correct)')
xlabel('Response to stimulus')

%% compute performance-response correlation for naive, middle and adanced


    SessNum = ceil(length(SessionList)/3);

    BinSessions= 1:SessNum:length(SessionList);
    
    if strcmp(animal_name,'ALK071')
        SessNum=4; BinSessions=[1 6 9]; end

     c= 1; 
    for iBin = BinSessions(1:3)

        if c < 3
        Response=StimResponse(iBin:SessNum+iBin-1,:);
        Performance=StimPerformance(iBin:SessNum+iBin-1,:);

        else
Response=StimResponse(iBin:end,:);
        Performance=StimPerformance(iBin:end,:);
            
        end            
        ResponsenoNan= Response(find(~isnan(Response)));
        PerformanceNoNan= Performance(find(~isnan(Performance)));
        
        
        CorEstimates (c)=corr(ResponsenoNan(:),PerformanceNoNan(:))
c=c+1;
    end

    figure
    bar(CorEstimates)

   



