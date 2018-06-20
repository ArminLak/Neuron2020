clear all
close all

%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 51
load('BehPhotoM_Exp23')

RTLimit = 10; % in s, excluding trials with RT longer than this


% ------------------ start stop times for task events ------------------

sstart = -0.2; %stimulus 
sstop = 0.8;

astart = -0.5; %action
astop = 0.5;

rstart = -0.2; %reward
rstop = 0.8;


%%


sessionz = 1:length(BehPhotoM(animal_ID).Session);


for iSession = sessionz
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStim;
    
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;
    
    TrialTimingData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    TempStimz = unique(TrialTimingData(:,2));
    

    % here you can use imagesc to make colourful plots
    % for each session lets make the following subplots
    
    f = figure('Position', [300 200 800 900]); hold on
     
    
    %------------------------------- plot psychometric curve------------------
   %% 
    subplot(4, 7, 1);
    
    c = 1;
    for istim = TempStimz'
    
        performance(c) = nanmean (TrialTimingData(TrialTimingData(:,2)==istim,3));
%         RT(c) = nanmean (ReactionTime(TrialTimingData(:,2)==istim));
    
        c=c+1; 
    end


    plot(TempStimz, performance,'k')
%%
    
    % ---------- raster for all trials in session in order of RT ---------


    
    % aling to stimulus, sorted based on RT
    
    % subplots for different levels of stimuli (average L and R)
    
    % subplots (outcome aligned) for large reward, small reward and no
    % reward
    
    % so you will write imagesc
    
    % then call colormap(bluewhitered)
end