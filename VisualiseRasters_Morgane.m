clear all
close all

%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 51
load('BehPhotoM_Exp23')

RTLimit = 10; % in s, excluding trials with RT longer than this



%%


sessionz = 1:length(BehPhotoM(animal_ID).Session);


for iSession = sessionz
    
    TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    
    TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
    
    TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStim;
    
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;
    

    % here you can use imagesc to make colourful plots
    % for each session lets make the following subplots
    
    % aling to stimulus, sorted based on RT
    
    % subplots for different levels of stimuli (average L and R)
    
    % subplots (outcome aligned) for large reward, small reward and no
    % reward
    
    % so you will write imagesc
    
    % then call colormap(bluewhitered)
end