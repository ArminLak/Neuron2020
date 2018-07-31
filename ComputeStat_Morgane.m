% this code does all the t-tests and linear regressions for morgane's masters dissertation 
% inspired by 'VisualiseMultipleSessions'

%% section 1: 2-sample t tests for each animal

clear all
close all

%VTA : [48, 50,51]  coresponding to ALK068, 70 and 71
% NAc : [57] coresponding to MMM001
% DMS : [53] coresponding to ALK074


% select animal
animal_ID = 51

% select database
load('BehPhotoM_Exp23')
%load('BehPhotoM_Exp23_NAc')
%load('BehPhotoM_Exp23_DMS')

% define implant
Implant = 'Un' 

RTLimit = 10; % in s, excluding trials with RT longer than this

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








