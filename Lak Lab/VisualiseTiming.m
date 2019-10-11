% Morgane October 2019
% Visualise histograms of time delays between trial events 

close all
clear all


computer_name = 'WIN-AL003';
subject = 'exampleSubject';
date = '2019-10-10';
session = '2';


path2beh = ['\\', computer_name, '\b_server\Data\subjects\', subject, '\', date, '\', session];

addpath(genpath(path2beh));

load([date, '_', session, '_', subject, '_Block.mat']);
load([date, '_', session, '_', subject, '_parameters.mat']);

TrialStartToStimulus = block.events.stimulusOnTimes - block.events.newTrialTimes;
StimulusToGoCue = block.events.interactiveOnTimes - block.events.stimulusOnTimes;
ResponseToReward = block.events.feedbackTimes - block.events.responseMadeTimes;


figure; 
plots(1) = subplot(1,3,1)
histogram(TrialStartToStimulus)
title('Time from Trial Start to Stimulus onset')
ylabel('Frequency')
xlabel('Time(s)')

plots(2) = subplot(1,3,2)
histogram(StimulusToGoCue)
title('Time from stimulus onset to go cue')
xlabel('Time(s)')


plots(3) = subplot(1,3,3)
histogram(ResponseToReward)
title('Time from response made to reward delivered')
xlabel('Time(s)')

set(gcf, 'Position', [100 300 1000 400])

