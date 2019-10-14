clear all


computer_name = 'WIN-AL003';
subject = 'exampleSubject';
date = '2019-10-14';
session = '2';


path2beh = ['\\', computer_name, '\b_server\Data\subjects\', subject, '\', date, '\', session];

addpath(genpath(path2beh));

load([date, '_', session, '_', subject, '_Block.mat']);
load([date, '_', session, '_', subject, '_parameters.mat']);


nTrials = length(block.events.feedbackTimes);
% trials with large reward on left 
LeftBlockTrials = block.events.highRewardSideValues(1:nTrials)  == -1;
RightBlockTrials = block.events.highRewardSideValues(1:nTrials) == 1;
c = 1;

contrasts = block.events.contrastRightValues(1:nTrials) - block.events.contrastLeftValues(1:nTrials);
block.events.contrastLeftValues = -block.events.contrastLeftValues; 
RTs = block.events.responseTimes(1:nTrials) - block.events.stimulusOnTimes(1:nTrials);
Stimz = sort(unique(contrasts), 'ascend')';

for iStim = Stimz'
    
    % mean 
    LeftBlockMeanRTs(c)   = mean(RTs(LeftBlockTrials & contrasts==iStim));
    LeftBlockErrors(c)    =  std(RTs(LeftBlockTrials & contrasts==iStim)) ./ sqrt(length(RTs(LeftBlockTrials & contrasts==iStim)));
 
    RightBlockMeanRTs(c)  = mean(RTs(RightBlockTrials & contrasts==iStim));
    RightBlockErrors(c)   = std(RTs(RightBlockTrials & contrasts==iStim)) ./ sqrt(length(RTs(RightBlockTrials & contrasts==iStim)));

    c = c + 1;
end


figure; 

errorbar(Stimz, LeftBlockMeanRTs, LeftBlockErrors, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.5 0.2 0.1])
hold on;
errorbar(Stimz, RightBlockMeanRTs, RightBlockErrors, 'LineStyle', '-', 'LineWidth', 2, 'Color', [1 0.6 0.2])

legend({'Large on Left', 'Large on Right'});
xlabel('Stimulus contrast')
ylabel('Reaction time (s)')


