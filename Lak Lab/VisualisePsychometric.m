% function [performance] = VisualisePsychometric(computer_name, subject, date, session)

clear all


computer_name = 'WIN-AL003';
subject = 'exampleSubject';
date = '2019-10-10';
session = '2';


path2beh = ['\\', computer_name, '\b_server\Data\subjects\', subject, '\', date, '\', session];

addpath(genpath(path2beh));

load([date, '_', session, '_', subject, '_Block.mat']);
load([date, '_', session, '_', subject, '_parameters.mat']);

% trials with large reward on left 
largeLeftTrials = block.events.highRewardSideValues(1:length(block.events.feedbackTimes))  == -1;
largeRightTrials = block.events.highRewardSideValues(1:length(block.events.feedbackTimes)) == 1;
c = 1;

contrasts = block.events.contrastRightValues - block.events.contrastLeftValues;
Stimz = sort(unique(contrasts), 'ascend')';

for iStim = Stimz'
    
    LeftBlockHit(c)         = sum(largeLeftTrials & contrasts(1:length(block.events.feedbackTimes))==iStim & block.events.responseValues==1);
    LeftBlockTotal(c)       = sum(largeLeftTrials & contrasts(1:length(block.events.feedbackTimes))==iStim);
    
    RightBlockHit(c)        = sum(largeRightTrials & contrasts(1:length(block.events.feedbackTimes))==iStim & block.events.responseValues==1);
    RightBlockTotal(c)      = sum(largeRightTrials & contrasts(1:length(block.events.feedbackTimes))==iStim);

    c = c + 1;
end

% remove contrasts with no trials 

toRemove = find(RightBlockTotal == 0);
RightBlockTotal(toRemove) = [];
RightBlockHit(toRemove) = [];

toRemove = find(LeftBlockTotal == 0);
LeftBlockTotal(toRemove) = [];
LeftBlockHit(toRemove) = [];

StimzAllowed = Stimz'

figure; 
plot(StimzAllowed, LeftBlockHit./LeftBlockTotal, 'color', [0.5 0.2 0.1], 'LineWidth', 2)
hold on; 
plot(StimzAllowed, RightBlockHit./RightBlockTotal, 'color', [1 0.6 0.2], 'LineWidth', 2)
set(gca, 'YTick', [0 0.5 1], 'XTick', Stimz)
ylabel('P(R)')
xticklabel(abs(Stimz))


% end
