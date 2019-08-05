% script generating psychometric curves from mice 2AFC choice

close all
clear

animal_ID = 48

%animals = [48 50 51]

%load('BehPhotoM_Exp23_VTA.mat')

animals = [53, 62, 63, 71,72] 
load('BehPhotoM_Exp23_DMS')

%animals = [56 57 59 66]
%load('BehPhotoM_Exp23_NAc')


plotCount = 1;
figure;

for animal_ID = animals 

BehData = [];
BeepData = [];
StimData = [];
ActionData = [];
RewardData = [];

   sessionz = 1:length(BehPhotoM(animal_ID).Session);

    for iSession = sessionz

        TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        BehData = [BehData; TempBehData];

    end

color=[127/255 75/255 31/255
    247/255 143/255 30/255];

BehData(BehData(:,3)==-1,3)=0; % change -1 and 1 to 0 and 1 indicating L and R choice

RTLimit = 6; % in s, excluding trials with RT longer than this
RT = BehData(:,10) - BehData(:,13);

toRemove = find ( RT > RTLimit); %remove trials with slow reaction time
BehData(toRemove,:) = [];


testedStims = unique(BehData(:, 2));  % find tested stimuli

leftIndexes = find(BehData(:,8)==1); % find index for large and right blox
rightIndexes = find(BehData(:,8)==2);

LeftData = BehData(leftIndexes, :);
RightData = BehData(rightIndexes, :);

% ------ 1) Left blocks ----------------------------------------------


stimCounter=1;                     % counter for istim loop
sem=nan(1, 9);                     % initialising sem

for istim = testedStims'
    
   Nhit(stimCounter)=sum(LeftData(LeftData(:,2)==istim,3));
   Ncorrect(stimCounter)=sum(LeftData(LeftData(:,2)==istim,10));
   
   stimCounter = stimCounter+1;
end

% figure;
ylim([0 1])

Noutof = histc(LeftData(:,2), testedStims)';   % very cool coding here! 
pR = Nhit ./ Noutof;                     % compute mean p_R per stim
[phat,pci] = binofit(Nhit, Noutof)% binomial dist C.I.
pci = (pci(:,2)-pci(:,1))./2
psuccess=Ncorrect ./ Noutof;

[M,V]=binostat(Noutof, psuccess)
std = sqrt(V)./100

subplot(1, length(animals), plotCount)
[paramsValues] = Fit_psych_fun_Armin(testedStims', Nhit, Noutof,[1 1 1 1],1,color(1,:))
ylim([0 1]), xlabel('Contrast'), ylabel('Rightward Choice (%)')
xlim([min(testedStims) max(testedStims)])
yticks([0 0.5 1])
xticks([-0.75 0 0.75])
set(gca, 'YTickLabel', [0 50 100])
hold on

%confidence intervals:
errorbar(testedStims', pR, pci, 'o', 'MarkerSize', 1,...
    'MarkerFaceColor', color(2,:), 'MarkerEdgeColor', color(1,:),...
    'Color', color(1,:), 'LineWidth', 2)
hold on
plot(testedStims,pR, 'o','LineWidth',2,...
                       'Color',color(1,:),...
                       'MarkerSize',5,...
                       'MarkerFaceColor', color(1,:),'MarkerEdgeColor', color(1,:) );
legend('Psych function', 'Unfitted', 'Location', 'northeastoutside'), 
if animal_ID == 48
    title('Mouse 68')
elseif animal_ID == 50
    title('Mouse 70')
elseif animal_ID == 51
    title('Mouse 71')
end


% ------ 2) Right blocks ----------------------------------------------


stimCounter=1;                     % counter for istim loop
sem=nan(1, 9);                     % initialising sem

for istim = testedStims'
    
   Nhit(stimCounter)=sum(RightData(RightData(:,2)==istim,3));
   Ncorrect(stimCounter)=sum(RightData(RightData(:,2)==istim,10));
   
   stimCounter = stimCounter+1;
end

hold on
ylim([0 1])

Noutof = histc(RightData(:,2), testedStims)';   % very cool coding here! 
pR = Nhit ./ Noutof;                     % compute mean p_R per stim
[phat,pci] = binofit(Nhit, Noutof)% binomial dist C.I.
pci = (pci(:,2)-pci(:,1))./2
psuccess=Ncorrect ./ Noutof;

[M,V]=binostat(Noutof, psuccess)
std = sqrt(V)./100


[paramsValues] = Fit_psych_fun_Armin(testedStims', Nhit, Noutof,[1 1 1 1],1,color(2,:))
% ylim([0 1]), xlabel('Contrast'), ylabel('% Rightward Choice')
% xlim([min(testedStims) max(testedStims)])
% yticks([0 0.5 1])
% xticks([-0.75 0 0.75])
% set(gca, 'YTickLabel', [0 50 100])
hold on

%confidence intervals:
errorbar(testedStims', pR, pci, 'o', 'MarkerSize', 1,...
    'MarkerFaceColor', color(2,:), 'MarkerEdgeColor', color(2,:),...
    'Color', color(2,:), 'LineWidth', 2)
hold on
plot(testedStims,pR, 'o','LineWidth',2,...
                       'Color',color(2,:),...
                       'MarkerSize',5,...
                       'MarkerFaceColor', color(2,:),'MarkerEdgeColor', color(2,:) );
% legend('Psych function', 'Unfitted', 'Location', 'northeastoutside'), 
% title('95% Confidence Intervals')

plotCount = plotCount + 1

end 