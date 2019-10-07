% This code is for visualisation of stimulus-aligned average raster, intended to cover the entire trial (~3s)

% Morgane September 2019

% close all
% clear all

%DMS 

      load('BehPhotoM_Exp23_DMS.mat')
 Animals = [53, 62, 63, 71,72]

 brain_region = 'DMS'

% NAC
 if exist('brain_region', 'var') && strcmp(brain_region, 'NAc')
     clearvars -except BehPhotoM
  else clear all
     load('BehPhotoM_Exp23_NAc.mat')
 end
Animals = [56 57 59 66]
brain_region = 'NAc'

% VTA: 
%  if exist('brain_region', 'var') && strcmp(brain_region, 'VTA')
%      clearvars -except BehPhotoM
%  else clear all
%       load('BehPhotoM_Exp23_VTA.mat')
%  end
% Animals = [48 50 51 64]
% brain_region = 'VTA'

stim2plot = 0.5;

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 150;

StartTime = 3700; % saved in the database.
RT_min = 0.2;
RT_max = 2.1; % reaction time range to include

FT_min = 0.1;
FT_max = 0.5; % outcome time range to include (action -> outcome)



sStart = -0.1; %stimulus
sStop = 3;

aStart = -2.3; 
aStop = 0.8;

oStart = -2.6;
oStop = 0.5;

downSample = 1200; % sampling rate AFTER DOWNSAMPLING in makeDataSummaryMatrix
eventOnset = 3700;  % in the saved matrix, the event occurs in this column

% --- get and organise data -----------------------------------------------
[LargeSmallErrorColor] = getColors();

LargeStimDataContra = [];
SmallStimDataContra = [];
ErrorStimDataContra = [];
LargeActionDataContra = [];
SmallActionDataContra = [];
ErrorActionDataContra = [];
LargeOutcomeDataContra = [];
SmallOutcomeDataContra = [];
ErrorOutcomeDataContra = [];

PopLargeStimDataContra = zeros(1,13100);
PopSmallStimDataContra = zeros(1,13100);
PopErrorStimDataContra = zeros(1,13100);
PopLargeActionDataContra = zeros(1,13100);
PopSmallActionDataContra = zeros(1,13100);
PopErrorActionDataContra = zeros(1,13100);
PopLargeOutcomeDataContra = zeros(1,13100);
PopSmallOutcomeDataContra = zeros(1,13100);
PopErrorOutcomeDataContra = zeros(1,13100);

SActionTimes_Contra = {NaN, NaN, NaN}; 
SOutcomeTimes_Contra = {NaN, NaN, NaN};
AStimulusTimes_Contra = {NaN, NaN, NaN}; 
AOutcomeTimes_Contra = {NaN, NaN, NaN};
OStimulusTimes_Contra = {NaN, NaN, NaN}; 
OActionTimes_Contra = {NaN, NaN, NaN};

chan_count = 0;
for iAnimal = Animals
    
    if isempty(BehPhotoM(iAnimal).GrandSummaryR)
        Chan = 1;
    elseif isempty(BehPhotoM(iAnimal).GrandSummaryL)
        Chan = 2;
    else
        Chan = [1 2];
    end
    
    
    
    
    for iSession = 1:length(BehPhotoM(iAnimal).Session)
        BehData = BehPhotoM(iAnimal).Session(iSession).TrialTimingData;
        RTs = BehData(:,10) - BehData(:,13);
        FTs = BehData(:,14) - BehData(:,10);
        OTs = BehData(:,14) - BehData(:,13);
        
        RTExcludeTrials = unique([find(isnan(RTs));(find(RTs < RT_min)) ; (find(RTs >RT_max)) ; (find(FTs < FT_min)) ; (find(FTs > FT_max))]);
        
        correctTrials = find(BehData(:,9)==1);
        errorTrials = find(BehData(:,9)==0);
        toLargeRewardTrials = [intersect(find(BehData(:,3)==-1), find(BehData(:,8)==1)); intersect(find(BehData(:,3)==1), find(BehData(:,8)==2))];
        toSmallRewardTrials = [intersect(find(BehData(:,3)==1), find(BehData(:,8)==1)); intersect(find(BehData(:,3)==-1), find(BehData(:,8)==2))];
        
        for iChan = Chan
            if iChan == 1
                contraTrials = setdiff(find(BehData(:,2)==stim2plot), RTExcludeTrials);
                

                LargeStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimL(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimL(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimL(intersect(contraTrials, errorTrials),:),1);
                
                LargeActionDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronActionL(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallActionDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronActionL(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorActionDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronActionL(intersect(contraTrials, errorTrials),:),1);

                LargeOutcomeDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronRewardL(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallOutcomeDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronRewardL(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorOutcomeDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronRewardL(intersect(contraTrials, errorTrials),:),1);            
            
            elseif iChan == 2
                contraTrials = setdiff(find(BehData(:,2)==-stim2plot), RTExcludeTrials);
                
                LargeStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimR(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimR(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimR(intersect(contraTrials, errorTrials),:),1);
                
                LargeActionDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronActionR(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallActionDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronActionR(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorActionDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronActionR(intersect(contraTrials, errorTrials),:),1);
                
                LargeOutcomeDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronRewardR(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallOutcomeDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronRewardR(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorOutcomeDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronRewardR(intersect(contraTrials, errorTrials),:),1);                

            end
            
            
            largesmallerrorindex = {mintersect(contraTrials, correctTrials, toLargeRewardTrials), mintersect(contraTrials, correctTrials, toSmallRewardTrials), intersect(contraTrials, errorTrials)};
            
            for i = 1:3
                
            SActionTimes_Contra{i} = [SActionTimes_Contra{i}; BehData(largesmallerrorindex{i},10) - BehData(largesmallerrorindex{i},13)]; % action times relative to stimulus
            SOutcomeTimes_Contra{i} = [SOutcomeTimes_Contra{i}; BehData(largesmallerrorindex{i},14) - BehData(largesmallerrorindex{i},13)]; % outcome times relative to stimulus
            
            AStimulusTimes_Contra{i} = [AStimulusTimes_Contra{i}; BehData(largesmallerrorindex{i},13) - BehData(largesmallerrorindex{i},10)];
            AOutcomeTimes_Contra{i} = [AOutcomeTimes_Contra{i}; BehData(largesmallerrorindex{i},14) - BehData(largesmallerrorindex{i},10)];
            
            OStimulusTimes_Contra{i} = [OStimulusTimes_Contra{i}; BehData(largesmallerrorindex{i},13) - BehData(largesmallerrorindex{i},14)];
            OActionTimes_Contra{i} = [OActionTimes_Contra{i}; BehData(largesmallerrorindex{i},10) - BehData(largesmallerrorindex{i},14)];
            
            end
            
            
            PopLargeStimDataContra = nansum([PopLargeStimDataContra ;LargeStimDataContra]);
            PopSmallStimDataContra = nansum([PopSmallStimDataContra ;SmallStimDataContra]);
            PopErrorStimDataContra = nansum([PopErrorStimDataContra ; ErrorStimDataContra]);
            
            PopLargeActionDataContra = nansum([PopLargeActionDataContra ; LargeActionDataContra]);
            PopSmallActionDataContra = nansum([PopSmallActionDataContra ; SmallActionDataContra]);
            PopErrorActionDataContra = nansum([PopErrorActionDataContra ; ErrorActionDataContra]);

            PopLargeOutcomeDataContra = nansum([PopLargeOutcomeDataContra ; LargeOutcomeDataContra]);
            PopSmallOutcomeDataContra = nansum([PopSmallOutcomeDataContra ; SmallOutcomeDataContra]);
            PopErrorOutcomeDataContra = nansum([PopErrorOutcomeDataContra ; ErrorOutcomeDataContra]);
        
        
        
        end
        
    end
    
end

maxDenomStim = max(max([PopLargeStimDataContra; PopSmallStimDataContra; PopErrorStimDataContra]));
PopLargeStimDataContra = PopLargeStimDataContra ./ maxDenomStim; 
PopSmallStimDataContra = PopSmallStimDataContra ./ maxDenomStim;
PopErrorStimDataContra = PopErrorStimDataContra ./ maxDenomStim; 

maxDenomAction = max(max([PopLargeActionDataContra; PopSmallActionDataContra; PopErrorActionDataContra]));
PopLargeActionDataContra = PopLargeActionDataContra ./ maxDenomAction; 
PopSmallActionDataContra = PopSmallActionDataContra ./ maxDenomAction;
PopErrorActionDataContra = PopErrorActionDataContra ./ maxDenomAction; 

maxDenomOutcome = max(max([PopLargeOutcomeDataContra; PopSmallOutcomeDataContra; PopErrorOutcomeDataContra]));
PopLargeOutcomeDataContra = PopLargeOutcomeDataContra ./ maxDenomOutcome; 
PopSmallOutcomeDataContra = PopSmallOutcomeDataContra ./ maxDenomOutcome;
PopErrorOutcomeDataContra = PopErrorOutcomeDataContra ./ maxDenomOutcome; 



for i = 1:3

    SActionTimes_Contra{i} = downSample*(abs(sStart)+SActionTimes_Contra{i});
    SOutcomeTimes_Contra{i} = downSample*(abs(sStart)+SOutcomeTimes_Contra{i}-0.07);
    AStimulusTimes_Contra{i} = downSample*(abs(aStart)+AStimulusTimes_Contra{i});
    AOutcomeTimes_Contra{i} = downSample*(abs(aStart)+AOutcomeTimes_Contra{i}-0.07);
    OActionTimes_Contra{i} = downSample*(abs(oStart)+OActionTimes_Contra{i});
    OStimulusTimes_Contra{i} = downSample*(abs(oStart)+OStimulusTimes_Contra{i}+0.07);
    
    
end


%%
figure; 


plots(1) = subplot(3,1,1); % stimulus-aligned 
plot(smooth(PopLargeStimDataContra(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', 'k')
xlabel('Time from stimulus(s)')
xticks([0 abs(sStart*downSample) abs(sStart-sStop)*downSample])
xticklabels([sStart 0 sStop])
line([abs(sStart*downSample) abs(sStart*downSample)], [0 1], 'Color', [0.5 0.5 0.5] ,'LineStyle','--', 'LineWidth', 2)
hold on;
yyaxis right
histogram(SActionTimes_Contra{1}, 'BinWidth', 90, 'FaceColor', [0 0.58 0.77] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
histogram(SOutcomeTimes_Contra{1}, 'BinWidth', 90, 'FaceColor', [0.27 0.74 0], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
ylim([0 400])
hold on; 
yyaxis left
txt = '\nabla';
text(nanmedian(SActionTimes_Contra{1}), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian(SActionTimes_Contra{1})))+0.2, txt)
txt2 = '\nabla';
text(nanmedian(SOutcomeTimes_Contra{1}), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian(SOutcomeTimes_Contra{1})))+0.2, txt2)



plots(2) = subplot(3,1,2); % action-aligned 
plot(smooth(PopLargeActionDataContra(eventOnset+(aStart*downSample):eventOnset+(aStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', 'k')
xticks([0 abs(aStart*downSample) abs(aStart-aStop)*downSample])
xticklabels([aStart 0 aStop])
line([abs(aStart*downSample) abs(aStart*downSample)], [0 1], 'Color', [0 0.58 0.77] ,'LineStyle','--', 'LineWidth', 2)
hold on;
xlabel('Time from Action(s)')
yyaxis right
histogram([AStimulusTimes_Contra{1}], 'BinWidth', 90, 'FaceColor', [0.5 0.5 0.5] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
histogram([AOutcomeTimes_Contra{1}], 'BinWidth', 90, 'FaceColor', [0.27 0.74 0], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
ylim([0 400])
hold on; 
yyaxis left
txt = '\nabla';
text(nanmedian([AStimulusTimes_Contra{1}]), PopLargeActionDataContra(eventOnset+(aStart*downSample)+floor(nanmedian([AStimulusTimes_Contra{1}])))+0.2, txt)
txt2 = '\nabla';
text(nanmedian([AOutcomeTimes_Contra{1}]), PopLargeActionDataContra(eventOnset+(aStart*downSample)+floor(nanmedian([AOutcomeTimes_Contra{1}])))+0.2, txt2)


plots(3) = subplot(3,1,3); %outcome-aligned 
plot(smooth(PopLargeOutcomeDataContra(eventOnset+(oStart*downSample):eventOnset+(oStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', 'k')
xticks([0 abs(oStart*downSample) abs(oStart-oStop)*downSample])
xticklabels([oStart 0 oStop])
line([abs(oStart*downSample) abs(oStart*downSample)], [0 1], 'Color', [0.27 0.74 0] ,'LineStyle','--', 'LineWidth', 2)
xlabel('Time from Outcome(s)')
hold on;
yyaxis right
histogram(OStimulusTimes_Contra{1}, 'BinWidth', 90, 'FaceColor', [0.5 0.5 0.5] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
histogram(OActionTimes_Contra{1}, 'BinWidth', 90, 'FaceColor', [0 0.58 0.77], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
ylim([0 400])
hold on; 
yyaxis left
txt = '\nabla';
text(nanmedian(OStimulusTimes_Contra{1}), PopLargeOutcomeDataContra(eventOnset+(sStart*downSample)+floor(nanmedian(OStimulusTimes_Contra{1})))+0.2, txt)
txt2 = '\nabla';
text(nanmedian(OActionTimes_Contra{1}), PopLargeOutcomeDataContra(eventOnset+(sStart*downSample)+floor(nanmedian(OActionTimes_Contra{1})))+0.2, txt2)



function [LargeSmallErrorColor] = getColors()

LargeSmallErrorColor = [0.15 0.4 0
                        0.15 0.75 0
                        1 0 0];                

end

function totalChannels = getTotalChanN(region)

if strcmpi(region, 'DMS')
    totalChannels = 7;
elseif strcmpi(region, 'NAC')
    totalChannels = 5;
elseif strcmpi(region, 'VTA')
    totalChannels = 4;
end
end
