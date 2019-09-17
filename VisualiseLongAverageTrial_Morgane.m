% This code is for visualisation of stimulus-aligned average raster, intended to cover the entire trial (~3s)

% Morgane September 2019

% close all
% clear all

%DMS 
 if exist('brain_region', 'var') && strcmp(brain_region, 'DMS')
     clearvars -except BehPhotoM
 else clear all
 end
 Animals = [53, 62, 63, 71,72]
 brain_region = 'DMS'

% NAC
%  if exist('brain_region', 'var') && strcmp(brain_region, 'NAc')
%      clearvars -except BehPhotoM
%   else clear all
%  end
% Animals = [56 57 59 66]
% brain_region = 'NAc'

% VTA: 
%  if exist('brain_region', 'var') && strcmp(brain_region, 'VTA')
%      clearvars -except BehPhotoM
%  else clear all
%  end
% Animals = [48 50 51 64]
% brain_region = 'VTA'

stim2plot = 0.5;

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 150;



StartTime = 3700; % saved in the database.
RT_min = 0.2;
RT_max = 3; % reaction time range to include

OT_min = 0.2;
OT_max = 3; % outcome time range to include

sStart = -0.1; %stimulus
sStop = 3;

downSample = 1200; % sampling rate AFTER DOWNSAMPLING in makeDataSummaryMatrix
eventOnset = 3700;  % in the saved matrix, the event occurs in this column

% --- get and organise data -----------------------------------------------
[LargeSmallErrorColor] = getColors();

LargeStimDataContra = [];
SmallStimDataContra = [];
ErrorStimDataContra = [];

PopLargeStimDataContra = zeros(1,13100);
PopSmallStimDataContra = zeros(1,13100);
PopErrorStimDataContra = zeros(1,13100);

ActionTimes_Contra = {NaN, NaN, NaN}; 
OutcomeTimes_Contra = {NaN, NaN, NaN};

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
        OTs = BehData(:,14) - BehData(:,13);
        
        RTExcludeTrials = unique([(find(RTs < RT_min)) ; (find(RTs >RT_max)) ; (find(OTs < OT_min)) ; (find(OTs > OT_max))]);
        
        correctTrials = find(BehData(:,9)==1);
        errorTrials = find(BehData(:,9)==0);
        toLargeRewardTrials = [intersect(find(BehData(:,3)==-1), find(BehData(:,8)==1)); intersect(find(BehData(:,3)==1), find(BehData(:,8)==2))];
        toSmallRewardTrials = [intersect(find(BehData(:,3)==1), find(BehData(:,8)==1)); intersect(find(BehData(:,3)==-1), find(BehData(:,8)==2))];
        
        for iChan = Chan
            if iChan == 1
                contraTrials = setdiff(find(BehData(:,2)==-stim2plot), RTExcludeTrials);
                
                LargeStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimL(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimL(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimL(intersect(contraTrials, errorTrials),:),1);
            elseif iChan == 2
                contraTrials = setdiff(find(BehData(:,2)==stim2plot), RTExcludeTrials);
                
                LargeStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimR(mintersect(contraTrials, correctTrials, toLargeRewardTrials),:),1);
                SmallStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimR(mintersect(contraTrials, correctTrials, toSmallRewardTrials),:),1);
                ErrorStimDataContra = mean(BehPhotoM(iAnimal).Session(iSession).NeuronStimR(intersect(contraTrials, errorTrials),:),1);
            end
            
            
            ActionTimes_Contra{1} = [ActionTimes_Contra{1}; BehData(mintersect(contraTrials, correctTrials, toLargeRewardTrials),14) - BehData(mintersect(contraTrials, correctTrials, toLargeRewardTrials),13)]; % large rwd column
            ActionTimes_Contra{2} = [ActionTimes_Contra{2}; BehData(mintersect(contraTrials, correctTrials, toSmallRewardTrials),14) - BehData(mintersect(contraTrials, correctTrials, toSmallRewardTrials),13)]; % small rwd column
            ActionTimes_Contra{3} = [ActionTimes_Contra{3}; BehData(intersect(contraTrials, errorTrials),14) - BehData(mintersect(contraTrials, errorTrials),13)]; % error column
            OutcomeTimes_Contra{1} = [OutcomeTimes_Contra{1}; BehData(mintersect(contraTrials, correctTrials, toLargeRewardTrials),10) - BehData(mintersect(contraTrials, correctTrials, toLargeRewardTrials),13)]; % large rwd column
            OutcomeTimes_Contra{2} = [OutcomeTimes_Contra{2}; BehData(mintersect(contraTrials, correctTrials, toSmallRewardTrials),10) - BehData(mintersect(contraTrials, correctTrials, toSmallRewardTrials),13)]; % small rwd column
            OutcomeTimes_Contra{3} = [OutcomeTimes_Contra{3}; BehData(intersect(contraTrials, errorTrials),10) - BehData(mintersect(contraTrials, errorTrials),13)]; % error column
            
            
            PopLargeStimDataContra = nansum([PopLargeStimDataContra ;LargeStimDataContra]);
            PopSmallStimDataContra = nansum([PopSmallStimDataContra ;SmallStimDataContra]);
            PopErrorStimDataContra = nansum([PopErrorStimDataContra ; ErrorStimDataContra]);
            
        end
        
    end
    
end

maxDenom = max(max([PopLargeStimDataContra; PopLargeStimDataContra; PopErrorStimDataContra]));
PopLargeStimDataContra = PopLargeStimDataContra ./ maxDenom; 
PopSmallStimDataContra = PopSmallStimDataContra ./ maxDenom;
PopErrorStimDataContra = PopErrorStimDataContra ./ maxDenom; 

for i = 1:3

    ActionTimes_Contra{i} = downSample*(abs(sStart)+ActionTimes_Contra{i});
    OutcomeTimes_Contra{i} = downSample*(abs(sStart)+OutcomeTimes_Contra{i}-0.07);
    
end


%%
figure; 

plots(1) = subplot(3,1,1); % large correct only 
plot(smooth(PopLargeStimDataContra(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', LargeSmallErrorColor(1,:))
hold on; 

yyaxis right
histogram(ActionTimes_Contra{1}, 'BinWidth', 90, 'FaceColor', [0.27 0.74 0] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
histogram(OutcomeTimes_Contra{1}, 'BinWidth', 90, 'FaceColor', [0 0.58 0.77], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
ylim([0 600])
hold on; 
yyaxis left
txt = '\nabla';
text(nanmedian(ActionTimes_Contra{1}), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian(ActionTimes_Contra{1})))+0.2, txt)
txt2 = '\nabla';
text(nanmedian(OutcomeTimes_Contra{1}), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian(OutcomeTimes_Contra{1})))+0.2, txt2)


plots(2) = subplot(3,1,2); % large and error
plot(smooth(PopErrorStimDataContra(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', LargeSmallErrorColor(3,:))
hold on;
plot(smooth(PopLargeStimDataContra(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', LargeSmallErrorColor(1,:))
ylabel('\Delta F / F')

yyaxis right
histogram([ActionTimes_Contra{1}; ActionTimes_Contra{3}], 'BinWidth', 90, 'FaceColor', [0.27 0.74 0] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
histogram([OutcomeTimes_Contra{1}; OutcomeTimes_Contra{3}], 'BinWidth', 90, 'FaceColor', [0 0.58 0.77], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
ylim([0 600])
hold on; 
yyaxis left
txt = '\nabla';
text(nanmedian([ActionTimes_Contra{1}; ActionTimes_Contra{3}]), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian([ActionTimes_Contra{1};ActionTimes_Contra{3}])))+0.2, txt)
txt2 = '\nabla';
text(nanmedian([OutcomeTimes_Contra{1}; OutcomeTimes_Contra{3}]), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian([OutcomeTimes_Contra{1}; OutcomeTimes_Contra{3}])))+0.2, txt2)



plots(3) = subplot(3,1,3); %large, small, and error
plot(smooth(PopErrorStimDataContra(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', LargeSmallErrorColor(3,:))
hold on;
plot(smooth(PopSmallStimDataContra(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', LargeSmallErrorColor(2,:))
plot(smooth(PopLargeStimDataContra(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample)), smooth_factor), 'LineWidth', 2, 'Color', LargeSmallErrorColor(1,:))

yyaxis right
histogram([ActionTimes_Contra{1}; ActionTimes_Contra{2}; ActionTimes_Contra{3}], 'BinWidth', 90, 'FaceColor', [0.27 0.74 0] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
histogram([OutcomeTimes_Contra{1}; OutcomeTimes_Contra{2}; OutcomeTimes_Contra{3}], 'BinWidth', 90, 'FaceColor', [0 0.58 0.77], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
ylim([0 600])
legend('Error', 'Small reward', 'Large reward')
xlabel('Time from stimulus(s)')
hold on; 
yyaxis left
txt = '\nabla';
text(nanmedian([ActionTimes_Contra{1}; ActionTimes_Contra{2}; ActionTimes_Contra{3}]), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian([ActionTimes_Contra{1}; ActionTimes_Contra{2}; ActionTimes_Contra{3}])))+0.2, txt)
txt2 = '\nabla';
text(nanmedian([OutcomeTimes_Contra{1}; OutcomeTimes_Contra{2}; OutcomeTimes_Contra{3}]), PopLargeStimDataContra(eventOnset+(sStart*downSample)+floor(nanmedian([OutcomeTimes_Contra{1}; OutcomeTimes_Contra{2}; OutcomeTimes_Contra{3}])))+0.2, txt2)




set(plots, 'XTick', ([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]), ...
           'XTickLabel', ([sStart 0 1 2 sStop]))




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
