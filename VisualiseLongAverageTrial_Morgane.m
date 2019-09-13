% creates 2 rasters for contra and ipsi (per animal or per session) ordered by
% Error/Small/Large and by RTs. Works for unihem and bihem animals
% Morgane April 2019

% can you see this? yes

% for information about normalisation, check end of this code. 

% close all
clear all



% ----- enter reqs --------------------------------------------------------

%DMS 
%  if exist('brain_region', 'var') && strcmp(brain_region, 'DMS')
%      clearvars -except BehPhotoM
%  else clear all
%  end
%  Animals = [53, 62, 63, 71,72]
%  brain_region = 'DMS'

% NAC
 if exist('brain_region', 'var') && strcmp(brain_region, 'NAc')
     clearvars -except BehPhotoM
  else clear all
 end
Animals = [56 57 59 66]
brain_region = 'NAc'

% VTA: 
%  if exist('brain_region', 'var') && strcmp(brain_region, 'VTA')
%      clearvars -except BehPhotoM
%  else clear all
%  end
% Animals = [48 50 51 64]
% brain_region = 'VTA'
% ------ -----------------------------------



exp_ID = '23'
plotRasters = 0;

concatenate = 1; %show all trials across all sessions for an animal
IpsiContra = 1; %separate rasters based on ipsi / contra? 

stim_2_plot = 0.5; %should be positive
smooth_factor = 100;

[ActionColor, IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors();
% colorRange= getColorRange(animal_ID, exp_ID); % color range for plotting imagesc data

RT_min = 0.2;
RT_max = 3; % reaction time range to include

OT_min = 0.2;
OT_max = 3; % outcome time range to include 

% ------------------ start stop times for task events in second ------------------

sStart = -0.1; %stimulus
sStop = 3;

aStart = -0.6; %action
aStop = 0.2;

rStart = -0.2; %reward
rStop = 0.6;

downSample = 1200; % sampling rate AFTER DOWNSAMPLING in makeDataSummaryMatrix
eventOnset = 3700;  % in the saved matrix, the event occurs in this column

% -------------------------------------------------------------------------


load(['BehPhotoM_Exp',exp_ID,'_',brain_region,'.mat'])

fullStim = [-stim_2_plot stim_2_plot];

animal_count = 0;

for animal_ID = Animals
    animal_count = animal_count +1;
nSessions = 1:length(BehPhotoM(animal_ID).Session);

for iSession = nSessions
    
    if animal_count ==1 && iSession ==1
        leftStimTrials = [];
        rightStimTrials = [];
        
        BehDataIpsi     = [];
        BehDataContra   = [];
        StimDataIpsi    = [];  
        StimDataContra  = [];     

        errorTrialsIpsi      = [];
        errorTrialsContra    = [];
        smallRewTrialsIpsi   = [];
        smallRewTrialsContra = [];
        largeRewTrialsIpsi   = [];
        largeRewTrialsContra = [];
        
        RTs = [];
        RTExcludeTrials = [];
        
    end
    

    % this is interval between outcome and stimulus 
    RTs = BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,10) - BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,13);
    OTs = BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,14) - BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,13);
    
    RTExcludeTrials = [(find(RTs < RT_min)) ; (find(RTs >RT_max)) ; (find(OTs < OT_min)) ; (find(OTs > OT_max))];

    leftStimTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2)== -(stim_2_plot)),RTExcludeTrials);
    rightStimTrials = setdiff(find(BehPhotoM(animal_ID).Session(iSession).TrialTimingData(:,2)==stim_2_plot), RTExcludeTrials);
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimL')%load hemispheric deltaFoverF. add Action/Reward to the following two 'if's if desired. 
        StimDataIpsi        = [StimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimL(leftStimTrials,:)];
        StimDataContra      = [StimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimL(rightStimTrials,:)];
        BehDataIpsi         = [BehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(leftStimTrials,:)];
        BehDataContra       = [BehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(rightStimTrials,:)];
    end
    
    if isfield(BehPhotoM(animal_ID).Session, 'NeuronStimR') 
        StimDataIpsi        = [StimDataIpsi;      BehPhotoM(animal_ID).Session(iSession).NeuronStimR(rightStimTrials,:)];
        StimDataContra      = [StimDataContra;    BehPhotoM(animal_ID).Session(iSession).NeuronStimR(leftStimTrials,:)];
        BehDataIpsi         = [BehDataIpsi;       BehPhotoM(animal_ID).Session(iSession).TrialTimingData(rightStimTrials,:)];
        BehDataContra       = [BehDataContra;     BehPhotoM(animal_ID).Session(iSession).TrialTimingData(leftStimTrials,:)];
    end
    
end
end
        
        %ipsi and contra error / small rwd / large rwd indexing
        errorTrialsIpsi(:,2) = find(BehDataIpsi(:,9)==0);
        errorTrialsContra(:,2) = find(BehDataContra(:,9)==0);
        
        if strcmp(exp_ID, '7')
            
            smallRewTrialsIpsi(:,2)     = [find(BehDataIpsi(:,9)==1)];
            
            smallRewTrialsContra(:,2)   = [find(BehDataContra(:,9)==1)];
            
            
        elseif strcmp(exp_ID, '23')
            
            largeRewTrialsIpsi(:,2)     = [intersect(find(BehDataIpsi(:,9)==1), intersect(find(BehDataIpsi(:,8)==1), find(BehDataIpsi(:,2)<0))); ...
                intersect(find(BehDataIpsi(:,9)==1), intersect(find(BehDataIpsi(:,8)==2), find(BehDataIpsi(:,2)>0)))];
            
            largeRewTrialsContra(:,2) = [intersect(find(BehDataContra(:,9)==1), intersect(find(BehDataContra(:,8)==1), find(BehDataContra(:,2)<0))); ...
                intersect(find(BehDataContra(:,9)==1), intersect(find(BehDataContra(:,8)==2), find(BehDataContra(:,2)>0)))];
            
            smallRewTrialsIpsi(:,2) = [intersect(find(BehDataIpsi(:,9)==1), intersect(find(BehDataIpsi(:,8)==2), find(BehDataIpsi(:,2)<0))); ...
                intersect(find(BehDataIpsi(:,9)==1), intersect(find(BehDataIpsi(:,8)==1), find(BehDataIpsi(:,2)>0)))];
            
            smallRewTrialsContra(:,2) = [intersect(find(BehDataContra(:,9)==1), intersect(find(BehDataContra(:,8)==2), find(BehDataContra(:,2)<0))); ...
                intersect(find(BehDataContra(:,9)==1), intersect(find(BehDataContra(:,8)==1), find(BehDataContra(:,2)>0)))];
        end
        
    
       
        %get RTs for different trial types and sort by this
        errorTrialsIpsi(:,1) = BehDataIpsi(errorTrialsIpsi(:,2), 10) - BehDataIpsi(errorTrialsIpsi(:,2), 13);
        errorTrialsIpsi = sortrows(errorTrialsIpsi);
                errorTrialsIpsiOutcome(:,1) = BehDataIpsi(errorTrialsIpsi(:,2), 14) - BehDataIpsi(errorTrialsIpsi(:,2), 13);

        errorTrialsContra(:,1) = BehDataContra(errorTrialsContra(:,2), 10) - BehDataContra(errorTrialsContra(:,2), 13);
        errorTrialsContra = sortrows(errorTrialsContra); 
                errorTrialsContraOutcome(:,1) = BehDataContra(errorTrialsContra(:,2), 14) - BehDataContra(errorTrialsContra(:,2), 13);

        
        smallRewTrialsIpsi(:,1) = BehDataIpsi(smallRewTrialsIpsi(:,2), 10) - BehDataIpsi(smallRewTrialsIpsi(:,2), 13);
        smallRewTrialsIpsi = sortrows(smallRewTrialsIpsi);
        smallRewTrialsIpsiOutcome(:,1) = BehDataIpsi(smallRewTrialsIpsi(:,2), 14) - BehDataIpsi(smallRewTrialsIpsi(:,2), 13);
        
        smallRewTrialsContra(:,1) = BehDataContra(smallRewTrialsContra(:,2), 10) - BehDataContra(smallRewTrialsContra(:,2), 13);
        smallRewTrialsContra = sortrows(smallRewTrialsContra);
                smallRewTrialsContraOutcome(:,1) = BehDataContra(smallRewTrialsContra(:,2), 14) - BehDataContra(smallRewTrialsContra(:,2), 13);

        
        if strcmp(exp_ID, '23')
            largeRewTrialsIpsi(:,1) = BehDataIpsi(largeRewTrialsIpsi(:,2), 10) - BehDataIpsi(largeRewTrialsIpsi(:,2), 13);
            largeRewTrialsIpsi = sortrows(largeRewTrialsIpsi);
            largeRewTrialsIpsiOutcome(:,1) = BehDataIpsi(largeRewTrialsIpsi(:,2), 14) - BehDataIpsi(largeRewTrialsIpsi(:,2), 13);

            
            largeRewTrialsContra(:,1) = BehDataContra(largeRewTrialsContra(:,2), 10) - BehDataContra(largeRewTrialsContra(:,2), 13);
            largeRewTrialsContra = sortrows(largeRewTrialsContra);
            largeRewTrialsContraOutcome(:,1) = BehDataContra(largeRewTrialsContra(:,2), 14) - BehDataContra(largeRewTrialsContra(:,2), 13);
            
            largeRewTrialsContraAction(:,1) = BehDataContra(largeRewTrialsContra(:,2), 10) - BehDataContra(largeRewTrialsContra(:,2), 13);

        end
        
        % plotting order
        if strcmp(exp_ID, '7')
            plotIndexIpsi = [errorTrialsIpsi(:,2); smallRewTrialsIpsi(:,2)];
            plotIndexContra = [errorTrialsContra(:,2); errorTrialsContra(:,2)];
            
            RTIpsi = [errorTrialsIpsi(:,1); smallRewTrialsIpsi(:,1)];
            RTContra = [errorTrialsContra(:,1); smallRewTrialsContra(:,1)];
            
            TIpsiOutcome = [errorTrialsIpsiOutcome(:,1); smallRewTrialsIpsiOutcome(:,1)];
            TContraOutcome = [errorTrialsContraOutcome(:,1); smallRewTrialsContraOutcome(:,1)];
            
        elseif strcmp(exp_ID, '23')
            plotIndexIpsi = [errorTrialsIpsi(:,2); smallRewTrialsIpsi(:,2); largeRewTrialsIpsi(:,2)];
            plotIndexContra = [errorTrialsContra(:,2); smallRewTrialsContra(:,2); largeRewTrialsContra(:,2)];
            
            RTIpsi = [errorTrialsIpsi(:,1); smallRewTrialsIpsi(:,1); largeRewTrialsIpsi(:,1)];
            RTContra = [errorTrialsContra(:,1); smallRewTrialsContra(:,1); largeRewTrialsContra(:,1)];
            
            TIpsiOutcome = [errorTrialsIpsiOutcome(:,1); smallRewTrialsIpsiOutcome(:,1); largeRewTrialsIpsiOutcome(:,1)];
            TContraOutcome = [errorTrialsContraOutcome(:,1); smallRewTrialsContraOutcome(:,1); largeRewTrialsContraOutcome(:,1)];
        end
        
        

        actionTimesIpsi = abs(sStart*downSample) + (RTIpsi*downSample);
        actionTimesContra = abs(sStart*downSample) + (RTContra*downSample);
        actionTimesContra = abs(sStart*downSample) + (largeRewTrialsContraAction(:,1)*downSample);
        
        outcomeTimesIpsi = abs(sStart*downSample) + (TIpsiOutcome*downSample);
        outcomeTimesContra = abs(sStart*downSample) + (TContraOutcome*downSample);
        outcomeTimesContra = abs(sStart*downSample) + (largeRewTrialsContraOutcome(:,1)*downSample);
        
        
        
    % ------------------- plot  ----------------------------------
%     if concatenate == 0 || iSession == max(nSessions)
        
        M = [];
        M = {StimDataContra, StimDataIpsi};
        A = [];
        A = {actionTimesContra, actionTimesIpsi};
        Ar = [];
        Ar = {outcomeTimesContra, outcomeTimesIpsi};
       
        P = [];
        P = {plotIndexContra, plotIndexIpsi};
        E = [];
        E = {errorTrialsContra, errorTrialsIpsi};
        R = [];
        R = {smallRewTrialsContra, smallRewTrialsIpsi};
   
        if plotRasters 
            
        figure;
        
        for iSubplot = 1:2
            data = M{iSubplot}; % looping through StimDataContra / Ipsi
            actionTimes = A{iSubplot};
            outcomeTimes = Ar{iSubplot};
            
            plotIndex = P{iSubplot};
            errorTrials = E{iSubplot};
            smallRewTrials = R{iSubplot};
            
            subplot(1,2,iSubplot)
            imagesc(smooth2a(data(plotIndex, (eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))), 0, smooth_factor), colorRange)
            colormap('bluewhitered')
            hold on
            for ievent = 1:length(actionTimes)
                if ievent <= length(errorTrials)
                    plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0.7 0 0.2 ], 'MarkerFaceColor', [0.7 0 0.2 ], 'MarkerSize', 3)
                elseif length(errorTrials) < ievent && ievent <= length(errorTrials)+length(smallRewTrials)
                    plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0 1 0 ], 'MarkerFaceColor', [0 1 0 ], 'MarkerSize', 3)
                elseif ievent > length(errorTrials)+length(smallRewTrials)
                    plot(actionTimes(ievent),ievent, 'o', 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 3)
                end
            end
            
              for ievent = 1:length(outcomeTimes)
                if ievent <= length(errorTrials)
                    plot(outcomeTimes(ievent),ievent, 'x', 'MarkerEdgeColor', [0.7 0 0.2 ], 'MarkerFaceColor', [0.7 0 0.2 ], 'MarkerSize', 5)
                elseif length(errorTrials) < ievent && ievent <= length(errorTrials)+length(smallRewTrials)
                    plot(outcomeTimes(ievent),ievent, 'x', 'MarkerEdgeColor', [0 1 0 ], 'MarkerFaceColor', [0 1 0 ], 'MarkerSize', 5)
                elseif ievent > length(errorTrials)+length(smallRewTrials)
                    plot(outcomeTimes(ievent),ievent, 'x', 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 5)
                end
            end
            xticks([1 abs(sStart*downSample) (sStop-sStart)*downSample-1])
            xticklabels([sStart 0 sStop])
            ylabel('Trials')
            xlabel('Time from stimulus (s)')
            line([abs(sStart*downSample) abs(sStart*downSample)], [0 length(plotIndex)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 1.5); % stim onset line
            
        end
%     end
    
    end


%% average signal through trial aligned  at stimulus onset. large reward only

        actionTimesContra = [];
        actionTimesContra = abs(sStart)*downSample + (largeRewTrialsContraAction(:,1))*downSample;
        outcomeTimesContra = [];
        outcomeTimesContra = abs(sStart)*downSample + (largeRewTrialsContraOutcome(:,1))*downSample;
        
        if strcmpi(brain_region, 'VTA')
            outcomeTimesContra = outcomeTimesContra - (0.07*downSample); % to account for feedback delivery delay 
        end
        
        
        figure; 

subplot(3,1,1)
data = smooth2a(M{1},0,smooth_factor);

yyaxis left
tempdata = mean(data(largeRewTrialsContra(:,2),(eventOnset+(sStart*downSample):eventOnset+(sStop*downSample))));
tempdata = tempdata./max(tempdata);
plot(tempdata, 'LineWidth', 1, 'Color', 'k');
xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) % 
xticklabels([sStart 0 1 2 sStop])
ylim([-0.2 1.1])
hold on; 
txt = '\nabla';
text(nanmedian(actionTimesContra), tempdata(floor(nanmedian(actionTimesContra)))+0.2, txt)
ax = gca; 
ax.TickDir = 'out';

hold on; 
txt2 = '\nabla';
text(nanmedian(outcomeTimesContra), tempdata(floor(nanmedian(outcomeTimesContra)))+0.2, txt2)

% subplot(3, 1, 2) % action time distribution
% title('Action times')
hold on;
yyaxis right
histogram(actionTimesContra, 'BinWidth', 90, 'FaceColor', [0 0.58 0.77] , 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
%xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) % 
%xticklabels([sStart 0 1 2 sStop])
ylim([0 400])

% subplot(3, 1, 3) % reward time distribution
% title('Reward times')
hold on; 
yyaxis right
histogram(outcomeTimesContra, 'BinWidth', 90, 'FaceColor', [0.27 0.74 0], 'FaceAlpha', 0.3 , 'EdgeAlpha', 0);
%xticks([1 abs(sStart*downSample) (abs(sStart)+1)*downSample (abs(sStart)+2)*downSample (sStop-sStart)*downSample-1]) % 
%xticklabels([sStart 0 1 2 sStop])


 linkaxes(get(gcf,'children'),'x')
 
 % alt action color = [0.07 0.62 1.00]
            
%% script-specficic functions

function [colorRange] = getColorRange(animal_ID, expID)

if animal_ID == 48
    if strcmp(expID, '23')
        [colorRange] = [-3 4]
    end
elseif animal_ID == 72
    if strcmp(expID, '23')
        [colorRange] = [-3 4]
    end
elseif animal_ID == 71
    if strcmp(expID, '23')
        [colorRange] = [-3 8]
    end
elseif animal_ID == 64
    if strcmp(expID, '23')
        [colorRange] = [-8 20]
    end
elseif animal_ID == 63
    if strcmp(expID, '23')
        [colorRange] = [-3 3]
    end
elseif animal_ID == 57
    if strcmp(expID, '23')
        [colorRange] = [-2 2]
    end
else [colorRange] = [-3 6]
    
end

end

function [ActionColor, IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors()
IpsiContraColor = [ 0 0 0
                    153/255 51/255 1
                    ];
                
ErrorCorrectColor = [   'r'
                        'g'
                        ];
                    
SmallLargeColor = [ 0.15 0.75 0
                    0.15 0.4 0
                    ];
                
ActionColor = [0.35 0.65 0.86];

end
