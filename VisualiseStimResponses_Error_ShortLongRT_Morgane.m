% Morgane October 2019
% This code visualises

% close all 
clear all

%DMS
brain_region = 'DMS';
Animals = [53, 62, 63, 71,72];
load('BehPhotoM_Exp23_DMS')


shortRT = 0.5; % short RTs are trials where the RT was less than this number of seconds 
longRT = 1.5; % long trials had a RT longer than this number of seconds
Stim2plot = 0.5;
smooth_factor = 200;
sampleRate = 1200;

if strcmp(getComputerName(), 'proxeddu') % morgane in new york
   addpath(genpath('D:\Morgane')); 
end



TimingVisualise = [-0.2 0.8
                   -0.7, 0.7
                   -0.2, 0.8]; % stim, action, reward in s
               
StartTime = 3700; % saved in the database.


%% --- get and organise data -----------------------------------------------
[IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors();

PopShortRTStimDataCorrect   = [];
PopLongRTStimDataCorrect    = [];
PopShortRTStimDataError     = [];
PopLongRTStimDataError      = [];

c = 1;
animal_count = 0;
chan_count = 0;
for iAnimal = Animals
    
    animal_count = animal_count + 1;
    
    if isempty(BehPhotoM(iAnimal).GrandSummaryR)        % uni hem animals L
        ChanN = 1;
    elseif isempty(BehPhotoM(iAnimal).GrandSummaryL)    % uni hem animals R
        ChanN = 2;  
    else                                                % bi hem animals
        ChanN = [1 2];
    end
    
    
    for iChan = ChanN
        chan_count = chan_count + 1;
        
        
        for iSession = 1:length(BehPhotoM(iAnimal).Session)
            
            tempBehData = BehPhotoM(iAnimal).Session(iSession).TrialTimingData;
            
            iShortRT = tempBehData(:,10)-tempBehData(:,13)< shortRT;
            iLongRT  = tempBehData(:,10)-tempBehData(:,13)> longRT;
            iLeftTrials = tempBehData(:,2) == -Stim2plot ;
            iRightTrials = tempBehData(:,2) == Stim2plot ;
            iCorrectTrials = (tempBehData(:,9)==1 & ( (tempBehData(:,8)==1 & tempBehData(:,3)==-1 ) | ( tempBehData(:,8)==2 & tempBehData(:,3)==1 )));
            iErrorTrials = tempBehData(:,9) == 0;
            

            if iChan == 1
                
                tempShortRTStimDataCorrect = BehPhotoM(iAnimal).Session(iSession).NeuronStimL(iShortRT & iRightTrials & iCorrectTrials,:);
                tempLongRTStimDataCorrect  = BehPhotoM(iAnimal).Session(iSession).NeuronStimL(iLongRT & iRightTrials & iCorrectTrials,:);
                
                tempShortRTStimDataError = BehPhotoM(iAnimal).Session(iSession).NeuronStimL(iShortRT & iRightTrials & iErrorTrials,:);
                tempLongRTStimDataError  = BehPhotoM(iAnimal).Session(iSession).NeuronStimL(iLongRT & iRightTrials & iErrorTrials,:);
                
            elseif iChan == 2
                
                tempShortRTStimDataCorrect = BehPhotoM(iAnimal).Session(iSession).NeuronStimR(iShortRT & iLeftTrials & iCorrectTrials,:);
                tempLongRTStimDataCorrect  = BehPhotoM(iAnimal).Session(iSession).NeuronStimR(iLongRT & iLeftTrials & iCorrectTrials,:);
                
                tempShortRTStimDataError = BehPhotoM(iAnimal).Session(iSession).NeuronStimR(iShortRT & iLeftTrials & iErrorTrials,:);
                tempLongRTStimDataError  = BehPhotoM(iAnimal).Session(iSession).NeuronStimR(iLongRT & iLeftTrials & iErrorTrials,:);
                
            end
            
            % normalise to 1. might be necessary to add minimum here - TBC
            tempDenom = max(max([tempShortRTStimDataCorrect; tempLongRTStimDataCorrect; tempShortRTStimDataError; tempLongRTStimDataError]));
            
            PopShortRTStimDataCorrect   = [PopShortRTStimDataCorrect; tempShortRTStimDataCorrect./ tempDenom];
            PopLongRTStimDataCorrect    = [PopLongRTStimDataCorrect; tempLongRTStimDataCorrect ./ tempDenom];
            PopShortRTStimDataError    = [PopShortRTStimDataError    ; tempShortRTStimDataError      ./ tempDenom];
            PopLongRTStimDataError     = [PopLongRTStimDataError; tempLongRTStimDataError       ./ tempDenom];
            
            
            
            
        end
         
    %average for each channel
    PopShortRTStimDataCorrect(chan_count, :) = nanmean(PopShortRTStimDataCorrect(chan_count:end, :));
    PopLongRTStimDataCorrect(chan_count, :) = nanmean(PopLongRTStimDataCorrect(chan_count:end, :));
    PopShortRTStimDataError(chan_count, :) = nanmean(PopShortRTStimDataError(chan_count:end, :));
    PopLongRTStimDataError(chan_count, :) = nanmean(PopLongRTStimDataError(chan_count:end, :));
    
    PopShortRTStimDataCorrect(chan_count+1:end, :) = [];
    PopLongRTStimDataCorrect(chan_count+1:end, :) = [];
    PopShortRTStimDataError(chan_count+1:end, :) = [];
    PopLongRTStimDataError(chan_count+1:end, :) = [];
    
    
    end
   
    
end
%%

figure;

plots(1) = subplot(1, 2, 1); % correct

shadedErrorBar_Morgane(1:13100, smooth(nanmean(PopShortRTStimDataCorrect), smooth_factor), smooth(nanstd(PopShortRTStimDataCorrect)/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', SmallLargeColor(2,:), 'LineWidth', 2}, 'patchSaturation', 0.12);
hold on;
shadedErrorBar_Morgane(1:13100, smooth(nanmean(PopLongRTStimDataCorrect), smooth_factor), smooth(nanstd(PopLongRTStimDataCorrect)/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', SmallLargeColor(2,:), 'LineWidth', 2, 'LineStyle', '--'}, 'patchSaturation', 0.12);
title('correct')
legend('Fast', 'Slow')
        

plots(3) = subplot(1, 2, 2); % error
shadedErrorBar_Morgane(1:13100, smooth(nanmean(PopShortRTStimDataError), smooth_factor), smooth(nanstd(PopShortRTStimDataError)/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', ErrorCorrectColor(1,:), 'LineWidth', 2}, 'patchSaturation', 0.12);
hold on; 
shadedErrorBar_Morgane(1:13100, smooth(nanmean(PopLongRTStimDataError), smooth_factor), smooth(nanstd(PopLongRTStimDataError)/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', ErrorCorrectColor(1,:), 'LineWidth', 2, 'LineStyle', '--'}, 'patchSaturation', 0.12);
title('error')
legend('Fast', 'Slow')


    set(plots, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)], 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off')


function [IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors()
IpsiContraColor = [ 0 0 0
                    0.2 0.2 0.2
                    0.4 0.4 0.4
                    0.6 0.6 0.7
                    0.6 0.4 0.8
                    0.6 0.2 0.9
                    0.6 0 1
                    ];
                
                IpsiContraColor2 = [ 0 0 0
                    0.6 0 1
                    ];
                
ErrorCorrectColor = [   'r'
                        'g'
                        ];
                    
SmallLargeColor = [ 0.15 0.75 0
                    0.15 0.4 0
                    ];

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


