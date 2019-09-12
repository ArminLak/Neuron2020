% This code is for visualisation of stimulus- and action-aligned responses,
% separated by error and correct. 

% Morgane September 2019: only deals with large reward trials 

% close all 
clear all

%DMS
% Animals = [53, 62, 63, 71,72];
% load('BehPhotoM_Exp23_DMS')

%VTA
Animals = [48 50 51 64]
load('BehPhotoM_Exp23_VTA')

AvgAcrossIpsiVsContra = 1; % make one line for ipsi and one line for contra, including all contrast levels? (ignores Stimz2plot)

Stimz2plot = 1:7; % or e.g. [1 7] (max level ipsi and contra) 
% Stimz2plot = [2 6];

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 200; 

TimingVisualise = [-0.2 0.8
                   -0.7, 0.7
                   -0.2, 0.8]; % stim, action, reward in s
               
StartTime = 3700; % saved in the database.               
               

% --- get and organise data -----------------------------------------------
[IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors();

GrandPopStimLargeError     = zeros(7,13100);
GrandPopStimLargeCorrect   = zeros(7,13100);
GrandPopActionLargeCorrect = zeros(7,13100);
GrandPopActionLargeError   = zeros(7,13100);

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
        BehPhotoM(iAnimal).GrandSummary     = [];
        SingleAnimalStimTraceLargeCorrect   = [];
        SingleAnimalStimTraceLargeError     = [];
        SingleAnimalActionTraceLargeCorrect = [];
        SingleAnimalActionTraceLargeError   = [];

        if iChan == 1
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
        elseif iChan == 2
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
        end
        
        SingleAnimalStimTraceLargeCorrect       = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeCorrect;
        SingleAnimalStimTraceLargeCorrect       = SingleAnimalStimTraceLargeCorrect ./ max(max(SingleAnimalStimTraceLargeCorrect));
        
        SingleAnimalStimTraceLargeError         = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeError;
        SingleAnimalStimTraceLargeError         = SingleAnimalStimTraceLargeError ./ max(max(SingleAnimalStimTraceLargeError));

        SingleAnimalActionTraceLargeCorrect     = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeCorrect;
        SingleAnimalActionTraceLargeCorrect     = SingleAnimalActionTraceLargeCorrect ./ max(max(SingleAnimalActionTraceLargeCorrect));
        
        SingleAnimalActionTraceLargeError       = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeError;
        SingleAnimalActionTraceLargeError       = SingleAnimalActionTraceLargeError ./ max(max(SingleAnimalActionTraceLargeError));
        
        if iChan == 1
            GrandPopStimLargeCorrect(:,:,chan_count)    = SingleAnimalStimTraceLargeCorrect;
            GrandPopStimLargeError(:,:,chan_count)      = SingleAnimalStimTraceLargeError;
            GrandPopActionLargeCorrect(:,:,chan_count)  = SingleAnimalActionTraceLargeCorrect;
            GrandPopActionLargeError(:,:,chan_count)    = SingleAnimalActionTraceLargeError;

        elseif iChan == 2
            GrandPopStimLargeCorrect(:,:,chan_count)    = flipud(SingleAnimalStimTraceLargeCorrect);
            GrandPopStimLargeError(:,:,chan_count)      = flipud(SingleAnimalStimTraceLargeError);
            GrandPopActionLargeCorrect(:,:,chan_count)  = flipud(SingleAnimalActionTraceLargeCorrect);
            GrandPopActionLargeError(:,:,chan_count)    = flipud(SingleAnimalActionTraceLargeError);            
               
        end
    end   
end

% --- divide by total number of channels ----------------------------------


AvgGrandPopStimLargeCorrect    = sum(GrandPopStimLargeCorrect,3) ./ chan_count;
AvgGrandPopStimLargeError      = sum(GrandPopStimLargeError,3) ./ chan_count;
AvgGrandPopActionLargeCorrect  = sum(GrandPopActionLargeCorrect,3) ./ chan_count;
AvgGrandPopActionLargeError    = sum(GrandPopActionLargeError,3) ./ chan_count;


%%
% --- plotting ------------------------------------------------------------
figure; 
% ROW 1: STIM RESPONSES 
stimplots(1) = subplot(2, 2, 1); % correct
for c = Stimz2plot
    hold on; 
    plot(smooth(AvgGrandPopStimLargeCorrect(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end
title('Large reward trials')
ylabel('Stimulus response')


stimplots(2) = subplot(2, 2, 2); % error 
for c = Stimz2plot
    hold on; 
    plot(smooth(AvgGrandPopStimLargeError(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end
title('Error trials')


actionplots(1) = subplot(2, 2, 3); % correct
for c = Stimz2plot
    hold on;
    plot(smooth(AvgGrandPopActionLargeCorrect(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end    
ylabel('Action response')


actionplots(2) = subplot(2, 2, 4); % correct
for c = Stimz2plot
    hold on;
    plot(smooth(AvgGrandPopActionLargeError(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end


    %%
figure; %figure 2: avg across all ipsi/contra 
plotindex = [1:3; 5:7];
% ROW 1: STIM RESPONSES 
stimplots(3) = subplot(2, 2, 1); % correct
for c = 1:2
    hold on; 
    shadedErrorBar(1:13100, smooth(mean(AvgGrandPopStimLargeCorrect(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopStimLargeCorrect(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
        'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.06)
end
title('Large reward trials')
ylabel('Stimulus response')

stimplots(4) = subplot(2, 2, 2); % error
for c = 1:2
    hold on; 
    shadedErrorBar(1:13100, smooth(mean(AvgGrandPopStimLargeError(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopStimLargeError(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
        'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.06)
end
title('Error trials')

actionplots(3) = subplot(2, 2, 3); % correct
for c = 1:2
    hold on;
    shadedErrorBar(1:13100, smooth(mean(AvgGrandPopActionLargeCorrect(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopActionLargeCorrect(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
        'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.06)
end    
ylabel('Action response')
xlabel('Time(s)')

actionplots(4) = subplot(2, 2, 4); % correct
for c = 1:2
    hold on;
    shadedErrorBar(1:13100, smooth(mean(AvgGrandPopActionLargeError(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopActionLargeError(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
        'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.06)
end
xlabel('Time (s)')

    set(stimplots, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)], ...
        'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)], 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off')
    set(actionplots, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTick', [StartTime + (TimingVisualise(2,1)*sampleRate), StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTickLabel', {'-0.7','0','0.7'},'TickDir','out','Box','off')

% --- script-specific functions -------------------------------------------

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
