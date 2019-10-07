% This code is for visualisation of stimulus- and action-aligned responses,
% separated by error and correct. 

% Morgane September 2019

% close all 
clear all

%DMS
Animals = [53, 62, 63, 71,72];
load('BehPhotoM_Exp23_DMS')

%NAC
% Animals = [56 57 59 66];
% load('BehPhotoM_Exp23_NAc')

%VTA
% % Animals = [48 50 51 64];
% Animals = [51];
% load('BehPhotoM_Exp23_VTA')

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
GrandPopStimSmallCorrect   = zeros(7,13100);
GrandPopActionLargeCorrect = zeros(7,13100);
GrandPopActionSmallCorrect = zeros(7,13100);
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
        SingleAnimalStimTraceSmallCorrect   = [];
        SingleAnimalStimTraceLargeError     = [];
        SingleAnimalActionTraceLargeCorrect = [];
        SingleAnimalActionTraceSmallCorrect = [];
        SingleAnimalActionTraceLargeError   = [];

        if iChan == 1
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
        elseif iChan == 2
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
        end
        
        SingleAnimalStimTraceLargeCorrect       = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeCorrect;
        SingleAnimalStimTraceLargeCorrect       = SingleAnimalStimTraceLargeCorrect ./ max(max(SingleAnimalStimTraceLargeCorrect));
        
        SingleAnimalStimTraceSmallCorrect       = BehPhotoM(iAnimal).GrandSummary.StimRasterSmallCorrect;
        SingleAnimalStimTraceSmallCorrect       = SingleAnimalStimTraceSmallCorrect ./ max(max(SingleAnimalStimTraceSmallCorrect));
        
        SingleAnimalStimTraceLargeError         = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeError;
        SingleAnimalStimTraceLargeError         = SingleAnimalStimTraceLargeError ./ max(max(SingleAnimalStimTraceLargeError));

        SingleAnimalActionTraceLargeCorrect     = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeCorrect;
        SingleAnimalActionTraceLargeCorrect     = SingleAnimalActionTraceLargeCorrect ./ max(max(SingleAnimalActionTraceLargeCorrect));

        SingleAnimalActionTraceSmallCorrect     = BehPhotoM(iAnimal).GrandSummary.ActionRasterSmallCorrect;
        SingleAnimalActionTraceSmallCorrect     = SingleAnimalActionTraceSmallCorrect ./ max(max(SingleAnimalActionTraceSmallCorrect));
        
        SingleAnimalActionTraceLargeError       = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeError;
        SingleAnimalActionTraceLargeError       = SingleAnimalActionTraceLargeError ./ max(max(SingleAnimalActionTraceLargeError));
        
        if iChan == 1
            GrandPopStimLargeCorrect(:,:,chan_count)    = SingleAnimalStimTraceLargeCorrect;
            GrandPopStimSmallCorrect(:,:,chan_count)    = SingleAnimalStimTraceSmallCorrect;
            GrandPopStimLargeError(:,:,chan_count)      = SingleAnimalStimTraceLargeError;
            GrandPopActionLargeCorrect(:,:,chan_count)  = SingleAnimalActionTraceLargeCorrect;
            GrandPopActionSmallCorrect(:,:,chan_count)  = SingleAnimalActionTraceSmallCorrect;
            GrandPopActionLargeError(:,:,chan_count)    = SingleAnimalActionTraceLargeError;

        elseif iChan == 2
            GrandPopStimLargeCorrect(:,:,chan_count)    = flipud(SingleAnimalStimTraceLargeCorrect);
            GrandPopStimSmallCorrect(:,:,chan_count)    = flipud(SingleAnimalStimTraceSmallCorrect);
            GrandPopStimLargeError(:,:,chan_count)      = flipud(SingleAnimalStimTraceLargeError);
            GrandPopActionLargeCorrect(:,:,chan_count)  = flipud(SingleAnimalActionTraceLargeCorrect);
            GrandPopActionSmallCorrect(:,:,chan_count)  = flipud(SingleAnimalActionTraceSmallCorrect);
            GrandPopActionLargeError(:,:,chan_count)    = flipud(SingleAnimalActionTraceLargeError);            
               
        end
    end   
end

% --- divide by total number of channels ----------------------------------


AvgGrandPopStimLargeCorrect    = sum(GrandPopStimLargeCorrect,3) ./ chan_count;
AvgGrandPopStimLargeError      = sum(GrandPopStimLargeError,3) ./ chan_count;
AvgGrandPopStimSmallCorrect    = sum(GrandPopStimSmallCorrect,3) ./ chan_count;
AvgGrandPopActionLargeCorrect  = sum(GrandPopActionLargeCorrect,3) ./ chan_count;
AvgGrandPopActionSmallCorrect  = sum(GrandPopActionSmallCorrect,3) ./ chan_count;
AvgGrandPopActionLargeError    = sum(GrandPopActionLargeError,3) ./ chan_count;


%%
% --- plotting ------------------------------------------------------------
figure; 
% ROW 1: STIM RESPONSES 
stimplots(1) = subplot(3,2, 1); % large reward 
for c = Stimz2plot
    hold on; 
    plot(smooth(AvgGrandPopStimLargeCorrect(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end
title('Stimulus response')
ylabel('Large reward trials')

stimplots(2) = subplot(3,2, 3); % small reward
for c = Stimz2plot
    hold on; 
    plot(smooth(AvgGrandPopStimSmallCorrect(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end
ylabel('Small reward trials')

stimplots(3) = subplot(3,2, 5); % error 
for c = Stimz2plot
    hold on; 
    plot(smooth(AvgGrandPopStimLargeError(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end
ylabel('Error trials')


actionplots(1) = subplot(3,2, 2); % large reward
for c = Stimz2plot
    hold on;
    plot(smooth(AvgGrandPopActionLargeCorrect(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end    
title('Action response')


actionplots(2) = subplot(3,2, 4); % small reward
for c = Stimz2plot
    hold on;
    plot(smooth(AvgGrandPopActionSmallCorrect(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end  


actionplots(3) = subplot(3,2, 6); % error
for c = Stimz2plot
    hold on;
    plot(smooth(AvgGrandPopActionLargeError(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
end
xlabel('Time (s)')

    %%
    
if length(Animals) > 1
    
    figure; %figure 2: avg across all ipsi/contra
    plotindex = [1:3; 5:7];
    % ROW 1: STIM RESPONSES
    stimplots(4) = subplot(3,2, 1); % large reward
    for c = 1:2
        hold on;
        shadedErrorBar_Morgane(1:13100, smooth(mean(AvgGrandPopStimLargeCorrect(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopStimLargeCorrect(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    end
    title('Stimulus response')
    ylabel('Large reward trials')

    stimplots(5) = subplot(3,2, 3); % small reward
    for c = 1:2
        hold on;
        shadedErrorBar_Morgane(1:13100, smooth(mean(AvgGrandPopStimSmallCorrect(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopStimSmallCorrect(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    end
    ylabel('Small reward trials')
    
    stimplots(6) = subplot(3, 2, 5); % error
    for c = 1:2
        hold on;
        shadedErrorBar_Morgane(1:13100, smooth(mean(AvgGrandPopStimLargeError(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopStimLargeError(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    end
     xlabel('Time(s)')
         ylabel('Error trials')

    actionplots(4) = subplot(3, 2, 2); % large reward
    for c = 1:2
        hold on;
        shadedErrorBar_Morgane(1:13100, smooth(mean(AvgGrandPopActionLargeCorrect(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopActionLargeCorrect(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    end
    title('Action response')

    actionplots(5) = subplot(3, 2, 4); % small reward
    for c = 1:2
        hold on;
        shadedErrorBar_Morgane(1:13100, smooth(mean(AvgGrandPopActionSmallCorrect(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopActionSmallCorrect(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    end
   

    actionplots(6) = subplot(3, 2, 6); % error
    for c = 1:2
        hold on;
        shadedErrorBar_Morgane(1:13100, smooth(mean(AvgGrandPopActionLargeError(plotindex(c,:),:)), smooth_factor), smooth(std(squeeze(mean(GrandPopActionLargeError(plotindex(c,:),:,:)))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', IpsiContraColor2(c,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    end
    xlabel('Time (s)')

end


    set(stimplots, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)], 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off')
    set(actionplots, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime + (TimingVisualise(2,1)*sampleRate), StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTickLabel', {'-0.7','0','0.7'},'TickDir','out','Box','off')
    if length(Animals)>1
%     set([stimplots([1 4]), actionplots([1 4])], 'YTickLabel', {'-0.4', '0', '0.6'});
%     set(actionplots(4:6), 'ylim', [-0.6 0.8])
%     set(stimplots(4:6), 'ylim', [-0.6 0.8])
    else
    set([stimplots(1), actionplots(1)], 'YTickLabel', {'-0.4', '0', '0.6'});
    end

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

