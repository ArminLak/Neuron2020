% This code is for visualisation of stimulus- and action-aligned responses,
% separated by error and correct. 

% Morgane September 2019

% close all 
clear all

%DMS
% Animals = [53, 62, 63, 71,72];
% load('BehPhotoM_Exp23_DMS')

%NAC
Animals = [56 57 59 66];
load('BehPhotoM_Exp23_NAc')

%VTA
% % Animals = [48 50 51 64];
% Animals = [51];
% load('BehPhotoM_Exp23_VTA')

AvgAcrossIpsiVsContra = 1; % make one line for ipsi and one line for contra, including all contrast levels? (ignores Stimz2plot)

Stimz2plot = [1 7]; % 2 and 6 are 0.25 contrast; 1 and 7 are 0.5; 

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 200; 

TimingVisualise = [-0.2 0.8
                   -0.7, 0.7
                   -0.2, 0.8]; % stim, action, reward in s
               
StartTime = 3700; % saved in the database.


if strcmp(getComputerName(), 'proxeddu')
   addpath(genpath('D:\Morgane')); 
end

%% --- get and organise data -----------------------------------------------
[IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors();

GrandPopStimLargeError     = zeros(2,13100);
GrandPopStimLargeCorrect   = zeros(2,13100);
GrandPopStimSmallCorrect   = zeros(2,13100);

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

        if iChan == 1
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
        elseif iChan == 2
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
        end
        
        SingleAnimalStimTraceLargeCorrect       = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeCorrect(Stimz2plot,:);
%         SingleAnimalStimTraceLargeCorrect       = SingleAnimalStimTraceLargeCorrect ./ max(max(SingleAnimalStimTraceLargeCorrect));
        
        SingleAnimalStimTraceSmallCorrect       = BehPhotoM(iAnimal).GrandSummary.StimRasterSmallCorrect(Stimz2plot,:);
%         SingleAnimalStimTraceSmallCorrect       = SingleAnimalStimTraceSmallCorrect ./ max(max(SingleAnimalStimTraceSmallCorrect));
        
        SingleAnimalStimTraceLargeError         = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeError(Stimz2plot,:);
%         SingleAnimalStimTraceLargeError         = SingleAnimalStimTraceLargeError ./ max(max(SingleAnimalStimTraceLargeError));
        
                SingleAnimalNormStimDenom = max( [max(SingleAnimalStimTraceLargeCorrect), max(SingleAnimalStimTraceSmallCorrect), max(SingleAnimalStimTraceLargeError)] );
        SingleAnimalStimTraceLargeError         = SingleAnimalStimTraceLargeError ./ SingleAnimalNormStimDenom;                
        SingleAnimalStimTraceSmallCorrect       = SingleAnimalStimTraceSmallCorrect ./ SingleAnimalNormStimDenom;
        SingleAnimalStimTraceLargeCorrect       = SingleAnimalStimTraceLargeCorrect ./ SingleAnimalNormStimDenom;

        
        if iChan == 1
            GrandPopStimLargeCorrect(:,:,chan_count)    = SingleAnimalStimTraceLargeCorrect;
            GrandPopStimSmallCorrect(:,:,chan_count)    = SingleAnimalStimTraceSmallCorrect;
            GrandPopStimLargeError(:,:,chan_count)      = SingleAnimalStimTraceLargeError;

        elseif iChan == 2
            GrandPopStimLargeCorrect(:,:,chan_count)    = flipud(SingleAnimalStimTraceLargeCorrect);
            GrandPopStimSmallCorrect(:,:,chan_count)    = flipud(SingleAnimalStimTraceSmallCorrect);
            GrandPopStimLargeError(:,:,chan_count)      = flipud(SingleAnimalStimTraceLargeError);       
               
        end
    end   
end

AvgGrandPopStimLargeCorrect    = mean(GrandPopStimLargeCorrect, 3);
AvgGrandPopStimSmallCorrect    = mean(GrandPopStimSmallCorrect, 3);
AvgGrandPopStimLargeError      = mean(GrandPopStimLargeError, 3);


%% plots
figure; 

for i = 2 % 1 = ipsi, 2 = contra
    stimplots(1) = subplot(1,2,1);
    shadedErrorBar_Morgane(1:13100, smooth(AvgGrandPopStimLargeCorrect(i,:), smooth_factor), smooth(std(squeeze(GrandPopStimLargeCorrect(i,:,:))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', SmallLargeColor(2,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    hold on;     
        shadedErrorBar_Morgane(1:13100, smooth(AvgGrandPopStimSmallCorrect(i,:), smooth_factor), smooth(std(squeeze(GrandPopStimLargeCorrect(i,:,:))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', SmallLargeColor(1,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
        
    stimplots(2) = subplot(1,2,2);
    shadedErrorBar_Morgane(1:13100, smooth(AvgGrandPopStimLargeCorrect(i,:), smooth_factor), smooth(std(squeeze(GrandPopStimLargeCorrect(i,:,:))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', SmallLargeColor(2,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    hold on; 
        shadedErrorBar_Morgane(1:13100, smooth(AvgGrandPopStimLargeError(i,:), smooth_factor), smooth(std(squeeze(GrandPopStimLargeCorrect(i,:,:))')/sqrt(chan_count), smooth_factor), ...
            'lineprops', {'color', ErrorCorrectColor(1,:), 'LineWidth', 2}, 'patchSaturation', 0.12)
    
    hold on; 
        
end



    set(stimplots, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)], 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off')



% script functions: 

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
