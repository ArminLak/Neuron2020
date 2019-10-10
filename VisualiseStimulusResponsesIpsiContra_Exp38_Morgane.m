% This code is for visualisation of stimulus- and action-aligned responses,
% separated by error and correct. 

% Morgane September 2019

% close all 
clear all

%DMS
Animals = [63, 71, 72];
load('BehPhotoM_Exp38_DMS')

%NAC
% Animals = [56 57 59 66];
% load('BehPhotoM_Exp23_NAc')

%VTA
% % Animals = [48 50 51 64];
% Animals = [51];
% load('BehPhotoM_Exp23_VTA')

AvgAcrossIpsiVsContra = 1; % make one line for ipsi and one line for contra, including all contrast levels? (ignores Stimz2plot)

Stimz2plot = 1:7; % or e.g. [1 7] (max level ipsi and contra) 
Stimz2plot = [1:6];

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 200; 

TimingVisualise = [-0.2 0.8
                   -0.7, 0.7
                   -0.2, 0.8]; % stim, action, reward in s
               
StartTime = 3700; % saved in the database.               
               

% --- get and organise data -----------------------------------------------
[IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors();

GrandPopStimLargeError     = zeros(length(Stimz2plot),13100);
GrandPopStimLargeCorrect   = zeros(length(Stimz2plot),13100);
GrandPopStimSmallCorrect   = zeros(length(Stimz2plot),13100);
GrandPopActionLargeCorrect = zeros(length(Stimz2plot),13100);
GrandPopActionSmallCorrect = zeros(length(Stimz2plot),13100);
GrandPopActionLargeError   = zeros(length(Stimz2plot),13100);

BehData_L               = [];
BehData_R               = [];
StimResponsesCorrect_L  = [];
StimResponsesCorrect_R  = [];

chan_count = 0;
for iAnimal = Animals
    
    
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
        
        for iSession = 1:length(BehPhotoM(iAnimal).Session)
            

            
            if iChan == 1
                            tempBehData = BehPhotoM(iAnimal).Session(iSession).TrialTimingData;
            BehData_L = [BehData_L; tempBehData(tempBehData(:,9)==1,:)];
            
                SingleAnimalStimCorrect       = BehPhotoM(iAnimal).Session(iSession).NeuronStimL(tempBehData(:,9)==1,:);
                SingleAnimalStimCorrect       = SingleAnimalStimCorrect ./ max(max(SingleAnimalStimCorrect));
                StimResponsesCorrect_L = [StimResponsesCorrect_L; SingleAnimalStimCorrect];
                
            elseif iChan == 2
                            tempBehData = BehPhotoM(iAnimal).Session(iSession).TrialTimingData;
            BehData_R = [BehData_R; tempBehData(tempBehData(:,9)==1,:)];
            
                SingleAnimalStimCorrect       = BehPhotoM(iAnimal).Session(iSession).NeuronStimR(tempBehData(:,9)==1,:);
                SingleAnimalStimCorrect       = SingleAnimalStimCorrect ./ max(max(SingleAnimalStimCorrect));
                StimResponsesCorrect_R = [StimResponsesCorrect_R; SingleAnimalStimCorrect];                
                
                
            end
            
%         SingleAnimalStimTraceLargeCorrect       = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeCorrect;
%         SingleAnimalStimTraceLargeCorrect       = SingleAnimalStimTraceLargeCorrect ./ max(max(SingleAnimalStimTraceLargeCorrect));
%         
%         SingleAnimalStimTraceLargeError         = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeError;
%         SingleAnimalStimTraceLargeError         = SingleAnimalStimTraceLargeError ./ max(max(SingleAnimalStimTraceLargeError));
% 
%         SingleAnimalActionTraceLargeCorrect     = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeCorrect;
%         SingleAnimalActionTraceLargeCorrect     = SingleAnimalActionTraceLargeCorrect ./ max(max(SingleAnimalActionTraceLargeCorrect));
%         
%         SingleAnimalActionTraceLargeError       = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeError;
%         SingleAnimalActionTraceLargeError       = SingleAnimalActionTraceLargeError ./ max(max(SingleAnimalActionTraceLargeError));
%         
%         if iChan == 1
%             GrandPopStimLargeCorrect(:,:,chan_count)    = SingleAnimalStimTraceLargeCorrect;
%             GrandPopStimLargeError(:,:,chan_count)      = SingleAnimalStimTraceLargeError;
%             GrandPopActionLargeCorrect(:,:,chan_count)  = SingleAnimalActionTraceLargeCorrect;
%             GrandPopActionLargeError(:,:,chan_count)    = SingleAnimalActionTraceLargeError;
% 
%         elseif iChan == 2
%                 
%             GrandPopStimLargeCorrect(:,:,chan_count)    = flipud(SingleAnimalStimTraceLargeCorrect);
%             GrandPopStimLargeError(:,:,chan_count)      = flipud(SingleAnimalStimTraceLargeError);
%             GrandPopActionLargeCorrect(:,:,chan_count)  = flipud(SingleAnimalActionTraceLargeCorrect);
%             GrandPopActionLargeError(:,:,chan_count)    = flipud(SingleAnimalActionTraceLargeError);            
%                
%             end
            
        end
    end   
end

% --- get averages for each stim level contra/ipsi  -----------------------

stimz = sort(unique([BehData_L(:,2) ; BehData_R(:,2)]));
c = 1
while c < 4
    
    StimResponsesIpsiContra(c,:) = mean( [StimResponsesCorrect_L(BehData_L(:,2)==stimz(c),:); StimResponsesCorrect_R(BehData_R(:,2)==-stimz(c),:)] );
    c = c+1;
end

while c < 7
    
    StimResponsesIpsiContra(c,:) = mean( [StimResponsesCorrect_L(BehData_L(:,2)==stimz(c),:); StimResponsesCorrect_R(BehData_R(:,2)==-stimz(c),:)] );
    c = c+1;
end

%%
% --- plotting ------------------------------------------------------------

figure; 

for c = 1: length(stimz)
    plot(smooth(StimResponsesIpsiContra(c,:), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
    hold on;
end

legend({'0.5 Ipsi', '0.25', '0.12', '0.12', '0.25', '0.5 Contra'})

 set(gca, 'ylim', [-0.1 0.2], 'xlim', [StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)], ...
        'YTick', [0 0.2], 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)], 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off')


% --- script-specific functions -------------------------------------------

function [IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors()
IpsiContraColor = [ 0 0 0
                    0.2 0.2 0.2
                    0.4 0.4 0.4
%                     0.6 0.6 0.7
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

