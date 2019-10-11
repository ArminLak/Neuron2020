% Morgane October 2019

% close all 
clear all




%DMS
Animals = [53, 62, 63, 71,72];
load('BehPhotoM_Exp23_DMS')

%NAC
% Animals = [56 57 59 66];
% load('BehPhotoM_Exp23_NAc')

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 200; 

TimingVisualise = [-0.2 0.8
                   -0.7, 0.7
                   -0.2, 0.8]; % stim, action, reward in s
               
StartTime = 3700; % saved in the database.               
               


% --- get and organise data -----------------------------------------------
[IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors();

GrandPopRewardLargeCorrect = zeros(7,13100);
GrandPopRewardSmallCorrect = zeros(7,13100);
GrandPopRewardError   = zeros(7,13100);


figure; 


for brain_region = 1:2
    
    if brain_region == 1
        
        Animals = [53, 62, 63, 71,72];
        load('BehPhotoM_Exp23_DMS')
        
    elseif brain_region == 2
        Animals = [56 57 59 66];
        load('BehPhotoM_Exp23_NAc')
 
    end


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
        BehPhotoM(iAnimal).GrandSummary       = [];
        SingleAnimalRewardTraceLargeCorrect   = [];
        SingleAnimalRewardTraceSmallCorrect   = [];
        SingleAnimalRewardTraceError          = [];

        if iChan == 1
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
        elseif iChan == 2
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
        end
        
        SingleAnimalRewardTraceLargeCorrect     = BehPhotoM(iAnimal).GrandSummary.Rew2DACorr;
        SingleAnimalRewardTraceLargeCorrect     = SingleAnimalRewardTraceLargeCorrect ./ max(max(SingleAnimalRewardTraceLargeCorrect));
        
        
        SingleAnimalRewardTraceSmallCorrect     = BehPhotoM(iAnimal).GrandSummary.RewAwayDACorr;
        SingleAnimalRewardTraceSmallCorrect     = SingleAnimalRewardTraceSmallCorrect ./ max(max(SingleAnimalRewardTraceSmallCorrect));
        
        SingleAnimalRewardTraceError            = mean([BehPhotoM(iAnimal).GrandSummary.RewAwayDAErr; BehPhotoM(iAnimal).GrandSummary.RewAwayDACorr]);
        SingleAnimalRewardTraceError            = SingleAnimalRewardTraceError ./ max(max(SingleAnimalRewardTraceError));
        
        
        GrandPopRewardLargeCorrect(chan_count,:)    = SingleAnimalRewardTraceLargeCorrect;
        GrandPopRewardSmallCorrect(chan_count,:)    = SingleAnimalRewardTraceSmallCorrect;
        GrandPopRewardError(chan_count,:)      = SingleAnimalRewardTraceError;
        
    end
end

GrandPopRewardLargeCorrect  = mean(GrandPopRewardLargeCorrect);
GrandPopRewardSmallCorrect  = mean(GrandPopRewardSmallCorrect);
GrandPopRewardError         = mean(GrandPopRewardError);





%% plot

plots(brain_region) = subplot(1,2,brain_region); % first VS then DS

plot(smooth(GrandPopRewardError, smooth_factor), 'color', 'r', 'LineWidth', 2)
hold on;
plot(smooth(GrandPopRewardSmallCorrect, smooth_factor), 'color', SmallLargeColor(1,:), 'LineWidth', 2)
plot(smooth(GrandPopRewardLargeCorrect, smooth_factor), 'color', SmallLargeColor(2,:), 'LineWidth', 2)






end

set(plots, 'ylim', [-0.8 2.5], 'xlim', [StartTime + (TimingVisualise(3,1)*sampleRate) StartTime + (TimingVisualise(3,2)*sampleRate)], ...
        'YTick', [-0.8 0 1 2], 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)], 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off')   

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


