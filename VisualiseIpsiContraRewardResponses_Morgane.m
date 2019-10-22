
% Morgane October 2019 
% Stimulus-aligned responses for striatal terminals. 6 lines (ipsi/contra,
% error/small/large). 

close all
clear all

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 200; 

TimingVisualise = [-0.2 0.8
                   -0.7, 0.7
                   -0.2, 0.8]; % stim, action, reward in s
               
StartTime = 3700; % saved in the database.               
               

% --- get and organise data -----------------------------------------------
[Colors] = getColors();

figure; 

for brain_region = 1:2
    
    if brain_region == 1
        
        Animals = [53, 62, 63, 71,72];
        %         Animals = [63];
        load('BehPhotoM_Exp23_DMS')
        
    elseif brain_region == 2
        
        Animals = [56 57 59 66];
        load('BehPhotoM_Exp23_NAc')
    end
    
    L_chan_count = 0;
    R_chan_count = 0;
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
            
            %rasters :
            
            SingleAnimalIpsiContraRewardLargeCorrect     = BehPhotoM(iAnimal).GrandSummary.RewardRasterIpsiContraLargeCorr;
            SingleAnimalIpsiContraRewardLargeCorrect     = SingleAnimalIpsiContraRewardLargeCorrect - min(min(SingleAnimalIpsiContraRewardLargeCorrect));
            
            SingleAnimalIpsiContraRewardSmallCorrect     = BehPhotoM(iAnimal).GrandSummary.RewardRasterIpsiContraSmallCorr;
            SingleAnimalIpsiContraRewardSmallCorrect     = SingleAnimalIpsiContraRewardSmallCorrect - min(min(SingleAnimalIpsiContraRewardSmallCorrect));
            
            SingleAnimalIpsiContraRewardError            = BehPhotoM(iAnimal).GrandSummary.RewardRasterIpsiContraErr;
            SingleAnimalIpsiContraRewardError            = SingleAnimalIpsiContraRewardError - min(min(SingleAnimalIpsiContraRewardError));
            
            SingleAnimalNormDenom                        = max( [max(SingleAnimalIpsiContraRewardLargeCorrect), max(SingleAnimalIpsiContraRewardSmallCorrect), ...
                                                         max(SingleAnimalIpsiContraRewardError)] );
            
            SingleAnimalIpsiContraRewardLargeCorrect     = SingleAnimalIpsiContraRewardLargeCorrect ./ SingleAnimalNormDenom;
            SingleAnimalIpsiContraRewardSmallCorrect     = SingleAnimalIpsiContraRewardSmallCorrect ./ SingleAnimalNormDenom;
            SingleAnimalIpsiContraRewardError            = SingleAnimalIpsiContraRewardError ./ SingleAnimalNormDenom;

            GrandPopRewardIpsiLargeSmallError(1,:, chan_count)      = SingleAnimalIpsiContraRewardLargeCorrect(1,:);
            GrandPopRewardIpsiLargeSmallError(2,:, chan_count)      = SingleAnimalIpsiContraRewardSmallCorrect(1,:);
            GrandPopRewardIpsiLargeSmallError(3,:, chan_count)      = SingleAnimalIpsiContraRewardError(1,:);
            
            GrandPopRewardContraLargeSmallError(1,:, chan_count)    = SingleAnimalIpsiContraRewardLargeCorrect(2,:);
            GrandPopRewardContraLargeSmallError(2,:, chan_count)    = SingleAnimalIpsiContraRewardSmallCorrect(2,:);
            GrandPopRewardContraLargeSmallError(3,:, chan_count)    = SingleAnimalIpsiContraRewardError(2,:);
            
            
            %binned: 
            
            SingleAnimalNormBinReward = BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardCorrect;
            SingleAnimalNormBinError = BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardError;
            SingleAnimalNormBinReward = SingleAnimalNormBinReward - min([min(min(SingleAnimalNormBinReward)), min(min(SingleAnimalNormBinError))]);
            SingleAnimalNormBinError = SingleAnimalNormBinError - min([min(min(SingleAnimalNormBinError)), min(min(SingleAnimalNormBinError))]);
            SingleAnimalNormBinReward = SingleAnimalNormBinReward ./max([max(max(SingleAnimalNormBinReward)), max(max(SingleAnimalNormBinError))]);
            SingleAnimalNormBinError = SingleAnimalNormBinError ./ max([max(max(SingleAnimalNormBinReward)), max(max(SingleAnimalNormBinError))]);
            

            GrandPopNormBinRewardLargeSmallError(1,:, chan_count) = SingleAnimalNormBinReward(2,:); %contra large
            GrandPopNormBinRewardLargeSmallError(2,:, chan_count) = SingleAnimalNormBinReward(1,:); % contra small
            GrandPopNormBinRewardLargeSmallError(3,:, chan_count) = SingleAnimalNormBinError(2,:); % contra error (2largeRew) 
            
              
            end
 
        end

    
%     GrandPopNormBinContraRewardLargeSmallError = GrandPopNormBinContraRewardLargeSmallError ./ chan_count;
 
    GrandPopRewardIpsiLargeSmallError   = mean(GrandPopRewardIpsiLargeSmallError, 3);
    GrandPopRewardContraLargeSmallError = mean(GrandPopRewardContraLargeSmallError, 3); 
    
    plots(brain_region) = subplot(2,2,brain_region); % left = DS, right = VS
    iColor = 1;
    for i = 1:3
        plot(smooth(GrandPopRewardIpsiLargeSmallError(i,:),200), '--', 'LineWidth', 2, 'color', Colors(iColor,:));
        
        hold on;
        plot(smooth(GrandPopRewardContraLargeSmallError(i,:),200), 'LineWidth', 2, 'color', Colors(iColor,:));
        iColor = iColor + 1;
    end
    xlabel('Time (s)')
    ylabel('\Delta F/F')
    
    plots2(brain_region + 2) = subplot(2, 2, brain_region + 2);
    iColor2 = 1;
    for i = 1:3
        
%         tempForSTDipsi = squeeze(GrandPopNormBinIpsiRewardLargeSmallError(i,:,:));
%         tempForSTDcontra = squeeze(GrandPopNormBinContraRewardLargeSmallError(i,:,:));
        
        errorbar(mean(GrandPopNormBinRewardLargeSmallError(i,:,:), 3), std(squeeze(GrandPopNormBinRewardLargeSmallError(i,:,:))')./sqrt(chan_count),'-', 'LineWidth', 2, 'color', Colors(iColor2,:));
        hold on;
        iColor2 = iColor2 + 1;
        
    end
end

set(plots, 'ylim', [0 1], 'xlim', [TimingVisualise(3,1)*sampleRate+StartTime TimingVisualise(3,2)*sampleRate+StartTime], ...
        'YTick', [0 1], 'XTick', [StartTime TimingVisualise(3,2)*sampleRate+StartTime], 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off')   
set(plots2, 'XTick', [1 2 3 4], 'XTickLabel', {'0', '0.12', '0.25', '0.5'}, 'TickDir', 'out', 'YTick', [0 1])
    legend({'Large reward', 'Small reward', 'Error'})
    
    
function [Colors] = getColors()
Colors = [  0.15 0.4 0  %contra large
            0.15 0.75 0 %contra small
            1 0 0       %contra error
            ];

end



