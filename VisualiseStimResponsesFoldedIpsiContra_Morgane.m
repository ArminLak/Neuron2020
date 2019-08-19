
% August 2019 Morgane Moss
% This code compares responses to ipsi- and contralateral stimuli, 
% depending on stimulus contrast, choice accuracy, and value of pending 
% reward. 


%close all
clear all

% VTA,
 %region = 'VTA';
 %Animals = [48 50 51 64]
 %load('BehPhotoM_Exp23_VTA')

% NAC
% region = 'NAC';
% Animals = [56 57 59 66]
% load('BehPhotoM_Exp23_NAc')

%DMS
region = 'DMS';
Animals = [53, 62, 63, 71,72]
load('BehPhotoM_Exp23_DMS')

%%

[IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors();
totalChannels = getTotalChanN(region);
[yaxes] = getAxes(region);

StimAllowed = [-0.5 -0.25 -0.12 0 0.12 0.25 0.5];

GrandPopNormBinStimNoFold = zeros(2,7); %(
GrandPopNormBinStimNoFoldCorrError = zeros(2,7);

StimAbsResponsesIpsiContra = zeros(2,4, totalChannels); % separated by ipsi/contra. Ipsi is row 1; dim2 = contrast 0-0.5

StimResponsesErrorCorrectContra = zeros(2,4, totalChannels); % separated by err/corr. Error is row 1.
StimResponsesErrorCorrectIpsi = zeros(2,4, totalChannels);

StimResponsesRewardSizeContra = zeros(2,4, totalChannels); % separated by large/small reward. Small reward is row 1
StimResponsesRewardSizeIpsi = zeros(2, 4, totalChannels);
c = 1;
animal_count = 0;
for iAnimal = Animals
    
    animal_count = animal_count + 1;
    ChanN = 0;
    
    if ~isempty(BehPhotoM(iAnimal).GrandSummaryL)
        ChanN = ChanN + 1;
    end
    
    if ~isempty(BehPhotoM(iAnimal).GrandSummaryR)
        ChanN = ChanN + 1;
    end
    
    for iChan = 1:ChanN
        
        if isempty(BehPhotoM(iAnimal).GrandSummaryR)
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
            hem = 'l';
        elseif isempty(BehPhotoM(iAnimal).GrandSummaryL)
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
            hem = 'r';
        elseif ~isempty(BehPhotoM(iAnimal).GrandSummaryL) && ~isempty(BehPhotoM(iAnimal).GrandSummaryR)
            if iChan == 1
                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
                hem = 'l';
            elseif iChan==2
                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
                hem = 'r';
            end
        end
        
        % tuning curve for reward size
        SingleAnimalTunningStim= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimNoFold;
        SingleAnimalNormTunningStim = SingleAnimalTunningStim ./ max(max(SingleAnimalTunningStim));
        GrandPopNormBinStimNoFold = SingleAnimalNormTunningStim + GrandPopNormBinStimNoFold ;
        GrandPopNormBinStimNoFold1(c,:)=SingleAnimalNormTunningStim(1,:);
        GrandPopNormBinStimNoFold2(c,:)=SingleAnimalNormTunningStim(2,:);
        
        % tuning curve for corr/error
        SingleAnimalTunningStimCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrectErrorNoFold;
        SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError ./ max(max(SingleAnimalTunningStimCorrError));
        GrandPopNormBinStimNoFoldCorrError = SingleAnimalNormTunningStimCorrError + GrandPopNormBinStimNoFoldCorrError ;
        GrandPopNormBinStimNoFoldCorrError1(c,:)=SingleAnimalNormTunningStimCorrError(1,:);
        GrandPopNormBinStimNoFoldCorrError2(c,:)=SingleAnimalNormTunningStimCorrError(2,:);
        
        SingleAnimalTunningRewardCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardCorrectErrorNoFold;
        
        if abs(max(max(SingleAnimalTunningRewardCorrError))) < abs( min(min(SingleAnimalTunningRewardCorrError))) % this is to avoid diving by a small number in one animal
            
            SingleAnimalTunningRewardCorrError = SingleAnimalTunningRewardCorrError + abs(min(min(SingleAnimalTunningRewardCorrError)));
        end
        
        
        %average across blocks
        %     GrandPopNormBinStimNoFold = mean(GrandPopNormBinStimNoFold);
        
        if strcmpi(hem, 'l')
            
            StimAbsResponsesIpsiContra(1,:, c) = (StimAbsResponsesIpsiContra(1,:,c) + ...
                fliplr(mean(GrandPopNormBinStimNoFold(:,1:4))))/iChan; %ipsi = contrast -0.5 to 0
            
            StimAbsResponsesIpsiContra(2,:,c) = (StimAbsResponsesIpsiContra(2,:,c) + ...
                mean(GrandPopNormBinStimNoFold(:,4:7)))/iChan; %contra = contrast 0 to 0.5
            
            StimResponsesErrorCorrectIpsi(:,:, c) = (StimResponsesErrorCorrectContra(:, :, c) + ...
                fliplr(GrandPopNormBinStimNoFoldCorrError(:,1:4)))/iChan;
            
            StimResponsesErrorCorrectContra(:,:, c) = (StimResponsesErrorCorrectIpsi(:,:,c) + ...
                GrandPopNormBinStimNoFoldCorrError(:,4:7))/iChan;
            
            StimResponsesRewardSizeIpsi(2,:, c) = (StimResponsesRewardSizeIpsi(2,:, c) +...
                fliplr(GrandPopNormBinStimNoFold(1,1:4)))/iChan; %large reward + ipsi (left)
            
            StimResponsesRewardSizeContra(1,:, c) = (StimResponsesRewardSizeIpsi(1,:,c) + ...
                GrandPopNormBinStimNoFold(1,4:7))/iChan; % small reward and contra
            
            StimResponsesRewardSizeIpsi(1,:, c) = (StimResponsesRewardSizeIpsi(1,:, c) +...
                fliplr(GrandPopNormBinStimNoFold(2,1:4)))/iChan; %small reward + ipsi (left)
            
            StimResponsesRewardSizeContra(2,:, c) = (StimResponsesRewardSizeIpsi(2,:,c) + ...
                GrandPopNormBinStimNoFold(2,4:7))/iChan; % large reward + contra
            
            
        elseif strcmpi(hem, 'r')
            
            StimAbsResponsesIpsiContra(1,:,c) = (StimAbsResponsesIpsiContra(1,:,c) + ...
                mean(GrandPopNormBinStimNoFold(:,4:7)))/iChan; % ipsi = 0 to 0.5
            
            StimAbsResponsesIpsiContra(2,:,c) = (StimAbsResponsesIpsiContra(2,:,c) + ...
                fliplr(mean(GrandPopNormBinStimNoFold(:,1:4))))/iChan;  % contra = 0 to -0.5
            
            StimResponsesErrorCorrectContra(:,:, c) = (StimResponsesErrorCorrectContra(:, :, c) + ...
                fliplr(GrandPopNormBinStimNoFoldCorrError(:,1:4)))/iChan;
            
            StimResponsesErrorCorrectIpsi(:,:, c) = (StimResponsesErrorCorrectIpsi(:,:,c) + ...
                GrandPopNormBinStimNoFoldCorrError(:,4:7))/iChan;
            
            
            
            StimResponsesRewardSizeIpsi(1,:, c) = (StimResponsesRewardSizeIpsi(1,:, c) +...
                GrandPopNormBinStimNoFold(1,4:7))/iChan; % small reward + ipsi (right)
            
            StimResponsesRewardSizeContra(2,:, c) = (StimResponsesRewardSizeContra(2,:,c) + ...
                fliplr(GrandPopNormBinStimNoFold(1, 1:4)))/iChan; % large reward and contra (left)
            
            StimResponsesRewardSizeIpsi(2,:, c) = (StimResponsesRewardSizeIpsi(2,:, c) +...
                GrandPopNormBinStimNoFold(2,4:7))/iChan; % large reward + ipsi (right)
            
            StimResponsesRewardSizeContra(1,:, c) = (StimResponsesRewardSizeContra(1,:,c) + ...
                fliplr(GrandPopNormBinStimNoFold(2, 1:4)))/iChan; % small reward and contra (left)
            
        end
        
        c = c+1;
    end
    
end


StimAbsResponsesIpsiContra = mean(StimAbsResponsesIpsiContra,3);
StimResponsesErrorCorrectIpsi = mean(StimResponsesErrorCorrectIpsi,3);
StimResponsesErrorCorrectContra = mean(StimResponsesErrorCorrectContra, 3);
StimResponsesRewardSizeIpsi = mean(StimResponsesRewardSizeIpsi, 3);
StimResponsesRewardSizeContra = mean(StimResponsesRewardSizeContra, 3);


figure; 
subplot(3, 2, 1) % responses to absolute stim contrast, broken by ipsi and contra

for iHem = 1:2
    hold on;
    plot(StimAbsResponsesIpsiContra(iHem,:),'color', IpsiContraColor(iHem,:),'LineWidth',2,'Marker','o','MarkerSize',5)
    
end
title('All stim responses')
legend('Ipsi', 'Contra')

subplot(3, 2, 3) % ERROR VS CORRECT, IPSI ONLY

for iErrCorr = 1:2
    hold on;
    plot(StimResponsesErrorCorrectIpsi(iErrCorr,:),'color', ErrorCorrectColor(iErrCorr), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5)
end
ylim(yaxes)
title('Ipsi stim responses')
xticklabels(num2str(unique(abs(StimAllowed))))

subplot(3, 2, 4) % ERROR VS CORRECT, CONTRA ONLY

for iErrCorr = 1:2
    hold on;
    plot(StimResponsesErrorCorrectContra(iErrCorr,:),'color', ErrorCorrectColor(iErrCorr), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5)
end    
title('Contra stim responses')
ylim(yaxes)
subplot(3, 2, 5) % LARGE VS SMALL, IPSI ONLY

for iSmallLarge = 1:2
    hold on;
    plot(StimResponsesRewardSizeIpsi(iSmallLarge,:),'color', SmallLargeColor(iSmallLarge,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5)
end
ylim(yaxes)

subplot(3, 2, 6) % LARGE VS SMALL, CONTRA ONLY
for iSmallLarge = 1:2
    hold on;
    plot(StimResponsesRewardSizeContra(iSmallLarge,:),'color', SmallLargeColor(iSmallLarge,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5)
end
ylim(yaxes)


% ------- colors and functions
function [yaxes] = getAxes(region)

    if strcmpi(region, 'DMS')
        yaxes = [-1 4];
    elseif strcmpi(region, 'VTA')
        yaxes = [-5 5];
    elseif strcmpi(region, 'NAC')
        yaxes = [0 4];
    end

end

function [IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors()
IpsiContraColor = [ 0 0 0
                    153/255 51/255 1
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

