% August 2019 Morgane Moss
% This code compares responses to ipsi- and contralateral stimuli, 
% depending on stimulus contrast, choice accuracy, and value of pending 
% reward. 


% close all
clear all

% VTA,
% 
 region = 'VTA';
 Animals = [48 50 51 64]
 load('BehPhotoM_Exp23_VTA')

% NAC
%  region = 'NAC';
%  Animals = [56 57 59 66]
%  load('BehPhotoM_Exp23_NAc')

%DMS
% region = 'DMS';
% Animals = [53, 62, 63, 71,72]
% load('BehPhotoM_Exp23_DMS')

%%

[IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors();
totalChannels = getTotalChanN(region);
[yaxes] = getAxes(region);

StimAllowed = [-0.5 -0.25 -0.12 0 0.12 0.25 0.5];

GrandPopNormBinStimNoFold_L = zeros(2,7); %(
GrandPopNormBinStimNoFoldCorrError_L = zeros(2,7);
GrandPopNormBinStimNoFold1_L = zeros(totalChannels, 7);
GrandPopNormBinStimNoFold1_R = zeros(totalChannels, 7);
GrandPopNormBinStimNoFold2_L = zeros(totalChannels, 7);
GrandPopNormBinStimNoFold2_R = zeros(totalChannels, 7);
GrandPopNormBinStimNoFoldCorrError1_L = zeros(totalChannels, 7);
GrandPopNormBinStimNoFoldCorrError2_L = zeros(totalChannels, 7);
GrandPopNormBinStimNoFoldCorrError1_R = zeros(totalChannels, 7);
GrandPopNormBinStimNoFoldCorrError2_R = zeros(totalChannels, 7);
GrandPopNormBinStimNoFold_R = zeros(2,7); %(
GrandPopNormBinStimNoFoldCorrError_R = zeros(2,7);

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
        
        
        % ---- subtract then divide ------
%         SingleAnimalTunningStim= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimNoFold;
%         SingleAnimalNormTunningStim = SingleAnimalTunningStim - min(min(SingleAnimalTunningStim));
%         SingleAnimalNormTunningStim = SingleAnimalNormTunningStim ./ max(max(SingleAnimalNormTunningStim));
%         
%         SingleAnimalTunningStimCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrectErrorNoFold;
%         SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError - min(min(SingleAnimalTunningStimCorrError));
%         SingleAnimalNormTunningStimCorrError = SingleAnimalNormTunningStimCorrError ./ max(max(SingleAnimalNormTunningStimCorrError));
        
        % ---- divide then subsract ------
        SingleAnimalTunningStim= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimNoFold;
        SingleAnimalNormTunningStim = SingleAnimalTunningStim ./ max(max(SingleAnimalTunningStim));
        SingleAnimalNormTunningStim = SingleAnimalTunningStim - min(min(SingleAnimalTunningStim));
        
        SingleAnimalTunningStimCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrectErrorNoFold;
        SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError ./ max(max(SingleAnimalTunningStimCorrError));
        SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError - min(min(SingleAnimalTunningStimCorrError));
               
        
         
        if strcmp(hem, 'l')
            GrandPopNormBinStimNoFold_L = GrandPopNormBinStimNoFold_L + SingleAnimalNormTunningStim; %+ GrandPopNormBinStimNoFold ;
            GrandPopNormBinStimNoFold1_L(c,:)=SingleAnimalNormTunningStim(1,:);
            GrandPopNormBinStimNoFold2_L(c,:)=SingleAnimalNormTunningStim(2,:);
            GrandPopNormBinStimNoFoldCorrError_L = GrandPopNormBinStimNoFoldCorrError_L + SingleAnimalNormTunningStimCorrError; %+ GrandPopNormBinStimNoFoldCorrError ;
            GrandPopNormBinStimNoFoldCorrError1_L(c,:)=SingleAnimalNormTunningStimCorrError(1,:);
            GrandPopNormBinStimNoFoldCorrError2_L(c,:)=SingleAnimalNormTunningStimCorrError(2,:);
            
        elseif strcmp(hem, 'r')
            GrandPopNormBinStimNoFold_R = GrandPopNormBinStimNoFold_R + SingleAnimalNormTunningStim; %+ GrandPopNormBinStimNoFold ;
            GrandPopNormBinStimNoFold1_R(c,:)=SingleAnimalNormTunningStim(2,:);
            GrandPopNormBinStimNoFold2_R(c,:)=SingleAnimalNormTunningStim(1,:);
            GrandPopNormBinStimNoFoldCorrError_R = GrandPopNormBinStimNoFoldCorrError_R + SingleAnimalNormTunningStimCorrError; %+ GrandPopNormBinStimNoFoldCorrError ;
            GrandPopNormBinStimNoFoldCorrError1_R(c,:)=SingleAnimalNormTunningStimCorrError(2,:);
            GrandPopNormBinStimNoFoldCorrError2_R(c,:)=SingleAnimalNormTunningStimCorrError(1,:);
            
            
        end

        c = c+1;
    end
    
end
     
    GrandPopNormBinStimNoFold = rot90(GrandPopNormBinStimNoFold_R,2) + GrandPopNormBinStimNoFold_L;
    GrandPopNormBinStimNoFold = GrandPopNormBinStimNoFold/totalChannels;
    
    GrandPopNormBinStimNoFold1 = GrandPopNormBinStimNoFold1_L + fliplr(GrandPopNormBinStimNoFold2_R);
    GrandPopNormBinStimNoFold2 = GrandPopNormBinStimNoFold2_L + fliplr(GrandPopNormBinStimNoFold1_R);

    GrandPopNormBinStimErrCorrNoFold = fliplr(GrandPopNormBinStimNoFoldCorrError_R) + GrandPopNormBinStimNoFoldCorrError_L;
    GrandPopNormBinStimErrCorrNoFold = GrandPopNormBinStimErrCorrNoFold/totalChannels;
    
    GrandPopNormBinStimNoFoldCorrError1 = GrandPopNormBinStimNoFoldCorrError1_L + fliplr(GrandPopNormBinStimNoFoldCorrError1_R); %error trials
    GrandPopNormBinStimNoFoldCorrError2 = GrandPopNormBinStimNoFoldCorrError2_L + fliplr(GrandPopNormBinStimNoFoldCorrError2_R); % correct trials
    
    % ----- plotting 
    
figure; 
    
subplot(3, 2, 1) % responses to absolute stim contrast, broken by ipsi and contra

    errorbar(fliplr(GrandPopNormBinStimNoFold(1,1:4)),nanstd([GrandPopNormBinStimNoFold1(:,1:4);GrandPopNormBinStimNoFold2(:,1:4)]) ./ sqrt(totalChannels), ...
        'color', IpsiContraColor(1,:),'LineWidth',2,'Marker','o','MarkerSize',5) %ipsi; large correct trials only. to average across all correct, average rows 1 and 2. 
    hold on;
    errorbar(GrandPopNormBinStimNoFold(2,4:7), nanstd([GrandPopNormBinStimNoFold2(:,4:7);GrandPopNormBinStimNoFold1(:,4:7)])./sqrt(totalChannels), ...
        'color', IpsiContraColor(2,:),'LineWidth',2,'Marker','o','MarkerSize',5) % large contra 
    title('All stim responses')
    legend('Ipsi', 'Contra')

subplot(3, 2, 3) % ERROR VS CORRECT, IPSI ONLY (col 1:4)
     errorbar(fliplr(GrandPopNormBinStimErrCorrNoFold(1,1:4)), nanstd(GrandPopNormBinStimNoFoldCorrError1(:,1:4)) ./ sqrt(totalChannels), ...
         'color', ErrorCorrectColor(1,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %ipsi error
     hold on;
     errorbar(fliplr(GrandPopNormBinStimNoFold(1,1:4)),nanstd([GrandPopNormBinStimNoFold1(:,1:4);GrandPopNormBinStimNoFold2(:,1:4)]) ./ sqrt(totalChannels), ...
         'color', SmallLargeColor(2,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %ipsi large
     title('Ipsi responses')
     xticks([])
     
subplot(3, 2, 4) % ERROR VS CORRECT, CONTRA ONLY (col 4:7)
     errorbar(GrandPopNormBinStimErrCorrNoFold(1,4:7), nanstd(GrandPopNormBinStimNoFoldCorrError1(:,4:7)) ./ sqrt(totalChannels), ...
         'color', ErrorCorrectColor(1,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %contra small
     hold on;
     errorbar(GrandPopNormBinStimNoFold(2,4:7), nanstd([GrandPopNormBinStimNoFold2(:,4:7);GrandPopNormBinStimNoFold1(:,4:7)])./sqrt(totalChannels), ...
         'color', SmallLargeColor(2,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %contra large
     title('Contra responses')
     xticks([])

subplot(3, 2, 5) % LARGE VS SMALL, IPSI ONLY (col 1:4)
     errorbar(fliplr(GrandPopNormBinStimNoFold(2,1:4)), nanstd(GrandPopNormBinStimNoFold2(:,1:4)) ./ sqrt(totalChannels), ...
         'color', SmallLargeColor(1,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %ipsi small
     hold on;
     errorbar(fliplr(GrandPopNormBinStimNoFold(1,1:4)),nanstd([GrandPopNormBinStimNoFold1(:,1:4);GrandPopNormBinStimNoFold2(:,1:4)]) ./ sqrt(totalChannels), ...
         'color', SmallLargeColor(2,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %ipsi large
     xticks([1 2 3 4])
     xticklabels([0 0.12 0.25 0.5])
     xlabel('|Contrast|')
     
subplot(3, 2, 6) % LARGE VS SMALL, CONTRA ONLY (col 4:7)
     errorbar(GrandPopNormBinStimNoFold(1,4:7), nanstd(GrandPopNormBinStimNoFold1(:,4:7)) ./ sqrt(totalChannels), ...
         'color', SmallLargeColor(1,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %contra small
     hold on;
     errorbar(GrandPopNormBinStimNoFold(2,4:7), nanstd([GrandPopNormBinStimNoFold2(:,4:7);GrandPopNormBinStimNoFold1(:,4:7)])./sqrt(totalChannels), ...
         'color', SmallLargeColor(2,:), 'LineWidth',2,'Marker','o','MarkerSize',5) %contra large
     xticks([1 2 3 4])
     xticklabels([0 0.12 0.25 0.5])
     xlabel('|Contrast|')

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
