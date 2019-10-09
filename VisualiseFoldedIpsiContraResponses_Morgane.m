% August 2019 Morgane Moss
% This code compares responses to ipsi- and contralateral stimuli, 
% depending on stimulus contrast, choice accuracy, and value of pending 
% reward. 


% close all
clear all

LargeSmallErrorOnSinglePlot = 0; %show large correct, small correct, and error  on same axes 
% VTA,
%  region = 'VTA';
%  Animals = [48 50 51 64]
%  load('BehPhotoM_Exp23_VTA')

% NAC
%   region = 'NAC';
%   Animals = [56 57 59 66]
%   load('BehPhotoM_Exp23_NAc')

% %DMS
  region = 'DMS';
Animals = [53, 62, 63, 71,72];
load('BehPhotoM_Exp23_DMS')


%%

[IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors();
totalChannels = getTotalChanN(region);
% totalChannels = 1;
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


GrandPopNormBinActionNoFold_L = zeros(2,7); %(
GrandPopNormBinActionNoFoldCorrError_L = zeros(2,7);
GrandPopNormBinActionNoFold1_L = zeros(totalChannels, 7);
GrandPopNormBinActionNoFold1_R = zeros(totalChannels, 7);
GrandPopNormBinActionNoFold2_L = zeros(totalChannels, 7);
GrandPopNormBinActionNoFold2_R = zeros(totalChannels, 7);
GrandPopNormBinActionNoFoldCorrError1_L = zeros(totalChannels, 7);
GrandPopNormBinActionNoFoldCorrError2_L = zeros(totalChannels, 7);
GrandPopNormBinActionNoFoldCorrError1_R = zeros(totalChannels, 7);
GrandPopNormBinActionNoFoldCorrError2_R = zeros(totalChannels, 7);
GrandPopNormBinActionNoFold_R = zeros(2,7); %(
GrandPopNormBinActionNoFoldCorrError_R = zeros(2,7);


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
         SingleAnimalTunningStim= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimNoFold;
       % SingleAnimalNormTunningStim = SingleAnimalTunningStim - min(min(SingleAnimalTunningStim));
       SingleAnimalNormTunningStim =  SingleAnimalTunningStim;
       SingleAnimalNormTunningStim = SingleAnimalNormTunningStim ./ max(max(SingleAnimalNormTunningStim));

        SingleAnimalTunningStimCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrectErrorNoFold;
        
        %SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError - min(min(SingleAnimalTunningStimCorrError));
        SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError;
        SingleAnimalNormTunningStimCorrError = SingleAnimalNormTunningStimCorrError ./ max(max(SingleAnimalNormTunningStimCorrError));
        
        SingleAnimalTuningAction = BehPhotoM(iAnimal).GrandSummary.PopNormBinActionNoFold;
        %SingleAnimalNomTuningAction = SingleAnimalTuningAction - min(min(SingleAnimalTuningAction));
         SingleAnimalNomTuningAction = SingleAnimalTuningAction ; 

        SingleAnimalNomTuningAction = SingleAnimalNomTuningAction ./ max(max(SingleAnimalNomTuningAction));
        
        
        SingleAnimalTuningActionCorrError = BehPhotoM(iAnimal).GrandSummary.PopNormBinActionCorrectErrorNoFold; 
%         SingleAnimalNormTuningActionCorrError = SingleAnimalTuningActionCorrError - min(min(SingleAnimalTuningActionCorrError));
         SingleAnimalNormTuningActionCorrError = SingleAnimalTuningActionCorrError;

SingleAnimalNormTuningActionCorrError = SingleAnimalNormTuningActionCorrError ./ max(max(SingleAnimalNormTuningActionCorrError));
%         
        
        SingleAnimalTuningActionZero = BehPhotoM(iAnimal).GrandSummary.PopNormBinActionZeroLargeSmallErrorNoFold;
        SingleAnimalNormTuningActionZero = SingleAnimalTuningActionZero - min(min(SingleAnimalTuningActionZero));
        SingleAnimalNormTuningActionZero = SingleAnimalNormTuningActionZero ./ max(max(SingleAnimalNormTuningActionZero)); % row 1 = ipsi, row 2 = contra. columns are large, small, error
         
        if strcmp(hem, 'l')
            GrandPopNormBinStimNoFold_L = GrandPopNormBinStimNoFold_L + SingleAnimalNormTunningStim; %+ GrandPopNormBinStimNoFold ;
            GrandPopNormBinStimNoFold1_L(c,:)=SingleAnimalNormTunningStim(1,:);
            GrandPopNormBinStimNoFold2_L(c,:)=SingleAnimalNormTunningStim(2,:);
            GrandPopNormBinStimNoFoldCorrError_L = GrandPopNormBinStimNoFoldCorrError_L + SingleAnimalNormTunningStimCorrError; %+ GrandPopNormBinStimNoFoldCorrError ;
            GrandPopNormBinStimNoFoldCorrError1_L(c,:)=SingleAnimalNormTunningStimCorrError(1,:);
            GrandPopNormBinStimNoFoldCorrError2_L(c,:)=SingleAnimalNormTunningStimCorrError(2,:);  
            
            GrandPopNormBinActionNoFold_L = GrandPopNormBinActionNoFold_L + SingleAnimalNomTuningAction;
            GrandPopNormBinActionNoFold1_L(c,:) = SingleAnimalNomTuningAction(1,:);
            GrandPopNormBinActionNoFold2_L(c,:) = SingleAnimalNomTuningAction(2,:);
            GrandPopNormBinActionNoFoldCorrError_L = GrandPopNormBinActionNoFoldCorrError_L + SingleAnimalNormTuningActionCorrError; %+ GrandPopNormBinStimNoFoldCorrError ;
            GrandPopNormBinActionNoFoldCorrError1_L(c,:) = SingleAnimalNormTuningActionCorrError(2,:);
            GrandPopNormBinActionNoFoldCorrError2_L(c,:) = SingleAnimalNormTuningActionCorrError(1,:);
            
        elseif strcmp(hem, 'r')
            GrandPopNormBinStimNoFold_R = GrandPopNormBinStimNoFold_R + SingleAnimalNormTunningStim; %+ GrandPopNormBinStimNoFold ;
            GrandPopNormBinStimNoFold1_R(c,:)=SingleAnimalNormTunningStim(2,:);
            GrandPopNormBinStimNoFold2_R(c,:)=SingleAnimalNormTunningStim(1,:);
            GrandPopNormBinStimNoFoldCorrError_R = GrandPopNormBinStimNoFoldCorrError_R + SingleAnimalNormTunningStimCorrError; %+ GrandPopNormBinStimNoFoldCorrError ;
            GrandPopNormBinStimNoFoldCorrError1_R(c,:)=SingleAnimalNormTunningStimCorrError(2,:);
            GrandPopNormBinStimNoFoldCorrError2_R(c,:)=SingleAnimalNormTunningStimCorrError(1,:);
            
            GrandPopNormBinActionNoFold_R = GrandPopNormBinActionNoFold_R + SingleAnimalNomTuningAction;
            GrandPopNormBinActionNoFold1_R(c,:) = SingleAnimalNomTuningAction(2,:);
            GrandPopNormBinActionNoFold2_R(c,:) = SingleAnimalNomTuningAction(1,:);
            GrandPopNormBinActionNoFoldCorrError_R = GrandPopNormBinActionNoFoldCorrError_R + SingleAnimalNormTuningActionCorrError; %+ GrandPopNormBinStimNoFoldCorrError ;
            GrandPopNormBinActionNoFoldCorrError1_R(c,:) = SingleAnimalNormTuningActionCorrError(2,:);
            GrandPopNormBinActionNoFoldCorrError2_R(c,:) = SingleAnimalNormTuningActionCorrError(1,:);
            

        end
            GrandPopNormActionZero_Ipsi(c,:) = SingleAnimalTuningActionZero(1,:);
            GrandPopNormActionZero_Contra(c,:) = SingleAnimalTuningActionZero(2,:);

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
    
    GrandPopNormBinActionNoFold = rot90(GrandPopNormBinActionNoFold_R,2) + GrandPopNormBinActionNoFold_L;
    GrandPopNormBinActionNoFold = GrandPopNormBinActionNoFold/totalChannels;
    
    GrandPopNormBinActionNoFold1 = GrandPopNormBinActionNoFold1_L + fliplr(GrandPopNormBinActionNoFold2_R);
    GrandPopNormBinActionNoFold2 = GrandPopNormBinActionNoFold2_L + fliplr(GrandPopNormBinActionNoFold1_R);
    
    GrandPopNormBinActionErrCorrNoFold = fliplr(GrandPopNormBinActionNoFoldCorrError_R) + GrandPopNormBinActionNoFoldCorrError_L;
    GrandPopNormBinActionErrCorrNoFold = GrandPopNormBinActionErrCorrNoFold/totalChannels;
    
    GrandPopNormBinActionNoFoldCorrError1 = GrandPopNormBinActionNoFoldCorrError1_L + fliplr(GrandPopNormBinActionNoFoldCorrError1_R);
    GrandPopNormBinActionNoFoldCorrError2 = GrandPopNormBinActionNoFoldCorrError2_L + fliplr(GrandPopNormBinActionNoFoldCorrError2_R);
    
    
    
    % ----- plotting 
    
figure; 
    
plots(1) = subplot(3, 2, 1); % responses to absolute stim contrast, broken by ipsi and contra

    errorbar(fliplr(GrandPopNormBinStimNoFold(1,1:4)),nanstd([GrandPopNormBinStimNoFold1(:,1:4);GrandPopNormBinStimNoFold2(:,1:4)]) ./ sqrt(totalChannels), ...
        'o', 'color', IpsiContraColor(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor', IpsiContraColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi; large correct trials only. to average across all correct, average rows 1 and 2. 
    hold on;
    errorbar(GrandPopNormBinStimNoFold(2,4:7), nanstd([GrandPopNormBinStimNoFold2(:,4:7);GrandPopNormBinStimNoFold1(:,4:7)])./sqrt(totalChannels), ...
        'o', 'color', IpsiContraColor(2,:),'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', IpsiContraColor(2,:), 'LineWidth',2, 'CapSize', 4) % large contra 
    title('All stim responses')
    legend('Ipsi', 'Contra')

plots(2) = subplot(3, 2, 3); % ERROR VS CORRECT, IPSI ONLY (col 1:4)
     errorbar(fliplr(GrandPopNormBinStimErrCorrNoFold(1,1:4)), nanstd(GrandPopNormBinStimNoFoldCorrError1(:,1:4)) ./ sqrt(totalChannels), ...
         'o', 'color', ErrorCorrectColor(1,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', ErrorCorrectColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi error
     hold on;
     if LargeSmallErrorOnSinglePlot
         errorbar(fliplr(GrandPopNormBinStimNoFold(2,1:4)), nanstd(GrandPopNormBinStimNoFold2(:,1:4)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi small
         hold on;
     end
     errorbar(fliplr(GrandPopNormBinStimNoFold(1,1:4)),nanstd([GrandPopNormBinStimNoFold1(:,1:4);GrandPopNormBinStimNoFold2(:,1:4)]) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(2,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %ipsi large
     title('Ipsi responses')
     xticks([])
     
plots(3) = subplot(3, 2, 4); % ERROR VS CORRECT, CONTRA ONLY (col 4:7)
     errorbar(GrandPopNormBinStimErrCorrNoFold(1,4:7), nanstd(GrandPopNormBinStimNoFoldCorrError1(:,4:7)) ./ sqrt(totalChannels), ...
         'o', 'color', ErrorCorrectColor(1,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', ErrorCorrectColor(1,:), 'LineWidth',2, 'CapSize', 4) %contra error
     hold on;
     if LargeSmallErrorOnSinglePlot
     errorbar(GrandPopNormBinStimNoFold(1,4:7), nanstd(GrandPopNormBinStimNoFold1(:,4:7)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %contra small
         hold on;
     end
     errorbar(GrandPopNormBinStimNoFold(2,4:7), nanstd([GrandPopNormBinStimNoFold2(:,4:7);GrandPopNormBinStimNoFold1(:,4:7)])./sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(2,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %contra large
     title('Contra responses')
     xticks([])

plots(4) = subplot(3, 2, 5); % LARGE VS SMALL, IPSI ONLY (col 1:4)
     errorbar(fliplr(GrandPopNormBinStimNoFold(2,1:4)), nanstd(GrandPopNormBinStimNoFold2(:,1:4)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi small
     hold on;
     errorbar(fliplr(GrandPopNormBinStimNoFold(1,1:4)),nanstd([GrandPopNormBinStimNoFold1(:,1:4);GrandPopNormBinStimNoFold2(:,1:4)]) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(2,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %ipsi large
     xlabel('|Contrast|')
     
plots(5) = subplot(3, 2, 6); % LARGE VS SMALL, CONTRA ONLY (col 4:7)
     errorbar(GrandPopNormBinStimNoFold(1,4:7), nanstd(GrandPopNormBinStimNoFold1(:,4:7)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %contra small
     hold on;
     errorbar(GrandPopNormBinStimNoFold(2,4:7), nanstd([GrandPopNormBinStimNoFold2(:,4:7);GrandPopNormBinStimNoFold1(:,4:7)])./sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(2,:), 'LineWidth',2,'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %contra large
     xlabel('|Contrast|')
     

%%     
figure; % action- responses ------------------------------------------------------------------------------------------------------------------
                 GrandPopNormActionZero1_Ipsi(c,:) = SingleAnimalTuningActionZero(1,:);
            GrandPopNormActionZero1_Contra(c,:) = SingleAnimalTuningActionZero(2,:);
            
plots(6) = subplot(3, 2, 1); % responses to absolute stim contrast, broken by ipsi and contra

    errorbar([2 3 4], fliplr(GrandPopNormBinActionNoFold(1,1:3)),nanstd([GrandPopNormBinActionNoFold1(:,1:3);GrandPopNormBinActionNoFold2(:,1:3)]) ./ sqrt(totalChannels), ...
        'o', 'color', IpsiContraColor(1,:),'MarkerSize',5, 'MarkerFaceColor', IpsiContraColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi; large correct trials only. to average across all correct, average rows 1 and 2. 
    hold on;
    errorbar([2 3 4], GrandPopNormBinActionNoFold(2,5:7), nanstd([GrandPopNormBinActionNoFold2(:,5:7);GrandPopNormBinActionNoFold1(:,5:7)])./sqrt(totalChannels), ...
        'o', 'color', IpsiContraColor(2,:),'MarkerSize',5, 'MarkerFaceColor', IpsiContraColor(2,:), 'LineWidth',2, 'CapSize', 4) % large contra 
    errorbar(1, nanmean(GrandPopNormActionZero_Ipsi(:,1)), nanstd(GrandPopNormActionZero_Ipsi(:,1))./sqrt(totalChannels), ...
        'o', 'color', IpsiContraColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', IpsiContraColor(1,:), 'LineWidth',2, 'CapSize', 4) % ipsi + large + 0
    errorbar(1, nanmean(GrandPopNormActionZero_Contra(:,1)), nanstd(GrandPopNormActionZero_Contra(:,1))./sqrt(totalChannels), ...
        'o', 'color', IpsiContraColor(2,:),'MarkerSize',5, 'MarkerFaceColor', IpsiContraColor(2,:), 'LineWidth',2, 'CapSize', 4) % contra + large + 0
    title('All action responses')
    legend('Ipsi', 'Contra')

plots(7) = subplot(3, 2, 3); % ERROR VS CORRECT, IPSI ONLY (col 1:4)
     errorbar([2 3 4], fliplr(GrandPopNormBinActionErrCorrNoFold(1,1:3)), nanstd(GrandPopNormBinActionNoFoldCorrError1(:,1:3)) ./ sqrt(totalChannels), ...
         'o', 'color', ErrorCorrectColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', ErrorCorrectColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi error
     hold on;
     errorbar(1, nanmean(GrandPopNormActionZero_Ipsi(:,3)), nanstd(GrandPopNormActionZero_Ipsi(:,3))./sqrt(totalChannels), ...
        'o', 'color', ErrorCorrectColor(1,:),'MarkerSize',5, 'MarkerFaceColor', ErrorCorrectColor(1,:), 'LineWidth',2, 'CapSize', 4) % ipsi + error + 0
     if LargeSmallErrorOnSinglePlot
         errorbar([2 3 4], fliplr(GrandPopNormBinActionNoFold(2,1:3)), nanstd(GrandPopNormBinActionNoFold2(:,1:3)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi small
         errorbar(1, nanmean(GrandPopNormActionZero_Ipsi(:,2)), nanstd(GrandPopNormActionZero_Ipsi(:,2))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) % ipsi + small + 0
         hold on;
     end
     errorbar([2 3 4], fliplr(GrandPopNormBinActionNoFold(1,1:3)),nanstd([GrandPopNormBinActionNoFold1(:,1:3);GrandPopNormBinActionNoFold2(:,1:3)]) ./ sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(2,:),'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %ipsi large
     errorbar(1, nanmean(GrandPopNormActionZero_Ipsi(:,1)), nanstd(GrandPopNormActionZero_Ipsi(:,1))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(2,:),'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) % ipsi + large + 0
     title('Ipsi responses')

     
plots(8) = subplot(3, 2, 4); % ERROR VS CORRECT, CONTRA ONLY (col 4:7)
     errorbar([2 3 4], GrandPopNormBinActionErrCorrNoFold(1,5:7), nanstd(GrandPopNormBinActionNoFoldCorrError1(:,5:7)) ./ sqrt(totalChannels), ...
         'o', 'color', ErrorCorrectColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', ErrorCorrectColor(1,:), 'LineWidth',2, 'CapSize', 4) %contra error
     hold on;
     errorbar(1, nanmean(GrandPopNormActionZero_Contra(:,3)), nanstd(GrandPopNormActionZero_Contra(:,3))./sqrt(totalChannels), ...
        'o', 'color', ErrorCorrectColor(1,:),'MarkerSize',5, 'MarkerFaceColor', ErrorCorrectColor(1,:), 'LineWidth',2, 'CapSize', 4) % contra + error + 0
     if LargeSmallErrorOnSinglePlot
        errorbar([2 3 4], GrandPopNormBinActionNoFold(1,5:7), nanstd(GrandPopNormBinActionNoFold1(:,5:7)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %contra small
         hold on;
         errorbar(1, nanmean(GrandPopNormActionZero_Contra(:,2)), nanstd(GrandPopNormActionZero_Contra(:,2))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) % contra + small + 0
         hold on;
     end
     errorbar([2 3 4], GrandPopNormBinActionNoFold(2,5:7), nanstd([GrandPopNormBinActionNoFold2(:,5:7);GrandPopNormBinActionNoFold1(:,5:7)])./sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(2,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %contra large
     errorbar(1, nanmean(GrandPopNormActionZero_Contra(:,1)), nanstd(GrandPopNormActionZero_Contra(:,1))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(2,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) % contra + large + 0
    
     title('Contra responses')

plots(9) = subplot(3, 2, 5); % LARGE VS SMALL, IPSI ONLY (col 1:4)
     errorbar([2 3 4], fliplr(GrandPopNormBinActionNoFold(2,1:3)), nanstd(GrandPopNormBinActionNoFold2(:,1:3)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %ipsi small
     hold on;
     errorbar(1, nanmean(GrandPopNormActionZero_Ipsi(:,2)), nanstd(GrandPopNormActionZero_Ipsi(:,2))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(1,:),'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) % ipsi + small + 0
     errorbar([2 3 4], fliplr(GrandPopNormBinActionNoFold(1,1:3)),nanstd([GrandPopNormBinActionNoFold1(:,1:3);GrandPopNormBinActionNoFold2(:,1:3)]) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(2,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %ipsi large
     errorbar(1, nanmean(GrandPopNormActionZero_Ipsi(:,1)), nanstd(GrandPopNormActionZero_Ipsi(:,1))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(2,:),'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) % ipsi + large + 0
     xlabel('|Contrast|')
     
plots(10) = subplot(3, 2, 6); % LARGE VS SMALL, CONTRA ONLY (col 4:7)
     errorbar([2 3 4], GrandPopNormBinActionNoFold(1,5:7), nanstd(GrandPopNormBinActionNoFold1(:,5:7)) ./ sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(1,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) %contra small
     hold on;
     errorbar(1, nanmean(GrandPopNormActionZero_Contra(:,2)), nanstd(GrandPopNormActionZero_Contra(:,2))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(1,:),'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(1,:), 'LineWidth',2, 'CapSize', 4) % contra + small + 0
     errorbar([2 3 4], GrandPopNormBinActionNoFold(2,5:7), nanstd([GrandPopNormBinActionNoFold2(:,5:7);GrandPopNormBinActionNoFold1(:,5:7)])./sqrt(totalChannels), ...
         'o', 'color', SmallLargeColor(2,:), 'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) %contra large
     errorbar(1, nanmean(GrandPopNormActionZero_Contra(:,1)), nanstd(GrandPopNormActionZero_Contra(:,1))./sqrt(totalChannels), ...
        'o', 'color', SmallLargeColor(2,:),'MarkerSize',5, 'MarkerFaceColor', SmallLargeColor(2,:), 'LineWidth',2, 'CapSize', 4) % contra + large + 0
     xlabel('|Contrast|')
     
     set(plots, 'YLim', ([0 1]))
     set(plots, 'XTick', ([]))
     set(plots([4 5 9 10]), 'XTick', ([1 2 3 4]), ...
           'XTickLabel', ([0 0.12 0.25 0.5]))
     

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
