% Morgane OCtober 2019 
% All tests for significance in DMS paper is done here. 

% to use this code:
% run section 1 and then run whichever section you want according to which tests you want to carry out


% CONTENTS:
% Section 2: STIM responses. ANOVAs on ipsi vs conta, and large vs correct, and correct vs error
% Section 3: ACTION responses. 



%% section 1: load and organise data 

clear all
close all


brain_region = 'DMS';
Animals = [53, 62, 63, 71,72];
load('BehPhotoM_Exp23_DMS')
sample_rate = 12000;   % photoM recording sampling rate
downSample = 1200;     % factor downsampling the Ca responses
eventOnset = 3700;     % in the saved matrix, the event time 
RTLimit = 6;           % remove very long trials 


GrandPopNormBinStimNoFold = zeros(2,7);
GrandPopNormBinStimErrCorrNoFold = zeros(2,7);
GrandPopNormBinActionNoFold = zeros(2,7);
GrandPopNormBinActionErrCorrNoFold = zeros(2,7);

chan_count = 0;
for iAnimal = Animals
    
    if isempty(BehPhotoM(iAnimal).GrandSummaryR)
        Chan = 1;
    elseif isempty(BehPhotoM(iAnimal).GrandSummaryL)
        Chan = 2;
    else
        Chan = [1 2];
    end
    
    for iChan = Chan
        
        chan_count = chan_count + 1;
        
        if iChan == 1
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
        elseif iChan == 2
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
        end
        
        % bin stim
        SingleAnimalTuningStim                 = BehPhotoM(iAnimal).GrandSummary.PopNormBinStimNoFold;
        SingleAnimalTuningStimCorrError        = BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrectErrorNoFold;
        
        SingleAnimalNormTuningStimDenom        = max( [ max(SingleAnimalTuningStim), max(SingleAnimalTuningStimCorrError)]);
        SingleAnimalNormTuningStimCorrError   = SingleAnimalTuningStimCorrError ./ SingleAnimalNormTuningStimDenom;
        SingleAnimalNormTuningStim            = SingleAnimalTuningStim ./ SingleAnimalNormTuningStimDenom;
        
        % bin action
        SingleAnimalTuningAction               = BehPhotoM(iAnimal).GrandSummary.PopNormBinActionNoFold;
        SingleAnimalTuningActionCorrError      = BehPhotoM(iAnimal).GrandSummary.PopNormBinActionCorrectErrorNoFold;
        
        SingleAnimalNormTuningActionDenom      = max( [ max(SingleAnimalTuningAction), max(SingleAnimalTuningActionCorrError)]);
        SingleAnimalNormTuningAction           = SingleAnimalTuningAction ./ SingleAnimalNormTuningActionDenom;
        SingleAnimalNormTuningActionCorrError  = SingleAnimalTuningActionCorrError ./ SingleAnimalNormTuningActionDenom;
        
        % bin action, stim = 0
        SingleAnimalTuningActionZero = BehPhotoM(iAnimal).GrandSummary.PopNormBinActionZeroLargeSmallErrorNoFold;
        SingleAnimalNormTuningActionZero = SingleAnimalTuningActionZero ./ max(max(SingleAnimalTuningActionZero));
        
        
        % bin reward
        SingleAnimalTuningReward               = BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardNoFold;
        SingleAnimalTuningRewardCorrError      = BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardCorrectErrorNoFold;
        
        SingleAnimalNormTuningRewardDenom      = max( [ max(SingleAnimalTuningReward), max(SingleAnimalTuningRewardCorrError)]);
        SingleAnimalNormTunningRewardCorrError = SingleAnimalTuningRewardCorrError ./ SingleAnimalNormTuningRewardDenom;
        SingleAnimalNormTunningReward          = SingleAnimalTuningReward ./ SingleAnimalNormTuningRewardDenom;
        
        
        if iChan == 1 % left hem
            
            GrandPopNormBinStimNoFold = GrandPopNormBinStimNoFold + SingleAnimalNormTuningStim;
            GrandPopNormBinStimErrCorrNoFold = GrandPopNormBinStimErrCorrNoFold + SingleAnimalNormTuningStimCorrError;
            GrandPopNormBinActionNoFold = GrandPopNormBinActionNoFold + SingleAnimalNormTuningAction;
            GrandPopNormBinActionErrCorrNoFold = GrandPopNormBinActionErrCorrNoFold + SingleAnimalNormTuningActionCorrError;
            
            GrandPopNormBinIpsiStim_Contrast(chan_count, :)   =  fliplr(SingleAnimalNormTuningStim(1,1:4));
            GrandPopNormBinContraStimContrast(chan_count, :) =  SingleAnimalNormTuningStim(2,4:7);
            
            GrandPopNormBinContraStim_SmallReward(chan_count, : ) = SingleAnimalNormTuningStim(1,4:7);
            GrandPopNormBinContraStim_LargeReward(chan_count, :) = SingleAnimalNormTuningStim(2,4:7);
            GrandPopNormBinContraStim_Error(chan_count, :) = SingleAnimalNormTuningStimCorrError(1, 4:7);
            
            
            GrandPopNormBinIpsiAction_Contrast(chan_count, :)   =  fliplr(SingleAnimalNormTuningAction(1,1:4)); % ipsi stimulus
            GrandPopNormBinContraActionContrast(chan_count, :) =  SingleAnimalNormTuningAction(2,4:7); % contra stimulus 
            
            GrandPopNormBinContraAction_SmallReward(chan_count, : ) = SingleAnimalNormTuningAction(1,4:7);
            GrandPopNormBinContraAction_LargeReward(chan_count, :) = SingleAnimalNormTuningAction(2,4:7);
            GrandPopNormBinContraAction_Error(chan_count, :) = SingleAnimalNormTuningActionCorrError(1, 4:7);
            
            
        elseif iChan == 2 % right hem
            
            GrandPopNormBinStimNoFold        = GrandPopNormBinStimNoFold + rot90(SingleAnimalNormTuningStim,2);
            GrandPopNormBinStimErrCorrNoFold = GrandPopNormBinStimErrCorrNoFold + fliplr(SingleAnimalNormTuningStimCorrError);
            GrandPopNormBinActionNoFold = GrandPopNormBinActionNoFold + rot90(SingleAnimalNormTuningAction, 2);
            GrandPopNormBinActionErrCorrNoFold = GrandPopNormBinActionErrCorrNoFold + fliplr(SingleAnimalNormTuningActionCorrError);
            
            GrandPopNormBinIpsiStimContrast(chan_count, :)   =  SingleAnimalNormTuningStim(2,4:7);
            GrandPopNormBinContraStimContrast(chan_count, :) =  fliplr(SingleAnimalNormTuningStim(1,1:4));
            
            GrandPopNormBinContraStim_SmallReward(chan_count, : ) = fliplr(SingleAnimalNormTuningStim(2, 1:4));
            GrandPopNormBinContraStim_LargeReward(chan_count, :) = fliplr(SingleAnimalNormTuningStim(1,1:4));
            GrandPopNormBinContraStim_Error(chan_count, :) = fliplr(SingleAnimalNormTuningStimCorrError(1, 1:4));
            
            
            GrandPopNormBinIpsiAction_Contrast(chan_count, :)   =  SingleAnimalNormTuningAction(2,4:7); % ipsi stimulus
            GrandPopNormBinContraActionContrast(chan_count, :) =  fliplr(SingleAnimalNormTuningAction(1,1:4)); % contra stimulus 
            
            GrandPopNormBinContraAction_SmallReward(chan_count, : ) = fliplr(SingleAnimalNormTuningAction(2,1:4));
            GrandPopNormBinContraAction_LargeReward(chan_count, :) = fliplr(SingleAnimalNormTuningAction(1,1:4));
            GrandPopNormBinContraAction_Error(chan_count, :) = fliplr(SingleAnimalNormTuningActionCorrError(1, 1:4));
            
        end
        
        
    end
    
end

%% section 2.1 : one-way ANOVA to test for ipsi stimulus responses 

[p, tbl, ipsiStimStats] = anova1(GrandPopNormBinIpsiStimContrast);

% continuous?

%% section 2.2 : two way ANOVA to test for contrast on ipsi/contra stim responses

stimResponsesContrast = [GrandPopNormBinIpsiStimContrast; GrandPopNormBinContraStimContrast]; % first 7 rows = ipsi, next 7 = contra (7 channels)
[~,~,stimContrastStats] = anova2(stimResponsesContrast,7);

% columns are contrast levels
% rows are ipsi / contra

%% section 2.3 : 2-way ANOVA to test reward size on contra stim responses

stimResponsesRewardSize = [GrandPopNormBinContraStim_LargeReward; GrandPopNormBinContraStim_SmallReward];
[~,~,stimRewardSizeStats] = anova2(stimResponsesRewardSize,7);
        
%% section 2.4 : 2-way ANOVA to test choice accuracy on contra stim responses

stimResponsesChoiceAccuracy = [GrandPopNormBinContraStim_LargeReward; GrandPopNormBinContraStim_Error];
[~,~,stimChoiceAccuracyStats] = anova2(stimResponsesChoiceAccuracy,7);
        

%% section 3.1 : one-way ANOVA to test for action responses in hemisphere ipsi to stimulus 

[p, tbl, ipsiActionStats] = anova1(GrandPopNormBinIpsiActionContrast);

%% section 3.2 : two-way ANOVA to test for contrast on action responses in ipsi/contra (excluding zero contrast) 

actionResponsesContrast = [GrandPopNormBinIpsiActionContrast; GrandPopNormBinContraActionContrast]; % first 7 rows = ipsi, next 7 = contra (7 channels)
[~,~,actionContrastStats] = anova2(actionResponsesContrast,7);

%% section 3.3 : 2-way ANOVA to test reward size 

actionResponsesRewardSize = [GrandPopNormBinContraAction_LargeReward(:,2:4); GrandPopNormBinContraAction_SmallReward(:,2:4)];
[~,~,actionRewardSizeStats] = anova2(actionResponsesRewardSize,7);






        
        
        
        
        