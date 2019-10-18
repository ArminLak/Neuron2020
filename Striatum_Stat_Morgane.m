% Morgane OCtober 2019 
% All tests for significance in DMS paper is done here. 

% to use this code:
% run section 1 and then run whichever section you want according to which tests you want to carry out


% CONTENTS:
% Section 2: STIM signals. ANOVAs on ipsi vs conta, and large vs correct, and correct vs error
% Section 3: ACTION signals. 
% Section 4: REWARD signals.


%% section 1: load and organise data 

clear all
close all


brain_region = 'DMS';
Animals = [53, 62, 63, 71,72];
load('BehPhotoM_Exp23_DMS')


%  load('BehPhotoM_Exp23_NAc')
%  Animals = [56 57 59 66];
%  load('MiceExpInfoPhotoM.mat')



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
        SingleAnimalNormTuningRewardCorrError = SingleAnimalTuningRewardCorrError ./ SingleAnimalNormTuningRewardDenom;
        SingleAnimalNormTuningReward          = SingleAnimalTuningReward ./ SingleAnimalNormTuningRewardDenom;
        
        
        if iChan == 1 % left hem
            
            GrandPopNormBinIpsiStim_Contrast(chan_count, :)   =  fliplr(SingleAnimalNormTuningStim(1,1:4));
            GrandPopNormBinContraStim_Contrast(chan_count, :) =  SingleAnimalNormTuningStim(2,4:7);
            
            GrandPopNormBinContraStim_SmallReward(chan_count, : ) = SingleAnimalNormTuningStim(1,4:7);
            GrandPopNormBinContraStim_LargeReward(chan_count, :) = SingleAnimalNormTuningStim(2,4:7);
            GrandPopNormBinContraStim_Error(chan_count, :) = SingleAnimalNormTuningStimCorrError(1, 4:7);
            
            
            GrandPopNormBinIpsiAction_Contrast(chan_count, :)   =  fliplr(SingleAnimalNormTuningAction(1,1:4)); % ipsi stimulus
            GrandPopNormBinContraAction_Contrast(chan_count, :) =  SingleAnimalNormTuningAction(2,4:7); % contra stimulus 
            
            GrandPopNormBinContraAction_SmallReward(chan_count, : ) = SingleAnimalNormTuningAction(1,4:7);
            GrandPopNormBinContraAction_LargeReward(chan_count, :) = SingleAnimalNormTuningAction(2,4:7);
            GrandPopNormBinContraAction_Error(chan_count, :) = SingleAnimalNormTuningActionCorrError(1, 4:7);
            
            
        elseif iChan == 2 % right hem
            
            % stim 
            GrandPopNormBinIpsiStim_Contrast(chan_count,: )  =  SingleAnimalNormTuningStim(2,4:7);
            GrandPopNormBinContraStim_Contrast(chan_count, :) =  fliplr(SingleAnimalNormTuningStim(1,1:4));
            
            GrandPopNormBinContraStim_SmallReward(chan_count, : ) = fliplr(SingleAnimalNormTuningStim(2, 1:4));
            GrandPopNormBinContraStim_LargeReward(chan_count, :) = fliplr(SingleAnimalNormTuningStim(1,1:4));
            GrandPopNormBinContraStim_Error(chan_count, :) = fliplr(SingleAnimalNormTuningStimCorrError(1, 1:4));
            
            % action 
            GrandPopNormBinIpsiAction_Contrast(chan_count, :)   =  SingleAnimalNormTuningAction(2,4:7); % ipsi stimulus
            GrandPopNormBinContraAction_Contrast(chan_count, :) =  fliplr(SingleAnimalNormTuningAction(1,1:4)); % contra stimulus 
            
            GrandPopNormBinContraAction_SmallReward(chan_count, : ) = fliplr(SingleAnimalNormTuningAction(2,1:4));
            GrandPopNormBinContraAction_LargeReward(chan_count, :) = fliplr(SingleAnimalNormTuningAction(1,1:4));
            GrandPopNormBinContraAction_Error(chan_count, :) = fliplr(SingleAnimalNormTuningActionCorrError(1, 1:4));
            
            % reward 
            GrandPopNormBinIpsiReward_Contrast(chan_count, :)   =  SingleAnimalNormTuningReward(2,4:7); % ipsi stimulus
            GrandPopNormBinContraReward_Contrast(chan_count, :) =  fliplr(SingleAnimalNormTuningReward(1,1:4)); % contra stimulus 
            
            GrandPopNormBinContraReward_SmallReward(chan_count, : ) = fliplr(SingleAnimalNormTuningReward(2,1:4));
            GrandPopNormBinContraReward_LargeReward(chan_count, :) = fliplr(SingleAnimalNormTuningReward(1,1:4));
            GrandPopNormBinContraReward_Error(chan_count, :) = fliplr(SingleAnimalNormTuningRewardCorrError(1, 1:4));
            
        end
        
        
    end
    
end

stimz = [0 0.12 0.25 0.5];

%% section 2.1 : linear regression to test for ipsi stimulus responses 

for i = 1:4
    contrasts(1:chan_count, i) = stimz(i);
end

contrasts = contrasts(:);
responses = GrandPopNormBinIpsiStim_Contrast(:);
ipsiStimResponseslm = fitlm(contrasts, responses);

%% section 2.2 : two way ANOVA to test for contrast on ipsi/contra stim responses

stimResponsesContrast = [GrandPopNormBinIpsiStim_Contrast; GrandPopNormBinContraStim_Contrast]; % first 7 rows = ipsi, next 7 = contra (7 channels)
[~,~,stimContrastStats] = anova2(stimResponsesContrast,chan_count);

% columns are contrast levels
% rows are ipsi / contra

%% section 2.3 : 2-way ANOVA to test reward size on contra stim responses

stimResponsesRewardSize = [GrandPopNormBinContraStim_LargeReward; GrandPopNormBinContraStim_SmallReward];
[~,~,stimRewardSizeStats] = anova2(stimResponsesRewardSize,chan_count);
        
%% section 2.4 : 2-way ANOVA to test choice accuracy on contra stim responses

stimResponsesChoiceAccuracy = [GrandPopNormBinContraStim_LargeReward(:, 1:3); GrandPopNormBinContraStim_Error(:, 1:3)];
[~,~,stimChoiceAccuracyStats] = anova2(stimResponsesChoiceAccuracy,chan_count);
        

%% section 3.1 : one-way ANOVA to test for action responses in hemisphere ipsi to stimulus 

[p, tbl, ipsiActionStats] = anova1(GrandPopNormBinIpsiAction_Contrast);

%% section 3.2 : two-way ANOVA to test for contrast on action responses in ipsi/contra (excluding zero contrast) 

actionResponsesContrast = [GrandPopNormBinIpsiAction_Contrast; GrandPopNormBinContraAction_Contrast]; % first n rows = ipsi, next n = contra (n = no. channels)
[~,~,actionContrastStats] = anova2(actionResponsesContrast,chan_count);

%% section 3.3 : 2-way ANOVA to test action wrt reward size 

actionResponsesRewardSize = [GrandPopNormBinContraAction_LargeReward(:,2:4); GrandPopNormBinContraAction_SmallReward(:,2:4)];
[~,~,actionRewardSizeStats] = anova2(actionResponsesRewardSize,chan_count);


%% section 4.1 : 2-way ANOVA to test reward response wrt contrast and  ipsi vs contra 

rewardResponsesContrast = [GrandPopNormBinIpsiReward_Contrast; GrandPopNormBinContraReward_Contrast]; % first n rows = ipsi, next n = contra (n = no. channels)
[~,~,stimContrastStats] = anova2(rewardResponsesContrast,chan_count);







        
        
        
        
        