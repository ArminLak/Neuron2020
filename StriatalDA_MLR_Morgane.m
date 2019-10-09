% one regression on all trials of all animals

 if exist('brain_region', 'var') && strcmp(brain_region, 'DMS')
     clearvars -except BehPhotoM
 else clear all
      load('BehPhotoM_Exp23_DMS.mat')
 end
Animals = [53, 62, 63, 71,72]

 brain_region = 'DMS'
 
% NAC
%  if exist('brain_region', 'var') && strcmp(brain_region, 'NAc')
%      clearvars -except BehPhotoM
%   else clear all
%      load('BehPhotoM_Exp23_NAc.mat')
%  end
% Animals = [56 57 59 66]
% brain_region = 'NAc'

% VTA: 
%  if exist('brain_region', 'var') && strcmp(brain_region, 'VTA')
%      clearvars -except BehPhotoM
%  else clear all
%       load('BehPhotoM_Exp23_VTA.mat')
%  end
% Animals = [48 50 51 64]
% brain_region = 'VTA'
 


 RTLimit = 6;
 
 c = 0;
 
 for animal = Animals
 
     c = c + 1;

if isempty(BehPhotoM(animal).GrandSummaryR)
    Chan = 1;
elseif isempty(BehPhotoM(animal).GrandSummaryL)
    Chan = 2;
else
    Chan = [1 2];
end
    
BehData = [];
StimData = [];
ActionData = [];
RewardData = [];

StimDataL   = [];
ActionDataL = [];
RewardDataL = [];

StimDataR   = [];
ActionDataR = [];
RewardDataR = [];


isContra = [];
NormBinStim = [];
NormBinAction = [];
Contrast = [];
CorrErr = [];
RewardSize = [];
Stimz = [];
    
    sessionz = 1:length(BehPhotoM(animal).Session);
    
    iter = 0;
    if isfield(BehPhotoM(animal).Session,'NeuronRewardL')
        
        for iSession = sessionz % left hem
            
            if ~isfield(BehPhotoM(animal).Session, 'NeuronRewardR')
            TempBehData = BehPhotoM(animal).Session(iSession).TrialTimingData;
            BehData = [BehData; TempBehData];
            end
                    
            TempStimData= BehPhotoM(animal).Session(iSession).NeuronStimL;            
            TempActionData= BehPhotoM(animal).Session(iSession).NeuronActionL;            
            TempRewardData= BehPhotoM(animal).Session(iSession).NeuronRewardL;            
      
            StimDataL = [StimDataL;TempStimData];            
            ActionDataL = [ActionDataL;TempActionData];            
            RewardDataL = [RewardDataL;TempRewardData];           
        end
    end
    
    
    if isfield(BehPhotoM(animal).Session,'NeuronRewardR')
        for iSession = sessionz
            
            TempBehData = BehPhotoM(animal).Session(iSession).TrialTimingData;
            BehData = [BehData; TempBehData];

            % right           
            TempStimData= BehPhotoM(animal).Session(iSession).NeuronStimR;
            TempActionData= BehPhotoM(animal).Session(iSession).NeuronActionR;           
            TempRewardData= BehPhotoM(animal).Session(iSession).NeuronRewardR;            
        
            StimDataR = [StimDataR;TempStimData];           
            ActionDataR = [ActionDataR;TempActionData];            
            RewardDataR = [RewardDataR;TempRewardData];
                        
        end
        
    end
   

    RT = BehData(:,10) - BehData(:,13);
    toRemove1 = find ( RT > RTLimit);
    toRemove2= find(BehData(:,1) < 20);
    toRemove = unique([toRemove1; toRemove2]);
    BehData(toRemove,:) = [];

    
for iChan = Chan

    if iChan == 1
        StimData = StimDataL;
        ActionData = ActionDataL;
        RewardData = RewardDataL;

        contraTrials = find(BehData(:,9)==1 & BehData(:,8)==2 & BehData(:,3)==1);
        

            isContra = [isContra; double((BehData(:,9)==1 & BehData(:,3)==1) | (BehData(:,9)==0 & BehData(:,3)==-1))];


    elseif iChan == 2
        StimData = StimDataR;
        ActionData = ActionDataR;
        RewardData = RewardDataR;


            isContra = [isContra; double((BehData(:,9)==1 & BehData(:,3)==-1) | (BehData(:,9)==0 & BehData(:,3)==1))];

        
    end
    
    StimData(toRemove,:) = [];
    ActionData(toRemove,:) = [];
    RewardData(toRemove,:) = [];
    
    
    
    [NormBinStim, NormBinAction] = getResponses(animal, NormBinStim, NormBinAction, StimData, ActionData); % see function at bottom
    
    
    
    Contrast = [Contrast; abs(BehData(:,2))];
    CorrErr = [CorrErr; BehData(:,9)];
    RewardSize = [RewardSize; double((BehData(:,8)==2 & ((BehData(:,9)==1 & BehData(:,3)==1) | (BehData(:,9)==0 & BehData(:,3)==-1))) | ...
                                     (BehData(:,8)==1 & ((BehData(:,9)==1 & BehData(:,3)==-1) | (BehData(:,9)==0 & BehData(:,3)==1))))+1]; 
    
    
end


%isContra(isContra==0)=-1;
%CorrErr(CorrErr==0)=-1;

RewardSize(CorrErr == 0) = 0;
RewardSize (RewardSize ==1) =0.5;
RewardSize (RewardSize ==2) =1;

Contrast = Contrast * 2;

%  StimRegTbl = horzcat(isContra, Contrast, CorrErr, RewardSize);
  StimRegTbl = horzcat(isContra, Contrast, RewardSize);
  
  stimreg(:,c) = regress(NormBinStim, StimRegTbl);




 end

stimregmean = mean(stimreg,2);
figure; barh(stimregmean)
yticklabels({'Contralateral', 'Contrast', 'Accuracy', 'Reward size'})



function [NormBinStim, NormBinAction] = getResponses(animal, NormBinStim, NormBinAction, StimData, ActionData)

if animal == 48
        NormBinStim = [NormBinStim; mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2)- mean(StimData(:,3400:3800),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3000:3700),2) - mean(ActionData(:,3700:4300),2)];
        
    elseif animal == 50
        NormBinStim = [NormBinStim; mean(StimData(:,4500:5000),2) - mean(StimData(:,3100:3400),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3000:3700),2) - mean(ActionData(:,3700:4300),2)];
        
    elseif animal == 51
        
        NormBinStim = [NormBinStim; mean(StimData(:,4500:5000),2) - mean(StimData(:,3400:3800),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3000:3700),2) - mean(ActionData(:,3700:4300),2)];
        
    elseif animal == 53

        NormBinStim = [NormBinStim; mean(StimData(:,4000:4500),2) - mean(StimData(:,3400:3800),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3700:4400),2) - mean(ActionData(:,2900:3600),2)];
        
                elseif animal == 64
        
        NormBinStim = [NormBinStim;mean(StimData(:,4500:5000),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3000:3700),2) - mean(ActionData(:,3700:4300),2)];
        
        
    elseif animal == 62
        
        NormBinStim = [NormBinStim; mean(StimData(:, 4300:5000),2) - mean(StimData(:,3400:3700),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,4000:4300),2) - mean(ActionData(:,3300:3700),2)];
        
    elseif animal == 63
        
        NormBinStim = [NormBinStim; mean(StimData(:, 4200:4800),2) - mean(StimData(:,3400:3700),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3900:4200),2) - mean(ActionData(:,3500:3700),2)];
        
            elseif animal == 56
        
        NormBinStim = [NormBinStim;mean(StimData(:,4600:5300),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3900:4300),2) - mean(ActionData(:,3500:3700),2)];
        
    elseif animal == 57
        
        NormBinStim = [NormBinStim;mean(StimData(:,4600:5300),2)- mean(StimData(:,3400:3800),2)];
                NormBinAction = [NormBinAction;mean(ActionData(:,3900:4300),2) - mean(ActionData(:,3500:3700),2)];

    elseif animal == 59
        
        NormBinStim = [NormBinStim; mean(StimData(:,4900:5400),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3900:4300),2) - mean(ActionData(:,3500:3700),2)];
        
    elseif animal == 71
        
        NormBinStim = [NormBinStim; mean(StimData(:, 4200:4500),2) - mean(StimData(:,3600:3700),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,4000:4500),2) - mean(ActionData(:,3500:3700),2)];
        
            elseif animal == 66
        
        NormBinStim = [NormBinStim; mean(StimData(:,5000:6000),2)];
        NormBinAction = [NormBinAction;mean(ActionData(:,3900:4300),2) - mean(ActionData(:,3500:3700),2)];
        
    elseif animal == 72
        
        NormBinStim = [NormBinStim; mean(StimData(:,4100:4500),2) - mean(StimData(:,3650:3750),2)];
        NormBinAction = [NormBinAction; mean(ActionData(:,3900:4300),2) - mean(ActionData(:,3500:3700),2)];    
        
    

end
    
end

