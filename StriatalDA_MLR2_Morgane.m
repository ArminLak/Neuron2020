


% close all 
clear all


% DMS

%  if exist('BehPhotoM', 'var') && strcmp(brain_region, 'DMS')
%      clearvars -except BehPhotoM
%  else clear all
%       load('BehPhotoM_Exp23_DMS.mat')
% %  end
% Animals = [53, 62, 63, 71,72]

% NAC
%   if exist('BehPhotoM', 'var') && strcmp(brain_region, 'NAc')
%       clearvars -except BehPhotoM
%    else clear all
%       load('BehPhotoM_Exp23_NAc.mat')
% 
%  Animals = [56 57 59 66]


brain_region = 'DMS';
load('BehPhotoM_Exp23_DMS.mat')
Animals = [53, 62, 63, 71,72]

% NAC
% brain_region = 'NAC'
% load('BehPhotoM_Exp23_NAc.mat')
% Animals = [56 57 59 66]
% 

% VTA 
%  if exist('BehPhotoM', 'var') && strcmp(brain_region, 'VTA')
%      clearvars -except BehPhotoM
%  else clear all
%       load('BehPhotoM_Exp23_VTA.mat')
%  end
% Animals = [48 50 51 64]




 RTLimit = 6;
 



isContra = [];
NormBinStim = [];
NormBinAction = [];
NormBinReward = [];
Contrast = [];
CorrErr = [];
RewardSize = [];
Stimz = [];
BehData = [];
 
 cstim = 0;
 
 for animal = Animals
     
TempBehData = [];
TempStimData = [];
TempActionData = [];
TempRewardData = [];

     cstim = cstim + 1;

if isempty(BehPhotoM(animal).GrandSummaryR)
    Chan = 1;
elseif isempty(BehPhotoM(animal).GrandSummaryL)
    Chan = 2;
else
    Chan = [1 2];
end


sessionz = 1:length(BehPhotoM(animal).Session);

iter = 0;
if isfield(BehPhotoM(animal).Session,'NeuronRewardL')
    
    for iSession = sessionz % left hem
        
        TempBehData = BehPhotoM(animal).Session(iSession).TrialTimingData;
        
        RT = TempBehData(:,10) - TempBehData(:,13);
        toRemove1 = find ( RT > RTLimit);
        toRemove2= find(TempBehData(:,1) < 20);
        toRemove = unique([toRemove1; toRemove2]);
        TempBehData(toRemove,:) = [];
        
        
        
        
        BehData = [BehData; TempBehData];
        
        TempStimData= BehPhotoM(animal).Session(iSession).NeuronStimL;
        TempActionData= BehPhotoM(animal).Session(iSession).NeuronActionL;
        TempRewardData= BehPhotoM(animal).Session(iSession).NeuronRewardL;
        
        TempStimData(toRemove,:) = [];
        TempActionData(toRemove,:) = [];
        TempRewardData(toRemove,:) = [];
        
        contraTrials = find(TempBehData(:,9)==1 & TempBehData(:,8)==2 & TempBehData(:,3)==1);
        
        
        isContra = [isContra; double((TempBehData(:,9)==1 & TempBehData(:,3)==1) | (TempBehData(:,9)==0 & TempBehData(:,3)==-1))];
        [NormBinStim, NormBinAction, NormBinReward] = getResponses(animal, NormBinStim, NormBinAction, NormBinReward, TempStimData, TempActionData, TempRewardData); % see function at bottom
        
        
    end
end


if isfield(BehPhotoM(animal).Session,'NeuronRewardR')
    for iSession = sessionz
        
        
        TempBehData = BehPhotoM(animal).Session(iSession).TrialTimingData;
        RT = TempBehData(:,10) -TempBehData(:,13);
        toRemove1 = find ( RT > RTLimit);
        toRemove2= find(TempBehData(:,1) < 20);
        toRemove = unique([toRemove1; toRemove2]);
        TempBehData(toRemove,:) = [];
        
        BehData = [BehData; TempBehData];
        
        % right
        TempStimData= BehPhotoM(animal).Session(iSession).NeuronStimR;
        TempActionData= BehPhotoM(animal).Session(iSession).NeuronActionR;
        TempRewardData= BehPhotoM(animal).Session(iSession).NeuronRewardR;
        
        TempStimData(toRemove,:) = [];
        TempActionData(toRemove,:) = [];
        TempRewardData(toRemove,:) = [];
        
        
        isContra = [isContra; double((TempBehData(:,9)==1 & TempBehData(:,3)==-1) | (TempBehData(:,9)==0 & TempBehData(:,3)==1))];
        
        [NormBinStim, NormBinAction, NormBinReward] = getResponses(animal, NormBinStim, NormBinAction, NormBinReward, TempStimData, TempActionData, TempRewardData); % see function at bottom
        
    end
    
end
   
 end
 
 Contrast = abs(BehData(:,2));
 CorrErr = BehData(:,9);
 RewardSize = double((BehData(:,8)==2 & ((BehData(:,9)==1 & BehData(:,3)==1) | (BehData(:,9)==0 & BehData(:,3)==-1))) | ...
                                     (BehData(:,8)==1 & ((BehData(:,9)==1 & BehData(:,3)==-1) | (BehData(:,9)==0 & BehData(:,3)==1))))+1;
 
RewardSize(CorrErr == 0) = 0;
RewardSize (RewardSize ==1) =0.5;
RewardSize (RewardSize ==2) =1;

Contrast = Contrast * 2;
 
 

 RegTbl = horzcat(isContra, Contrast, RewardSize);
figure; barh(regress(NormBinStim, RegTbl))
yticklabels({'Contralateral', 'Contrast', 'Reward size'})


 [b,bint,r,rint,stats] = regress(NormBinStim, RegTbl,0.01);  % last number refers to alpha
 
%% Cross-Validation 5 fold 

cvEvalFunc = @(pred, actual)1-mean(mean((pred-actual).^2))/mean(mean(actual.^2));
cstim = cvpartition(NormBinStim, 'KFold', 5);
cact = cvpartition(NormBinAction, 'Kfold', 5);
crwd = cvpartition(NormBinReward, 'Kfold', 5);

for i = 1:5
    
    [cvbstim, cvbintstim]= regress(NormBinStim(test(cstim,i)), RegTbl(test(cstim,i),:), 0.01);
    pred = cvbstim(1)*RegTbl(test(cstim,i),1) + cvbstim(2)*RegTbl(test(cstim,i),2) + cvbstim(3)*RegTbl(test(cstim,i),3);
    actual = NormBinStim(test(cstim,i));
    
    cvEvalStim(i) = cvEvalFunc(pred,actual);
    
    [cvbact, cvbintact]= regress(NormBinAction(test(cact,i)), RegTbl(test(cact,i),:), 0.01);
    pred = cvbact(1)*RegTbl(test(cact,i),1) + cvbact(2)*RegTbl(test(cact,i),2) + cvbact(3)*RegTbl(test(cact,i),3);
    actual = NormBinAction(test(cact,i));
    
    cvEvalAct(i) = cvEvalFunc(pred,actual);
    
    [cvbrwd, cvbintrwd]= regress(NormBinReward(test(crwd,i)), RegTbl(test(crwd,i),:), 0.01);
    pred = cvbrwd(1)*RegTbl(test(crwd,i),1) + cvbrwd(2)*RegTbl(test(crwd,i),2) + cvbrwd(3)*RegTbl(test(crwd,i),3);
    actual = NormBinReward(test(crwd,i));
    
    cvEvalRwd(i) = cvEvalFunc(pred,actual);
    
    
end




for i = 1:5     % testing w/o each beta for stim regression  
    actual = NormBinStim(test(cstim,i));
    
    [cvbstim, cvbintstim]= regress(NormBinStim(test(cstim,i)), RegTbl(test(cstim,i),[2 3]), 0.01);
    pred = cvbstim(1)*RegTbl(test(cstim,i),2) + cvbstim(2)*RegTbl(test(cstim,i),3);
    cvEvalStim_noIpsiContra(i) = cvEvalFunc(pred,actual);
    
    [cvbstim, cvbintstim]= regress(NormBinStim(test(cstim,i)), RegTbl(test(cstim,i),[1 3]), 0.01);
    pred = cvbstim(1)*RegTbl(test(cstim,i),1) + cvbstim(2)*RegTbl(test(cstim,i),3);
    cvEvalStim_noContrast(i) = cvEvalFunc(pred,actual);
    
    [cvbstim, cvbintstim]= regress(NormBinStim(test(cstim,i)), RegTbl(test(cstim,i),[1 2]), 0.01);
    pred = cvbstim(1)*RegTbl(test(cstim,i),1) + cvbstim(2)*RegTbl(test(cstim,i),2);
    cvEvalStim_noRewardSize(i) = cvEvalFunc(pred,actual);
    
    [cvbstim, cvbintstim]= regress(NormBinStim(test(cstim,i)), RegTbl(test(cstim,i),1), 0.01);
    pred = cvbstim*RegTbl(test(cstim,i),1);
    cvEvalStim_onlyIpsiContra(i) = cvEvalFunc(pred,actual);
    
    [cvbstim, cvbintstim]= regress(NormBinStim(test(cstim,i)), RegTbl(test(cstim,i),2), 0.01);
    pred = cvbstim*RegTbl(test(cstim,i),2);
    cvEvalStim_onlyContrast(i) = cvEvalFunc(pred,actual);
    
    [cvbstim, cvbintstim]= regress(NormBinStim(test(cstim,i)), RegTbl(test(cstim,i),3), 0.01);
    pred = cvbstim*RegTbl(test(cstim,i),3);
    cvEvalStim_onlyRewardSize(i) = cvEvalFunc(pred,actual);
    
    
end


for i = 1:5      % testing reduced models for action 
    actual = NormBinAction(test(cact,i));
    
    [cvbact, cvbintact]= regress(NormBinAction(test(cact,i)), RegTbl(test(cact,i),[2 3]), 0.01);
    pred = cvbact(1)*RegTbl(test(cact,i),2) + cvbact(2)*RegTbl(test(cact,i),3);
    cvEvalAct_noIpsiContra(i) = cvEvalFunc(pred,actual);
    
    [cvbact, cvbintact]= regress(NormBinAction(test(cact,i)), RegTbl(test(cact,i),[1 3]), 0.01);
    pred = cvbact(1)*RegTbl(test(cact,i),1) + cvbact(2)*RegTbl(test(cact,i),3);
    cvEvalAct_noContrast(i) = cvEvalFunc(pred,actual);
    
    [cvbact, cvbintact]= regress(NormBinAction(test(cact,i)), RegTbl(test(cact,i),[1 2]), 0.01);
    pred = cvbact(1)*RegTbl(test(cact,i),1) + cvbact(2)*RegTbl(test(cact,i),2);
    cvEvalAct_noRewardSize(i) = cvEvalFunc(pred,actual);
    
    [cvbact, cvbintact]= regress(NormBinAction(test(cact,i)), RegTbl(test(cact,i),1), 0.01);
    pred = cvbact*RegTbl(test(cact,i),1);
    cvEvalAct_onlyIpsiContra(i) = cvEvalFunc(pred,actual);
    
    [cvbact, cvbintact]= regress(NormBinAction(test(cact,i)), RegTbl(test(cact,i),2), 0.01);
    pred = cvbact*RegTbl(test(cact,i),2);
    cvEvalAct_onlyContrast(i) = cvEvalFunc(pred,actual);
    
    [cvbact, cvbintact]= regress(NormBinAction(test(cact,i)), RegTbl(test(cact,i),3), 0.01);
    pred = cvbact*RegTbl(test(cact,i),3);
    cvEvalAct_onlyRewardSize(i) = cvEvalFunc(pred,actual);
    
    
end

for i = 1:5      % testing reduced models for reward
    actual = NormBinReward(test(crwd,i));
    
    [cvbrwd, cvbintact]= regress(NormBinReward(test(crwd,i)), RegTbl(test(crwd,i),[2 3]), 0.01);
    pred = cvbrwd(1)*RegTbl(test(crwd,i),2) + cvbrwd(2)*RegTbl(test(crwd,i),3);
    cvEvalRwd_noIpsiContra(i) = cvEvalFunc(pred,actual);
    
    [cvbrwd, cvbintact]= regress(NormBinReward(test(crwd,i)), RegTbl(test(crwd,i),[1 3]), 0.01);
    pred = cvbrwd(1)*RegTbl(test(crwd,i),1) + cvbrwd(2)*RegTbl(test(crwd,i),3);
    cvEvalRwd_noContrast(i) = cvEvalFunc(pred,actual);
    
    [cvbrwd, cvbintact]= regress(NormBinReward(test(crwd,i)), RegTbl(test(crwd,i),[1 2]), 0.01);
    pred = cvbrwd(1)*RegTbl(test(crwd,i),1) + cvbrwd(2)*RegTbl(test(crwd,i),2);
    cvEvalRwd_noRewardSize(i) = cvEvalFunc(pred,actual);
    
    [cvbrwd, cvbintact]= regress(NormBinReward(test(crwd,i)), RegTbl(test(crwd,i),1), 0.01);
    pred = cvbrwd*RegTbl(test(crwd,i),1);
    cvEvalRwd_onlyIpsiContra(i) = cvEvalFunc(pred,actual);
    
    [cvbrwd, cvbintact]= regress(NormBinReward(test(crwd,i)), RegTbl(test(crwd,i),2), 0.01);
    pred = cvbrwd*RegTbl(test(crwd,i),2);
    cvEvalRwd_onlyContrast(i) = cvEvalFunc(pred,actual);
    
    [cvbrwd, cvbintact]= regress(NormBinReward(test(crwd,i)), RegTbl(test(crwd,i),3), 0.01);
    pred = cvbrwd*RegTbl(test(crwd,i),3);
    cvEvalRwd_onlyRewardSize(i) = cvEvalFunc(pred,actual);
    
    
end


figure; 
subplot(1,3,1)
title('Stimulus regrssions')
barh([mean(cvEvalStim_noContrast), mean(cvEvalStim_noIpsiContra), mean(cvEvalStim_noRewardSize); mean(cvEvalStim_onlyContrast), mean(cvEvalStim_onlyIpsiContra), mean(cvEvalStim_onlyRewardSize)]); 
hold on
% barh([mean(cvEval_onlyContrast), mean(cvEval_onlyIpsiContra), mean(cvEval_onlyRewardSize)]); 
if strcmpi(brain_region, 'DMS')
    xlim([0.1 round(mean(cvEvalStim), 1)])
elseif strcmpi(brain_region, 'NAC')
    xlim([0.03 round(mean(cvEvalStim), 1)])
else 
    xlim([0 round(mean(cvEvalStim), 1)])
end
yticklabels({'without', 'only'})
line([mean(cvEvalStim) mean(cvEvalStim)], ylim, 'LineWidth', 2, 'LineStyle', '--', 'color', [0.5 0.5 0.5])
legend({'Contrast', 'Ipsi/Contra', 'Reward size', 'Full model'})


subplot(1,3,2)
title('Action regressions')
barh([mean(cvEvalAct_noContrast), mean(cvEvalAct_noIpsiContra), mean(cvEvalAct_noRewardSize); mean(cvEvalAct_onlyContrast), mean(cvEvalAct_onlyIpsiContra), mean(cvEvalAct_onlyRewardSize)]); 
hold on
if strcmpi(brain_region, 'DMS')
    xlim([0.04 0.07])
elseif strcmpi(brain_region, 'NAC')
    xlim([0 round(mean(cvEvalAct), 1)])
else 
    xlim([0 round(mean(cvEvalAct), 1)])
end
yticklabels({'without', 'only'})
line([mean(cvEvalAct) mean(cvEvalAct)], ylim, 'LineWidth', 2, 'LineStyle', '--', 'color', [0.5 0.5 0.5])
legend({'Contrast', 'Ipsi/Contra', 'Reward size', 'Full model'})


subplot(1,3,3)
title('Reward regressions')
barh([mean(cvEvalRwd_noContrast), mean(cvEvalRwd_noIpsiContra), mean(cvEvalRwd_noRewardSize); mean(cvEvalRwd_onlyContrast), mean(cvEvalRwd_onlyIpsiContra), mean(cvEvalRwd_onlyRewardSize)]); 
hold on
if strcmpi(brain_region, 'DMS')
    xlim([0.04 round(mean(cvEvalRwd), 1)])
elseif strcmpi(brain_region, 'NAC')
    xlim([0 round(mean(cvEvalRwd), 1)])
else 
    xlim([0 round(mean(cvEvalRwd), 1)])
end
yticklabels({'without', 'only'})
line([mean(cvEvalRwd) mean(cvEvalRwd)], ylim, 'LineWidth', 2, 'LineStyle', '--', 'color', [0.5 0.5 0.5])
legend({'Contrast', 'Ipsi/Contra', 'Reward size', 'Full model'})






%%

function [NormBinStim, NormBinAction, NormBinReward] = getResponses(animal, NormBinStim, NormBinAction, NormBinReward, StimData, ActionData, RewardData)


NormBinReward = [NormBinReward; mean(RewardData(:,4000:5300),2) - mean(RewardData(:,3400:3700),2)];


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
