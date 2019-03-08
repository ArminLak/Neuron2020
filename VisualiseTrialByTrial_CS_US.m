clear all
close all

% plot trial-by-trial scatter plots of Stim resposne vs outcome response 

%VTA : [48, 50,51, 64, 69]  coresponding to ALK068, 70 and 71, ALK084,
%ALK085, 

% select animal

 Animals = [48 50 51 64 69]
 
BrainStrucutre = 'VTA'
ExpID = '23'
load(['BehPhotoM_Exp', ExpID, '_', BrainStrucutre]);

animalN = 1
for animal_ID = Animals
    
    


RTLimit = 6; % in s, excluding trials with RT longer than this

color = [
    1 0 0         % red
    1 0.5 0       % orange
    1 1 0         % yellow
    0.5 1 0.5     % light green
    0 1 1         % light blue
    0  0.5 1      % medium blue
    0 0 1];       % blue

colorGray = [ 0.8 0.8 0.8
    0.6 0.6 0.6
    0.4 0.4 0.4
    0 0 0];

colorRed = [ 1 0.8 0.8
    1 0.6 0.6
    1 0.4 0.4
    1 0 0];


%%
BehData = [];
BeepData = [];
StimData = [];
ActionData = [];
RewardData = [];

BeepDataL   = [];
StimDataL   = [];
ActionDataL = [];
RewardDataL = [];

BeepDataR   = [];
StimDataR   = [];
ActionDataR = [];
RewardDataR = [];


sessionz = 1:length(BehPhotoM(animal_ID).Session);

iter = 0;
if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')
    
    for iSession = sessionz % left hem
        
        TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        BehData = [BehData; TempBehData];
        
        % left
        TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepL;
        
        TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimL;
        
        TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronActionL;
        
        TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardL;
        
        BeepDataL = [BeepDataL;TempBeepData];
        
        StimDataL = [StimDataL;TempStimData];
        
        ActionDataL = [ActionDataL;TempActionData];
        
        RewardDataL = [RewardDataL;TempRewardData];
        
    end
    
    iter = iter + 1;
end


if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardR')
    BehData = [];
    
    for iSession = sessionz % left hem
        
        TempBehData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        BehData = [BehData; TempBehData];
        
        % right
        TempBeepData= BehPhotoM(animal_ID).Session(iSession).NeuronBeepR;
        
        TempStimData= BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
        
        TempActionData= BehPhotoM(animal_ID).Session(iSession).NeuronActionR;
        
        TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronRewardR;
        
        BeepDataR = [BeepDataR;TempBeepData];
        
        StimDataR = [StimDataR;TempStimData];
        
        ActionDataR = [ActionDataR;TempActionData];
        
        RewardDataR = [RewardDataR;TempRewardData];
        
        
    end
    
    iter = iter + 1;
    
end


RT = BehData(:,10) - BehData(:,13);

toRemove = find ( RT > RTLimit);

toRemove2= find(BehData(:,1) < 50);
toRemove = unique([toRemove; toRemove2]);

BehData(toRemove,:) = [];

for HemIter = 1:iter
    
    if iter == 2 && HemIter ==1
        
        BeepData = BeepDataL;
        StimData = StimDataL;
        ActionData = ActionDataL;
        RewardData = RewardDataL;
        
    end
    
    if iter == 2 && HemIter ==2
        
        BeepData = BeepDataR;
        StimData = StimDataR;
        ActionData = ActionDataR;
        RewardData = RewardDataR;
        
    end
    
    if iter == 1
        
        if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')
            
            BeepData = BeepDataL;
            StimData = StimDataL;
            ActionData = ActionDataL;
            RewardData = RewardDataL;
            
            
        elseif isfield(BehPhotoM(animal_ID).Session,'NeuronRewardR')
            
            BeepData = BeepDataR;
            StimData = StimDataR;
            ActionData = ActionDataR;
            RewardData = RewardDataR;
            
        end
        
    end
    
    
    BeepData(toRemove,:) = [];
    StimData(toRemove,:) = [];
    ActionData(toRemove,:) = [];
    RewardData(toRemove,:) = [];
    
    
    ToLargeR = find((BehData(:,3)==-1 & BehData(:,8)==1)  | ...
        (BehData(:,3)==1 & BehData(:,8)==2));
    
    BehData(ToLargeR,16)=1;
    
    ToSmallR = setdiff(1:size(BehData,1),ToLargeR)';
    
    BehData(ToSmallR, 16)=-1;
    
    abzStim = unique(abs(BehData(:,2)))';
    
  
    
  
    if animal_ID == 48
        NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2)- mean(StimData(:,3400:3800),2);
        
       % NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2);
        
        
    elseif animal_ID == 50
        NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,3100:3400),2);
        
        
    elseif animal_ID == 51
        
        %NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,4000:4200),2)- mean(StimData(:,3400:3800),2);
        NormBinStim = mean(StimData(:,4500:5000),2) - mean(StimData(:,3400:3800),2);
        
        
    elseif animal_ID == 56
        
        NormBinStim = mean(StimData(:,4500:5000),2)- mean(StimData(:,3100:3500),2);
        NormBinStim = mean(StimData(:,4600:5300),2);
        
        
    elseif animal_ID == 57
        
        NormBinStim = mean(StimData(:,4300:5000),2)- mean(StimData(:,3400:3800),2);

        NormBinStim = mean(StimData(:,4600:5300),2)- mean(StimData(:,3400:3800),2);

    elseif animal_ID == 59
        
        NormBinStim = mean(StimData(:,4400:4800),2);
        NormBinStim = mean(StimData(:,4900:5400),2);
        
    
        elseif animal_ID == 64 && strcmp(BrainStrucutre,'DMS')
        
        NormBinStim = mean(StimData(:,4500:5000),2);
        
    elseif animal_ID == 66
        
        NormBinStim = mean(StimData(:,5000:6000),2);
        
                
    elseif animal_ID == 69
        
        NormBinStim = mean(StimData(:,3900:5200),2)- mean(StimData(:,3400:3800),2);
             
        
    elseif animal_ID == 70
        
        NormBinStim = mean(StimData(:,4000:5000),2);
        
    else
        
        NormBinStim = mean(StimData(:,4500:5000),2)- mean(StimData(:,3400:3800),2);
        
    end
    
  
    
    
    
    NormBinReward = mean(RewardData(:,4000:5300),2) - mean(RewardData(:,3400:3800),2);
    
    if animal_ID == 48
        NormBinReward = mean(RewardData(:,4100:5000),2) - mean(RewardData(:,3800:3900),2);
    end
    
    if animal_ID == 50
        NormBinReward = mean(RewardData(:,4100:5000),2) - mean(RewardData(:,3800:3900),2);
    end
    
    if animal_ID == 56
        NormBinReward = mean(RewardData(:,4000:4800),2) - mean(RewardData(:,3400:3800),2);
    end
    
    if animal_ID == 57
        NormBinReward = mean(RewardData(:,4000:4800),2) - mean(RewardData(:,3200:3700),2);
    end
    
    if animal_ID == 59
        NormBinReward = mean(RewardData(:,4000:4800),2);
    end
    
        if animal_ID == 66 &&  strcmp(ExpID , '38')
 
        NormBinReward = mean(RewardData(:,4100:4600),2) - mean(RewardData(:,3400:3800),2);
        end

    
    
end

NormBinStim = zscore(NormBinStim);
NormBinReward = zscore(NormBinReward);



%%
figure(1)
for istim = abzStim

ax=gca
%scatter(NormBinStim(intersect(ToLargeR,find(BehData(:,9)==1))),NormBinReward(intersect(ToLargeR,find(BehData(:,9)==1))))
fit_ellipse( NormBinStim(mintersect(ToLargeR,find(BehData(:,9)==1),find(BehData(:,2)==istim))),NormBinReward(mintersect(ToLargeR,find(BehData(:,9)==1),find(BehData(:,2)==istim))),ax,[0 0 1])

hold on

%scatter(NormBinStim(intersect(ToSmallR,find(BehData(:,9)==1))),NormBinReward(intersect(ToSmallR,find(BehData(:,9)==1))))
fit_ellipse( NormBinStim(mintersect(ToSmallR,find(BehData(:,9)==1),find(BehData(:,2)==istim))),NormBinReward(mintersect(ToSmallR,find(BehData(:,9)==1),find(BehData(:,2)==istim))),ax,[0 1 0])

%scatter(NormBinStim(find(BehData(:,9)==0)),NormBinReward(find(BehData(:,9)==0)))
fit_ellipse( NormBinStim(mintersect(find(BehData(:,9)==0),find(BehData(:,2)==istim))),NormBinReward(mintersect(find(BehData(:,9)==0),find(BehData(:,2)==istim))),ax,[1 0 0])

end



%%

StimLevel=linspace(min(NormBinStim),max(NormBinStim),10);
StimInt =StimLevel(2) - StimLevel(1);

c = 1;
for StimLevelz = StimLevel(2:9)

   Resp(animalN,c)= nanmean(NormBinReward(intersect(find(NormBinStim > StimLevelz) ,find( NormBinStim < StimLevelz+StimInt)) ));

   c= c+1;
end

animalN = animalN +1;

end


figure 

plot(nanmean(Resp))