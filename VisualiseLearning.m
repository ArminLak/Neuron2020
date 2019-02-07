% Visualise avg responses to stimulus and reward in specified sessions

% Armin Oct 2018

% 14 jan 2019; to do: make the code work with L and R hem
% to do instead of z-scoring normalise to max response

clear all
close all



animal_list= [{'ALK068', 'ALK070', 'ALK071','ALK084'}];


%animal_list= [{'ALK078', 'MMM001', 'MMM002'}];


animal_list= [{'ALK084'}];

BrainStrucutre = 'VTA'

%animal_list= [{'ALK068'}];


% load the dataset
load(['BehPhotoM_Exp7', '_', BrainStrucutre]);


% plot colours
colorGray4 = [0.8 0.8 0.8 %lightest
    0.6 0.6 0.6
    0.4 0.4 0.4
    0.2 0.2 0.2
    0 0 0]; % black
colorGray3 = [0.8 0.8 0.8
    0.4 0.4 0.4
    0 0 0];
colorGray2 = [0.7 0.7 0.7
    0 0 0];
colorGreen = [0 1 0
    0 0.8 0
    0 0.6 0
    0  0.3 0];

% pre-allocation
NeuronRewardErr=nan(length(animal_list),20,5,13100);
NeuronRewardCor=nan(length(animal_list),20,5,13100);

NeuronStimCor =nan(length(animal_list),20,5,13100);
NeuronStimErr =nan(length(animal_list),20,5,13100);

NeuronActionCor =nan(length(animal_list),20,5,13100);
NeuronActionErr =nan(length(animal_list),20,5,13100);

Performance = nan(length(animal_list),20,5);
CorStimRew = nan(4,20);

for animalcount = 1:length(animal_list)
    
    animal_name = animal_list(animalcount);
    [animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);
    
    
    if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardR') && isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')
        hem = 2;
    else hem = 1;
    end
       
    SessionList = 1:length(BehPhotoM(animal_ID).Session);
    
    for ihem = 1:hem
            
        for iSession = SessionList
            
            if  hem == 1
                
                if isfield(BehPhotoM(animal_ID).Session,'NeuronRewardL')
                    TrialTimingData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
                    NeuronStim = BehPhotoM(animal_ID).Session(iSession).NeuronStimL;
                    NeuronAction = BehPhotoM(animal_ID).Session(iSession).NeuronActionL;
                    NeuronReward = BehPhotoM(animal_ID).Session(iSession).NeuronRewardL;
                else
                    TrialTimingData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
                    NeuronStim = BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
                    NeuronAction = BehPhotoM(animal_ID).Session(iSession).NeuronActionR;
                    NeuronReward = BehPhotoM(animal_ID).Session(iSession).NeuronRewardR;
                    
                end
                
            end
            
            if  hem == 2
                
                if ihem ==1
                    TrialTimingData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
                    NeuronStim = BehPhotoM(animal_ID).Session(iSession).NeuronStimL;
                    NeuronAction = BehPhotoM(animal_ID).Session(iSession).NeuronActionL;
                    NeuronReward = BehPhotoM(animal_ID).Session(iSession).NeuronRewardL;
                elseif ihem ==2
                    TrialTimingData = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
                    NeuronStim = BehPhotoM(animal_ID).Session(iSession).NeuronStimR;
                    NeuronAction = BehPhotoM(animal_ID).Session(iSession).NeuronActionR;
                    NeuronReward = BehPhotoM(animal_ID).Session(iSession).NeuronRewardR;
                    
                end
                
            end
            
            
            NeuronStimTurned=NeuronStim';
            NeuronStimVector=NeuronStimTurned(:);
            NeuronStimVector = zscore(NeuronStimVector);
            NeuronStimNorm = reshape(NeuronStimVector,size(NeuronStim,2),size(NeuronStim,1));
            NeuronStimNorm = NeuronStimNorm';
            
            NeuronStim = NeuronStimNorm;
            
            
            NeuronActionTurned=NeuronAction';
            NeuronActionVector=NeuronActionTurned(:);
            NeuronActionVector = zscore(NeuronActionVector);
            NeuronActionNorm = reshape(NeuronActionVector,size(NeuronAction,2),size(NeuronAction,1));
            NeuronActionNorm = NeuronActionNorm';
               
            NeuronAction = NeuronActionNorm;
            
            
            NeuronRewardTurned=NeuronReward';
            NeuronRewardVector=NeuronRewardTurned(:);
            NeuronRewardVector = zscore(NeuronRewardVector);
            NeuronRewardNorm = reshape(NeuronRewardVector,size(NeuronReward,2),size(NeuronReward,1));
            NeuronRewardNorm = NeuronRewardNorm';
            
            NeuronReward = NeuronRewardNorm;
            
            % we need better stim 
            c = 5;
            for istim = fliplr(unique(abs(TrialTimingData(:,2))'))
                
                NeuronStimCor(animalcount,iSession,c,:)   = mean(NeuronStim(TrialTimingData(:,9)==1 & abs(TrialTimingData(:,2))==istim, :));
                NeuronActionCor(animalcount,iSession,c,:) = mean(NeuronAction(TrialTimingData(:,9)==1 & abs(TrialTimingData(:,2))==istim, :));
                
                NeuronRewardCor(animalcount,iSession,c,:) = mean(NeuronReward(TrialTimingData(:,9)==1 & abs(TrialTimingData(:,2))==istim, :));
                
                NeuronStimErr(animalcount,iSession,c,:)   = mean(NeuronStim(TrialTimingData(:,9)==0 & abs(TrialTimingData(:,2))==istim, :));
                NeuronActionErr(animalcount,iSession,c,:) = mean(NeuronAction(TrialTimingData(:,9)==0 & abs(TrialTimingData(:,2))==istim, :));
                
                NeuronRewardErr(animalcount,iSession,c,:) = mean(NeuronReward(TrialTimingData(:,9)==0 & abs(TrialTimingData(:,2))==istim, :));
                
                Performance(animalcount,iSession,c) = mean(TrialTimingData(abs(TrialTimingData(:,2))==istim, 9));
                
                
                c=c-1;
                
            end
           
            RespStim=mean(NeuronStim(:,3900:4400)')-mean(NeuronStim(:,3600:3800)');
            RespOutcome=mean(NeuronReward(:,3900:4600)')-mean(NeuronReward(:,3400:3700)');
         
            figure(animalcount)
     hold on
     subplot(length(SessionList),1,iSession)
     
     scatter(RespStim,RespOutcome)
     
    CorStimRew(animalcount,iSession) = corr(RespStim',RespOutcome');
    
     
     xlim([-2, 4])
     ylim([-2,4])
        end
     
     
    end
    
    NSessionPerAnimal(animalcount) =  iSession;
end


%%
figure
  popStimResp = nan(9,5); % arbit size
    popRewResp = nan(9,5); % arbit size
    
if length(animal_list)==1
    
  
    
    for i=1:NSessionPerAnimal(1)
        
        subplot(NSessionPerAnimal(1),3,3*i-2); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronStimCor(:,i,istim,:)),2),'color',colorGray4(istim,:))    ; hold on
            
            TemppopStimResp=nanmean(squeeze(NeuronStimCor(:,i,istim,:)),2);
            popStimResp(i,istim) =nanmean(TemppopStimResp(3900:4400))-nanmean(TemppopStimResp(3600:3800));
            
            
        end
        xlim([3700 4500])
        ylim([-1 2])
        
        subplot(NSessionPerAnimal(1),3,3*i-1); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronActionCor(:,i,istim,:)),2),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        xlim([3500 4500])
        ylim([-1 2])
        
        subplot(NSessionPerAnimal(1),3,3*i); hold on
        
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronRewardCor(:,i,istim,:)),2),'color',colorGray4(istim,:))    ; hold on
            
            TemppopRewResp=nanmean(squeeze(NeuronRewardCor(:,i,istim,:)),2);
            popRewResp(i,istim) =nanmean(TemppopRewResp(3900:4600)) - nanmean(TemppopRewResp(3400:3700));
            
            
        end
        ylim([-1 2])
        
        xlim([3500 4500])
    end
    
    
    figure
    hold on
    
    plot(popStimResp(:,5))
    plot(popStimResp(:,4))
    plot(popStimResp(:,3))
    plot(popStimResp(:,2))
    
    
    figure
    hold on
    
    plot(popRewResp(:,5))
    plot(popRewResp(:,4))
    plot(popRewResp(:,3))
    plot(popRewResp(:,2))
    
else
    
    for i=1:13
        
        subplot(13,3,3*i-2); hold on
        for istim =1:5
            
           plot(nanmean(squeeze(NeuronStimCor(:,i,istim,:))),'color',colorGray4(istim,:))    ; hold on
            TemppopStimResp=nanmean(squeeze(NeuronStimCor(:,i,istim,:)));
            popStimResp(i,istim) =nanmean(TemppopStimResp(3900:4400))-nanmean(TemppopStimResp(3600:3800));
            
        end
        xlim([3700 4500])
        ylim([-1 2])
        
        subplot(13,3,3*i-1); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronActionCor(:,i,istim,:))),'color',colorGray4(istim,:))    ; hold on
            
            
        end
        xlim([3500 4500])
        ylim([-1 2])
        
        subplot(13,3,3*i); hold on
        for istim =1:5
            
            plot(nanmean(squeeze(NeuronRewardCor(:,i,istim,:))),'color',colorGray4(istim,:))    ; hold on
            
              TemppopRewResp=nanmean(squeeze(NeuronRewardCor(:,i,istim,:)));
            popRewResp(i,istim) =nanmean(TemppopRewResp(3900:4600)) - nanmean(TemppopRewResp(3400:3700));
         
            
        end
        ylim([-1 2])
        
        xlim([3500 4500])
    end
    
    
    figure
    hold on
    
     popStimResp(3,4)=nan;
%     popStimResp(10,:)=nan;
%     
     popRewResp (3,4) = nan;
%     popRewResp (10,:) = nan;
%     
     Performance(1,3,4) = nan;
%     
    
    plot(popStimResp(:,5))
    plot(popStimResp(:,4))
    plot(popStimResp(:,3))
    plot(popStimResp(:,2))
    
    figure
    hold on
    
    plot(popRewResp(:,5))
    plot(popRewResp(:,4))
    plot(popRewResp(:,3))
    plot(popRewResp(:,2))
    
end


figure
plot(nanmean(Performance(:,:,5)))
hold on
plot(nanmean(Performance(:,:,4)))
plot(nanmean(Performance(:,:,3)))
plot(nanmean(Performance(:,:,2)))


% figure
% for animalcount = 1:length(animal_list)
% 
% 
% a  = squeeze(NeuronStimCor(animalcount,:,5,:));
% b =     max(mean(a(:,4000:4500)'));
% 
%  subplot(9,3,3*i-2); hold on
%   for istim =1:5
%             
%            plot(nanmean(squeeze(NeuronStimCor(:,i,istim,:))),'color',colorGray4(istim,:))    ; hold on
%             TemppopStimResp=nanmean(squeeze(NeuronStimCor(:,i,istim,:)));
%             popStimResp(i,istim) =nanmean(TemppopStimResp(3900:4400))-nanmean(TemppopStimResp(3600:3800));
%             
%             popStimResp(i,istim) = popStimResp(i,istim) ./ b;
%         end
% 
% 
% end
%  
% figure
%     hold on
%     
%     plot(popStimResp(:,5))
%     plot(popStimResp(:,4))
%     plot(popStimResp(:,3))
%     plot(popStimResp(:,2))




