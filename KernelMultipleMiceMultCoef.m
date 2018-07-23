clear all
close all

% list of animals
Animals = [48 50 51]


load('BehPhotoM_Exp23')
%%

 ModelArrangement=11

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
    0 0 0
    ];

StimAllowed = [-0.5 -0.25 -0.12 0 0.12 0.25 0.5];
absStim = unique(abs(StimAllowed));

%
c=1;
for iAnimal = Animals

StimKernel(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).Kernels{1} ;
RewKernel(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).Kernels{2} ;


% resp_CorrStimLarge(c,:)=BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrStimLarge ./ max(BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrStimLarge) ;
% resp_ErrStimLarge(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_ErrStimLarge ./ max(BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrStimLarge);
% resp_CorrStimSmall(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrStimSmall./ max(BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrStimLarge);
% 
% 
% resp_CorrOutcomeLarge(c,:)=BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrOutcomeLarge ./ max(BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrOutcomeLarge) ;
% resp_ErrOutcomeLarge(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_ErrOutcomeLarge ./ max(BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrOutcomeLarge);
% resp_CorrOutcomeSmall(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrOutcomeSmall./ max(BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrOutcomeLarge);

resp_CorrStimLarge(c,:)=BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrStimLarge  ;
resp_ErrStimLarge(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_ErrStimLarge ;
resp_CorrStimSmall(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrStimSmall;


resp_CorrOutcomeLarge(c,:)=BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrOutcomeLarge  ;
resp_ErrOutcomeLarge(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_ErrOutcomeLarge ;
resp_CorrOutcomeSmall(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).resp_CorrOutcomeSmall;


PredictedStim1(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopStimAlign(1,:);
PredictedStim2(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopStimAlign(2,:);
PredictedStim3(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopStimAlign(3,:);
PredictedStim4(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopStimAlign(4,:);

PredictedAction1(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopActionAlign(1,:);
PredictedAction2(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopActionAlign(2,:);
PredictedAction3(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopActionAlign(3,:);
PredictedAction4(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopActionAlign(4,:);


PredictedRewL(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopRewAlignLarge;
PredictedRewS(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopRewAlignSmall;
PredictedRewNo(c,:) = BehPhotoM(iAnimal).KernelSummary(ModelArrangement).PredictedPopRewAlignNoRew;



c=c+1;
end

%%
figure % generate pop predicted responses
subplot(1,3,1)
plot(smooth(nanmean(PredictedStim1)))
hold on
plot(smooth(nanmean(PredictedStim2)))
plot(smooth(nanmean(PredictedStim3)))
plot(smooth(nanmean(PredictedStim4)))
xlim([7 26])

subplot(1,3,2)
plot(smooth(nanmean(PredictedAction1)))
hold on
plot(smooth(nanmean(PredictedAction2)))
plot(smooth(nanmean(PredictedAction3)))
plot(smooth(nanmean(PredictedAction4)))
xlim([4 23])

subplot(1,3,3)
plot(smooth(nanmean(PredictedRewL)))
hold on
plot(smooth(nanmean(PredictedRewS)))
plot(smooth(nanmean(PredictedRewNo)))
xlim([7 26])



% figure; hold on
% 
% 
%     plot(absStim,nanmean(resp_CorrStimLarge,1),'color',[0 1 0],...
%         'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
%         'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
% plot(absStim,nanmean(resp_ErrStimLarge,1),'color',[1 0 0],...
%         'marker','o','markersize',7,'markeredgecolor',[1 0 0],...
%         'markerfacecolor',[1 0 0],'linestyle','none','linewidth',1.2);
%    
%      plot(absStim,nanmean(resp_CorrStimSmall,1),'color',[0 1 0],...
%         'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
%         'markerfacecolor',[0 1 0],'linestyle','--','linewidth',1.2);
% 
%     
%     figure; hold on
%     
%       plot(absStim,nanmean(resp_CorrOutcomeLarge,1),'color',[0 1 0],...
%         'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
%         'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
%    
%   plot(absStim,nanmean(resp_ErrOutcomeLarge,1),'color',[1 0 0],...
%         'marker','o','markersize',7,'markeredgecolor',[1 0 0],...
%         'markerfacecolor',[1 0 0],'linestyle','none','linewidth',1.2);
%    
%      plot(absStim,nanmean(resp_CorrOutcomeSmall,1),'color',[0 1 0],...
%         'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
%         'markerfacecolor',[0 1 0],'linestyle','--','linewidth',1.2);
    
figure

subplot(2,3,1); hold on

ylim([-0.5 1.5])





    errorbar(absStim,nanmean(resp_CorrStimLarge,1),nanstd(resp_CorrStimLarge,1)/2,'color',[0 1 0],...
        'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
        'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
   
  errorbar(absStim,nanmean(resp_ErrStimLarge,1),nanstd(resp_ErrStimLarge,1)/2,'color',[1 0 0],...
        'marker','o','markersize',7,'markeredgecolor',[1 0 0],...
        'markerfacecolor',[1 0 0],'linestyle','none','linewidth',1.2);
   
     errorbar(absStim,nanmean(resp_CorrStimSmall,1),nanstd(resp_CorrStimSmall,1)/2,'color',[0 1 0],...
        'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
        'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
 
    
    subplot(2,3,2); hold on

      errorbar(absStim,nanmean(resp_CorrOutcomeLarge,1),nanstd(resp_CorrOutcomeLarge,1)/2,'color',[0 1 0],...
        'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
        'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
   
  errorbar(absStim,nanmean(resp_ErrOutcomeLarge,1),nanstd(resp_ErrOutcomeLarge,1)/2,'color',[1 0 0],...
        'marker','o','markersize',7,'markeredgecolor',[1 0 0],...
        'markerfacecolor',[1 0 0],'linestyle','none','linewidth',1.2);
   
     errorbar(absStim,nanmean(resp_CorrOutcomeSmall,1),nanstd(resp_CorrOutcomeSmall,1)/2,'color',[0 1 0],...
        'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
        'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
    %%
    
    
subplot(2,3,2)

plot(StimAllowed,nanmean(resp_CorrAct,1),'g')
hold on
plot(StimAllowed,nanmean(resp_ErrAct,1),'r')
ylim([-0.5 1.5])


subplot(2,3,3); hold on

ylim([-0.5 1.5])

 plot(StimAllowed,nanmean(resp_CorrRew,1),'color',[0 1 0],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_CorrRew,1),nanstd(resp_CorrRew,1)/2,'color',[0 1 0],...
        'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
        'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
   
    plot(StimAllowed,nanmean(resp_ErrRew,1),'color',[1 0 0],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_ErrRew,1),nanstd(resp_ErrRew,1)/2,'color',[1 0 0],...
        'marker','o','markersize',7,'markeredgecolor',[1 0 0],...
        'markerfacecolor',[1 0 0],'linestyle','none','linewidth',1.2);
    

subplot(2,3,4)
%errorbar(StimAllowed,nanmean(resp_Block1StimCor,1),nanstd(resp_Block1StimCor,1)/2,'color',[0.5 0.3 0.12])
hold on
%errorbar(StimAllowed,nanmean(resp_Block2StimCor,1),nanstd(resp_Block2StimCor,1)/2,'color',[1 0.56 0.11])
ylim([-0.5 1.5])

 plot(StimAllowed,nanmean(resp_Block1StimCor,1),'color',[0.5 0.3 0.12],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_Block1StimCor,1),nanstd(resp_Block1StimCor,1)/2,'color',[0.5 0.3 0.12],...
        'marker','o','markersize',7,'markeredgecolor',[0.5 0.3 0.12],...
        'markerfacecolor',[0.5 0.3 0.12],'linestyle','none','linewidth',1.2);
   
    plot(StimAllowed,nanmean(resp_Block2StimCor,1),'color',[1 0.56 0.11],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_Block2StimCor,1),nanstd(resp_Block2StimCor,1)/2,'color',[1 0.56 0.11],...
        'marker','o','markersize',7,'markeredgecolor',[1 0.56 0.11],...
        'markerfacecolor',[1 0.56 0.11],'linestyle','none','linewidth',1.2);
 
    
    
subplot(2,3,5); hold on


plot(StimAllowed,nanmean(resp_Block1ActCor,1),'color',[0.5 0.3 0.12])
hold on
plot(StimAllowed,nanmean(resp_Block2ActCor,1),'color',[1 0.56 0.11])
ylim([-0.5 1.5])



subplot(2,3,6); hold on

plot(StimAllowed,nanmean(resp_Block1RewCor,1),'color',[0.5 0.3 0.12],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_Block1RewCor,1),nanstd(resp_Block1RewCor,1)/2,'color',[0.5 0.3 0.12],...
        'marker','o','markersize',7,'markeredgecolor',[0.5 0.3 0.12],...
        'markerfacecolor',[0.5 0.3 0.12],'linestyle','none','linewidth',1.2);
   
    plot(StimAllowed,nanmean(resp_Block2RewCor,1),'color',[1 0.56 0.11],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_Block2RewCor,1),nanstd(resp_Block2RewCor,1)/2,'color',[1 0.56 0.11],...
        'marker','o','markersize',7,'markeredgecolor',[1 0.56 0.11],...
        'markerfacecolor',[1 0.56 0.11],'linestyle','none','linewidth',1.2);
    

ylim([-0.5 2.5])

