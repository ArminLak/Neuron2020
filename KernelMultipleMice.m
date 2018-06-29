thiclear all
close all

% list of animals
Animals = [48 50 51]

%Animals = [51]

load('BehPhotoM_Exp23')

 

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

%%
c=1;
for iAnimal = Animals

BeepKernel(c,:) = BehPhotoM(iAnimal).KernelSummary.Kernels{1} ;
StimKernel(c,:) = BehPhotoM(iAnimal).KernelSummary.Kernels{2} ;
ActKernel(c,:) = BehPhotoM(iAnimal).KernelSummary.Kernels{3} ;
RewKernel(c,:) = BehPhotoM(iAnimal).KernelSummary.Kernels{4} ;


resp_CorrStim(c,:)=BehPhotoM(iAnimal).KernelSummary.resp_CorrStim ./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim) ;
resp_ErrStim(c,:) = BehPhotoM(iAnimal).KernelSummary.resp_ErrStim ./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);


resp_CorrAct(c,:) = BehPhotoM(iAnimal).KernelSummary.resp_CorrAct./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);
resp_ErrAct(c,:) = BehPhotoM(iAnimal).KernelSummary.resp_ErrAct./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);



resp_CorrRew(c,:) = BehPhotoM(iAnimal).KernelSummary.resp_CorrRew ./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);
resp_ErrRew(c,:)= BehPhotoM(iAnimal).KernelSummary.resp_ErrRew./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim) ;

%resp_CorrRew(c,:) = BehPhotoM(iAnimal).KernelSummary.resp_CorrRew ./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrRew);
%resp_ErrRew(c,:)= BehPhotoM(iAnimal).KernelSummary.resp_ErrRew ./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrRew);


resp_Block1StimCor(c,:)=BehPhotoM(iAnimal).KernelSummary.resp_Block1StimCor./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);
resp_Block2StimCor(c,:)=BehPhotoM(iAnimal).KernelSummary.resp_Block2StimCor./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);

resp_Block1ActCor(c,:)=BehPhotoM(iAnimal).KernelSummary.resp_Block1ActCor./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);
resp_Block2ActCor(c,:)=BehPhotoM(iAnimal).KernelSummary.resp_Block2ActCor./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);

resp_Block1RewCor(c,:)=BehPhotoM(iAnimal).KernelSummary.resp_Block1RewCor./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);
resp_Block2RewCor(c,:)=BehPhotoM(iAnimal).KernelSummary.resp_Block2RewCor./ max(BehPhotoM(iAnimal).KernelSummary.resp_CorrStim);

PredictedStim1(c,:) = BehPhotoM(iAnimal).KernelSummary.PredictedPopStimAlign(1,:);
PredictedStim2(c,:) = BehPhotoM(iAnimal).KernelSummary.PredictedPopStimAlign(2,:);
PredictedStim3(c,:) = BehPhotoM(iAnimal).KernelSummary.PredictedPopStimAlign(3,:);
PredictedStim4(c,:) = BehPhotoM(iAnimal).KernelSummary.PredictedPopStimAlign(4,:);


PredictedRewL(c,:) = BehPhotoM(iAnimal).KernelSummary.PredictedPopRewAlignLarge;
PredictedRewS(c,:) = BehPhotoM(iAnimal).KernelSummary.PredictedPopRewAlignSmall;
PredictedRewNo(c,:) = BehPhotoM(iAnimal).KernelSummary.PredictedPopRewAlignNoRew;



c=c+1;
end

%%
figure % generate pop predicted responses
subplot(1,2,1)
plot(smooth(nanmean(PredictedStim1)))
hold on
plot(smooth(nanmean(PredictedStim2)))
plot(smooth(nanmean(PredictedStim3)))
plot(smooth(nanmean(PredictedStim4)))
xlim([4 33])

subplot(1,2,2)
plot(smooth(nanmean(PredictedRewL)))
hold on
plot(smooth(nanmean(PredictedRewS)))
plot(smooth(nanmean(PredictedRewNo)))
xlim([6 38])



resp_Block2ActCor(3,end) = nan;
resp_ErrAct(3,end) = nan;

resp_CorrRew(1,4) = nan;

resp_Block1RewCor(1,4) = nan;

resp_Block1RewCor(2,4) = nan;

resp_Block2RewCor(1,4) = nan;



figure

subplot(1,4,1)
plot(nanmean(BeepKernel(:,1:end-1),1))

subplot(1,4,2)
plot(nanmean(StimKernel(:,2:end),1))

subplot(1,4,3)
plot(nanmean(ActKernel,1))

subplot(1,4,4)
plot(nanmean(RewKernel,1))


% figure
% 
% windows = {[-400 800];[-200 1600];[-600 200];[200 2000]}; %
% 
% subplot(1,4,1)
% plot(-300 :50: 800, nanmean(BeepKernel(1,1:end-1),1))
% 
% subplot(1,4,2)
% plot(-100 :50: 1600,nanmean(StimKernel(1,2:end),1))
% 
% subplot(1,4,3)
% plot(-500 :50: 250, nanmean(ActKernel(1,:),1))
% 
% subplot(1,4,4)
% plot(0 :50: 1750, nanmean(RewKernel(1,:),1))



figure

subplot(2,3,1); hold on

ylim([-0.5 1.5])


 plot(StimAllowed,nanmean(resp_CorrStim,1),'color',[0 1 0],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_CorrStim,1),nanstd(resp_CorrStim,1)/2,'color',[0 1 0],...
        'marker','o','markersize',7,'markeredgecolor',[0 1 0],...
        'markerfacecolor',[0 1 0],'linestyle','none','linewidth',1.2);
   
    plot(StimAllowed,nanmean(resp_ErrStim,1),'color',[1 0 0],...
        'linestyle','-','linewidth',1.2);
    
    errorbar(StimAllowed,nanmean(resp_ErrStim,1),nanstd(resp_ErrStim,1)/2,'color',[1 0 0],...
        'marker','o','markersize',7,'markeredgecolor',[1 0 0],...
        'markerfacecolor',[1 0 0],'linestyle','none','linewidth',1.2);
    
    
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

