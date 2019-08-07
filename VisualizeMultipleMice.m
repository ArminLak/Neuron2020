clear all
close all

Hem2show = 'both' % 'L' 'R' or 'both'

% list of animals

% VTA,
Animals = [48 50 51 64]
load('BehPhotoM_Exp23_VTA')

% NAC
% Animals = [56 57 59 66]
% 
% load('BehPhotoM_Exp23_NAc')

% DMS
%Animals = [53, 62, 63, 71,72]  % 55 has 6 stimuli. so I will need to make some changes to be able to add this
%          53, 55,62, 63,64, 68, 70, 71, 72 
% 68 and 70 signals are not good, 64 the signal is ok but looks very
% strange

%load('BehPhotoM_Exp23_DMS')


TimingVisualise = [-0.2 0.8
                   -0.8, 0.2
                   -0.2, 0.8]; % stim, action, reward in s


sampleRate = 1200;
StartTime = 3700; % saved in the database.


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

GrandPopAbsStimResp = zeros(4,13100);
GrandPopAbsActionResp = zeros(4,13100);

GrandPopAbsStimRespMiddleStimCorrect= zeros(4,13100);
GrandPopAbsStimRespMiddleStimError= zeros(4,13100);


GrandPopAbsActionRespMiddleStimCorrect= zeros(4,13100);
GrandPopAbsActionRespMiddleStimError= zeros(4,13100);

GrandPopAbsStimRespMiddleStimCorLarge= zeros(1,13100);
GrandPopAbsStimRespMiddleStimCorSmall= zeros(1,13100);

GrandPopAbsActionRespMiddleStimCorLarge= zeros(1,13100);
GrandPopAbsActionRespMiddleStimCorSmall= zeros(1,13100);


GrandPopNormBinStimNoFold = zeros(2,7);
GrandPopNormBinRewardNoFold = zeros(2,7);
GrandPopNormBinStimNoFoldCorrError = zeros(2,7);
GrandPopNormBinRewardNoFoldCorrError = zeros(2,7);


GrandPopStim2AwayDACorError = zeros(4,13100);

GrandPopBeep2AwayDACorError = zeros(4,13100);

GrandPopRew2AwayDACorError = zeros(4,13100);

GrandPopStimBin = nan(4,4,length(Animals));

GrandPopRewBin = nan(4,4,length(Animals));

GrandPopStimLargeCorrect = zeros(7,13100);
GrandPopStimSmallCorrect = zeros(7,13100);
GrandPopStimLargeError   = zeros(7,13100);

GrandPopActionLargeCorrect = zeros(7,13100);
GrandPopActionSmallCorrect = zeros(7,13100);
GrandPopActionLargeError   = zeros(7,13100);


c=1;
animalCount = 1;
for iAnimal = Animals
    
    
    ChanN = 0;
    
    if ~isempty(BehPhotoM(iAnimal).GrandSummaryL)
        ChanN = ChanN + 1;
    end
    
    if ~isempty(BehPhotoM(iAnimal).GrandSummaryR)
        ChanN = ChanN + 1;
    end
    
    if ~strcmpi(Hem2show,'both')
        
        ChanN = 1;
    end
    
    for iChan = 1:ChanN
        
        BehPhotoM(iAnimal).GrandSummary=[]; 
        
        if ~isempty(BehPhotoM(iAnimal).GrandSummaryL) && strcmpi(Hem2show,'L')
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
        elseif  ~isempty(BehPhotoM(iAnimal).GrandSummaryR) && strcmpi(Hem2show,'R')
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
       
            
        end
        
        if strcmpi(Hem2show,'both')
            
            if ~isempty(BehPhotoM(iAnimal).GrandSummaryL)
                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
            end
            
            if ~isempty(BehPhotoM(iAnimal).GrandSummaryR)
                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
            end
            
            
            if ~isempty(BehPhotoM(iAnimal).GrandSummaryL) && ~isempty(BehPhotoM(iAnimal).GrandSummaryR)

            if iChan == 1
                                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;

            elseif iChan==2
                                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
            end
                 
            end
             
        end
        
        if            ~isempty(BehPhotoM(iAnimal).GrandSummary)
            
            PerBlock1(c,:)=BehPhotoM(iAnimal).GrandSummary.Performance(1,:);
            RTBlock1(c,:)=BehPhotoM(iAnimal).GrandSummary.RT(1,:);
            
            PerBlock2(c,:)=BehPhotoM(iAnimal).GrandSummary.Performance(2,:);
            RTBlock2(c,:)=BehPhotoM(iAnimal).GrandSummary.RT(2,:);
         
               
            % Stim-align rastersCor Large full stim
            SingleAnimalStimTraceLargeCorrect     = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeCorrect;
            SingleAnimalNormStimTraceLargeCorrect = SingleAnimalStimTraceLargeCorrect ./ max(max(SingleAnimalStimTraceLargeCorrect));                        
            GrandPopStimLargeCorrect              = SingleAnimalNormStimTraceLargeCorrect + GrandPopStimLargeCorrect ;
            
            % Stim-align rastersErr Large full stim
            SingleAnimalStimTraceLargeError      = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeError;
            SingleAnimalNormStimTraceLargeError  = SingleAnimalStimTraceLargeError ./ max(max(SingleAnimalStimTraceLargeCorrect));                        
            GrandPopStimLargeError               = SingleAnimalNormStimTraceLargeError + GrandPopStimLargeError ;
            
             % Stim-align rastersCor Small full stim
            SingleAnimalStimTraceSmallCorrect     = BehPhotoM(iAnimal).GrandSummary.StimRasterSmallCorrect;
            SingleAnimalNormStimTraceSmallCorrect = SingleAnimalStimTraceSmallCorrect ./ max(max(SingleAnimalStimTraceLargeCorrect));            
            GrandPopStimSmallCorrect              = SingleAnimalNormStimTraceSmallCorrect + GrandPopStimSmallCorrect ;
           
            
            % Stim-align rastersCor Small full action
            SingleAnimalActionTraceSmallCorrect     = BehPhotoM(iAnimal).GrandSummary.ActionRasterSmallCorrect;
            SingleAnimalNormActionTraceSmallCorrect = SingleAnimalActionTraceSmallCorrect ./ max(max(SingleAnimalActionTraceSmallCorrect));            
            GrandPopActionSmallCorrect              = SingleAnimalNormActionTraceSmallCorrect + GrandPopActionSmallCorrect ;
            
                        % Stim-align rastersCor Large full action
            SingleAnimalActionTraceLargeCorrect     = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeCorrect;
            SingleAnimalNormActionTraceLargeCorrect = SingleAnimalActionTraceLargeCorrect ./ max(max(SingleAnimalActionTraceLargeCorrect));                        
            GrandPopActionLargeCorrect              = SingleAnimalNormActionTraceLargeCorrect + GrandPopActionLargeCorrect ;
            
            % Stim-align rastersErr Large full action
            SingleAnimalActionTraceLargeError      = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeError;
            SingleAnimalNormActionTraceLargeError  = SingleAnimalActionTraceLargeError ./ max(max(SingleAnimalActionTraceLargeCorrect));                        
            GrandPopActionLargeError               = SingleAnimalNormActionTraceLargeError + GrandPopActionLargeError ;
            
            
            
            
            % Stim-align rastersCor
            SingleAnimalStimTrace= BehPhotoM(iAnimal).GrandSummary.AbsStimRaster;
            SingleAnimalNormStimTrace = SingleAnimalStimTrace ./ max(max(SingleAnimalStimTrace));
           
            GrandPopAbsStimResp = SingleAnimalNormStimTrace + GrandPopAbsStimResp ;
            
            % Action-align rastersCor
            SingleAnimalActionTrace= BehPhotoM(iAnimal).GrandSummary.AbsActionRaster;
            SingleAnimalNormActionTrace = SingleAnimalActionTrace ./ max(max(SingleAnimalActionTrace));
            
            if iAnimal==63
            SingleAnimalNormActionTrace = SingleAnimalActionTrace; end
            GrandPopAbsActionResp = SingleAnimalNormActionTrace + GrandPopAbsActionResp ;
            
            % Stim-align rasters(correct/error)
            SingleAnimalStimTraceMiddleStimCorrect(1,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterCorrect(1,:);
            SingleAnimalStimTraceMiddleStimCorrect(2,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterCorrect(2,:);
            SingleAnimalStimTraceMiddleStimCorrect(3,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterCorrect(3,:);
            SingleAnimalStimTraceMiddleStimCorrect(4,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterCorrect(4,:);
            
            SingleAnimalNormStimTraceMiddleStimCorrect = SingleAnimalStimTraceMiddleStimCorrect ./ max(max(SingleAnimalStimTraceMiddleStimCorrect));
            GrandPopAbsStimRespMiddleStimCorrect = SingleAnimalNormStimTraceMiddleStimCorrect + GrandPopAbsStimRespMiddleStimCorrect ;
            
            SingleAnimalStimTraceMiddleStimError(1,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterError(1,:);
            SingleAnimalStimTraceMiddleStimError(2,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterError(2,:);
            SingleAnimalStimTraceMiddleStimError(3,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterError(3,:);
            SingleAnimalStimTraceMiddleStimError(4,:)= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterError(4,:);
            
            
            SingleAnimalNormStimTraceMiddleStimError = SingleAnimalStimTraceMiddleStimError ./ max(max(SingleAnimalStimTraceMiddleStimCorrect));
            GrandPopAbsStimRespMiddleStimError = SingleAnimalNormStimTraceMiddleStimError + GrandPopAbsStimRespMiddleStimError ;
            
            
            % Action-align rasters(correct/error)
            SingleAnimalActionTraceMiddleStimCorrect(1,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterCorrect(1,:);
            SingleAnimalActionTraceMiddleStimCorrect(2,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterCorrect(2,:);
            SingleAnimalActionTraceMiddleStimCorrect(3,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterCorrect(3,:);
            SingleAnimalActionTraceMiddleStimCorrect(4,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterCorrect(4,:);
            
            SingleAnimalNormActionTraceMiddleStimCorrect = SingleAnimalActionTraceMiddleStimCorrect ./ max(max(SingleAnimalActionTraceMiddleStimCorrect));
            GrandPopAbsActionRespMiddleStimCorrect = SingleAnimalNormActionTraceMiddleStimCorrect + GrandPopAbsActionRespMiddleStimCorrect ;
            
            SingleAnimalActionTraceMiddleStimError(1,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterError(1,:);
            SingleAnimalActionTraceMiddleStimError(2,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterError(2,:);
            SingleAnimalActionTraceMiddleStimError(3,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterError(3,:);
            SingleAnimalActionTraceMiddleStimError(4,:)= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterError(4,:);
            
            
            SingleAnimalNormActionTraceMiddleStimError = SingleAnimalActionTraceMiddleStimError ./ max(max(SingleAnimalActionTraceMiddleStimCorrect));
            GrandPopAbsActionRespMiddleStimError = SingleAnimalNormActionTraceMiddleStimError + GrandPopAbsActionRespMiddleStimError ;
            
            
            
            % Stim-align rasters(middle stimuli, correct/ small/large reward)
            SingleAnimalStimTraceMiddleStimLargeCorrect= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterLargeCorrect(2,:);
            SingleAnimalNormStimTraceMiddleStimLargeCorrect = SingleAnimalStimTraceMiddleStimLargeCorrect ./ max(max(SingleAnimalStimTraceMiddleStimLargeCorrect));
            GrandPopAbsStimRespMiddleStimCorLarge = SingleAnimalNormStimTraceMiddleStimLargeCorrect + GrandPopAbsStimRespMiddleStimCorLarge ;
            
            SingleAnimalStimTraceMiddleStimSmallCorrect= BehPhotoM(iAnimal).GrandSummary.AbsStimRasterSmallCorrect(2,:);
            SingleAnimalNormStimTraceMiddleStimSmallCorrect = SingleAnimalStimTraceMiddleStimSmallCorrect ./ max(max(SingleAnimalStimTraceMiddleStimLargeCorrect));
            GrandPopAbsStimRespMiddleStimCorSmall = SingleAnimalNormStimTraceMiddleStimSmallCorrect + GrandPopAbsStimRespMiddleStimCorSmall ;
            
            
            
            % Action-align rasters(middle stimuli, correct/ small/large reward)
            SingleAnimalActionTraceMiddleStimLargeCorrect= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterLargeCorrect(2,:);
            SingleAnimalNormActionTraceMiddleStimLargeCorrect = SingleAnimalActionTraceMiddleStimLargeCorrect ./ max(max(SingleAnimalActionTraceMiddleStimLargeCorrect));
            GrandPopAbsActionRespMiddleStimCorLarge = SingleAnimalNormActionTraceMiddleStimLargeCorrect + GrandPopAbsActionRespMiddleStimCorLarge ;
            
            SingleAnimalActionTraceMiddleStimSmallCorrect= BehPhotoM(iAnimal).GrandSummary.AbsActionRasterSmallCorrect(2,:);
            SingleAnimalNormActionTraceMiddleStimSmallCorrect = SingleAnimalActionTraceMiddleStimSmallCorrect ./ max(max(SingleAnimalActionTraceMiddleStimLargeCorrect));
            GrandPopAbsActionRespMiddleStimCorSmall = SingleAnimalNormActionTraceMiddleStimSmallCorrect + GrandPopAbsActionRespMiddleStimCorSmall ;
            
            
            
            % tuning curve for different blocks
            SingleAnimalTunningStim= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimNoFold;
            SingleAnimalNormTunningStim = SingleAnimalTunningStim ./ max(max(SingleAnimalTunningStim));
            GrandPopNormBinStimNoFold = SingleAnimalNormTunningStim + GrandPopNormBinStimNoFold ;
            GrandPopNormBinStimNoFold1(c,:)=SingleAnimalNormTunningStim(1,:);
            GrandPopNormBinStimNoFold2(c,:)=SingleAnimalNormTunningStim(2,:);
            
            
            SingleAnimalTunningReward = BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardNoFold;
            
          
            
            if abs(max(max(SingleAnimalTunningReward))) < abs( min(min(SingleAnimalTunningReward))) % this is to avoid diving by a small number in one animal
           
                SingleAnimalTunningReward = SingleAnimalTunningReward + abs(min(min(SingleAnimalTunningReward)));      
            end

             SingleAnimalNormTunningRew = SingleAnimalTunningReward ./ max(max(SingleAnimalTunningReward));

            GrandPopNormBinRewardNoFold = SingleAnimalNormTunningRew + GrandPopNormBinRewardNoFold ;
            GrandPopNormBinRewardNoFold1(c,:)=SingleAnimalNormTunningRew(1,:);
            GrandPopNormBinRewardNoFold2(c,:)=SingleAnimalNormTunningRew(2,:);
            
            
            % tuning curve for corr/error
            SingleAnimalTunningStimCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrectErrorNoFold;
            SingleAnimalNormTunningStimCorrError = SingleAnimalTunningStimCorrError ./ max(max(SingleAnimalTunningStimCorrError));
            GrandPopNormBinStimNoFoldCorrError = SingleAnimalNormTunningStimCorrError + GrandPopNormBinStimNoFoldCorrError ;
            GrandPopNormBinStimNoFoldCorrError1(c,:)=SingleAnimalNormTunningStimCorrError(1,:);
            GrandPopNormBinStimNoFoldCorrError2(c,:)=SingleAnimalNormTunningStimCorrError(2,:);
            
            
            SingleAnimalTunningRewardCorrError= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardCorrectErrorNoFold;
            
           
              if abs(max(max(SingleAnimalTunningRewardCorrError))) < abs( min(min(SingleAnimalTunningRewardCorrError))) % this is to avoid diving by a small number in one animal
           
                SingleAnimalTunningRewardCorrError = SingleAnimalTunningRewardCorrError + abs(min(min(SingleAnimalTunningRewardCorrError)));      
              end
            
            
            SingleAnimalNormTunningRewardCorrError = SingleAnimalTunningRewardCorrError ./ max(max(SingleAnimalTunningRewardCorrError));

            GrandPopNormBinRewardNoFoldCorrError = SingleAnimalNormTunningRewardCorrError + GrandPopNormBinRewardNoFoldCorrError ;
            GrandPopNormBinRewardNoFoldCorrError1(c,:)=GrandPopNormBinRewardNoFoldCorrError(1,:);
            GrandPopNormBinRewardNoFoldCorrError2(c,:)=GrandPopNormBinRewardNoFoldCorrError(2,:);
            
            SingleAnimalBeep2AwayDACorrErr(1,:)= BehPhotoM(iAnimal).GrandSummary.Beep2DACorr;
            SingleAnimalBeep2AwayDACorrErr(2,:)= BehPhotoM(iAnimal).GrandSummary.BeepAwayDACorr;
            SingleAnimalBeep2AwayDACorrErr(3,:)= BehPhotoM(iAnimal).GrandSummary.Beep2DAErr;
            SingleAnimalBeep2AwayDACorrErr(4,:)= BehPhotoM(iAnimal).GrandSummary.BeepAwayDAErr;
            
            
            SingleAnimalBeep2AwayDACorrErrNorm = SingleAnimalBeep2AwayDACorrErr ./ max(max(SingleAnimalBeep2AwayDACorrErr));
            
            GrandPopBeep2AwayDACorError = SingleAnimalBeep2AwayDACorrErrNorm + GrandPopBeep2AwayDACorError ;
            
            SingleAnimalStim2AwayDACorrErr(1,:)= BehPhotoM(iAnimal).GrandSummary.Stim2DACorr;
            SingleAnimalStim2AwayDACorrErr(2,:)= BehPhotoM(iAnimal).GrandSummary.StimAwayDACorr;
            SingleAnimalStim2AwayDACorrErr(3,:)= BehPhotoM(iAnimal).GrandSummary.Stim2DAErr;
            SingleAnimalStim2AwayDACorrErr(4,:)= BehPhotoM(iAnimal).GrandSummary.StimAwayDAErr;
            
            
            SingleAnimalStim2AwayDACorrErrNorm = SingleAnimalStim2AwayDACorrErr ./ max(max(SingleAnimalStim2AwayDACorrErr));
            
            GrandPopStim2AwayDACorError = SingleAnimalStim2AwayDACorrErrNorm + GrandPopStim2AwayDACorError ;
            
            
            SingleBinStimCorrError(1:2,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimCorrect;
            SingleBinStimCorrError(3:4,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinStimError;
            
            SingleBinStimCorrErrorNorm = SingleBinStimCorrError ./ max(max(SingleBinStimCorrError));
            
            
            %GrandPopStimBin = SingleBinStimCorrErrorNorm + GrandPopStimBin ;
            
            GrandPopStimBin(:,:,c) = SingleBinStimCorrErrorNorm  ;
            
            SingleAnimalRew2AwayDACorrErr(1,:)= BehPhotoM(iAnimal).GrandSummary.Rew2DACorr;
            SingleAnimalRew2AwayDACorrErr(2,:)= BehPhotoM(iAnimal).GrandSummary.RewAwayDACorr;
            SingleAnimalRew2AwayDACorrErr(3,:)= BehPhotoM(iAnimal).GrandSummary.Rew2DAErr;
            SingleAnimalRew2AwayDACorrErr(4,:)= BehPhotoM(iAnimal).GrandSummary.RewAwayDAErr;
            
            
            SingleAnimalRew2AwayDACorrErrNorm = SingleAnimalRew2AwayDACorrErr ./ max(max(SingleAnimalRew2AwayDACorrErr));
            
            GrandPopRew2AwayDACorError = SingleAnimalRew2AwayDACorrErrNorm + GrandPopRew2AwayDACorError ;
            
            
            SingleBinRewCorrError(1:2,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardCorrect;
            SingleBinRewCorrError(3:4,:)= BehPhotoM(iAnimal).GrandSummary.PopNormBinRewardError;
            
            SingleBinRewCorrErrorNorm = SingleBinRewCorrError ./ max(max(SingleBinRewCorrError));
            
            
            %GrandPopRewBin = SingleBinRewCorrErrorNorm + GrandPopRewBin ;
            
            GrandPopRewBin(:,:,c) = SingleBinRewCorrErrorNorm  ;
            c=c+1;
            
        end
    end
    animalCount = animalCount +1; 
end

%%
figure

subplot(6,3,1); hold on
xlabel('Contrast')
ylabel('P(R)')
title( 'Psychometric curves')
set(gca,'TickDir','out','Box','off');



plot(StimAllowed,nanmean(PerBlock1),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
plot(StimAllowed,nanmean(PerBlock2),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)


legend('LargeRew@L','LargeRew@R','Location','southeast')


subplot(6,3,2); hold on
xlabel('Contrast')
ylabel('Norm RT')
title( 'Reaction Time')
set(gca,'TickDir','out','Box','off');


plot(StimAllowed,nanmean(RTBlock1),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)

plot(StimAllowed,nanmean(RTBlock2),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)



% subplot(5,3,4); hold on
% plot(StimAllowed,GrandPopNormBinStimNoFold(1,:)./ length(Animals),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinStimNoFold(2,:)./ length(Animals),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% title('Stimulus Align')
% xlabel('Contrast')


subplot(6,3,4); hold on % errorbars reflecting across sessions
errorbar(StimAllowed,GrandPopNormBinStimNoFold(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFold1) ./ sqrt(12),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinStimNoFold(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFold2) ./ sqrt(12),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
title('Stimulus Align')
xlabel('Contrast')




%subplot(5,3,5); hold on
% plot(StimAllowed,GrandPopNormBinRewardNoFold(1,:)./ length(Animals),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinRewardNoFold(2,:)./ length(Animals),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% xlabel('Contrast')
% title('Outcome Align')

subplot(6,3,5); hold on % errorbars reflecting across sessions
errorbar(StimAllowed,GrandPopNormBinRewardNoFold(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFold1) ./ sqrt(12),'color',[0.5 0.2 0.1],'LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinRewardNoFold(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFold2) ./ sqrt(12),'color',[1 0.6 0.2],'LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Contrast')


% subplot(5,3,7); hold on
% plot(StimAllowed,GrandPopNormBinStimNoFoldCorrError(1,:)./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinStimNoFoldCorrError(2,:)./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% title('Stimulus Align')
% xlabel('Contrast')

subplot(6,3,7); hold on
errorbar(StimAllowed,GrandPopNormBinStimNoFoldCorrError(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFoldCorrError1) ./ sqrt(14),'color','r','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinStimNoFoldCorrError(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinStimNoFoldCorrError2) ./ sqrt(14),'color','g','LineWidth',2,'Marker','o','MarkerSize',5)
set(gca,'TickDir','out','Box','off');
title('Stimulus Align')
xlabel('Contrast')




% subplot(5,3,8); hold on
% plot(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(1,:)./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)
% plot(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(2,:)./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)
% set(gca,'TickDir','out','Box','off');
% xlabel('Contrast')
% title('Outcome Align')

subplot(6,3,8); hold on

errorbar(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(1,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFoldCorrError1) ./ sqrt(14),'color','r','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(StimAllowed,GrandPopNormBinRewardNoFoldCorrError(2,:)./ length(Animals),...
    nanstd(GrandPopNormBinRewardNoFoldCorrError2) ./ sqrt(14),'color','g','LineWidth',2,'Marker','o','MarkerSize',5)

set(gca,'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Contrast')


subplot(6,3,10); hold on


errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopStimBin(1,:,:))'),nanstd(squeeze(GrandPopStimBin(1,:,:))')./ length(Animals),'--g','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopStimBin(2,:,:))'),nanstd(squeeze(GrandPopStimBin(2,:,:))')./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopStimBin(3,:,:))'),nanstd(squeeze(GrandPopStimBin(3,:,:))')./ length(Animals),'--r','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopStimBin(4,:,:))'),nanstd(squeeze(GrandPopStimBin(4,:,:))')./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)

set(gca,'TickDir','out','Box','off');
xlabel('Contrast')
title('Stimulus Align')


subplot(6,3,11); hold on


errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopRewBin(1,:,:))'),nanstd(squeeze(GrandPopRewBin(1,:,:))')./ length(Animals),'--g','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopRewBin(2,:,:))'),nanstd(squeeze(GrandPopRewBin(2,:,:))')./ length(Animals),'g','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopRewBin(3,:,:))'),nanstd(squeeze(GrandPopRewBin(3,:,:))')./ length(Animals),'--r','LineWidth',2,'Marker','o','MarkerSize',5)
errorbar(unique(abs(StimAllowed)),nanmean(squeeze(GrandPopRewBin(4,:,:))'),nanstd(squeeze(GrandPopRewBin(4,:,:))')./ length(Animals),'r','LineWidth',2,'Marker','o','MarkerSize',5)


set(gca,'TickDir','out','Box','off');
xlabel('Contrast')
title('Outcome Align')


subplot(6,3,12); hold on

plot(GrandPopStim2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopStim2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopStim2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopStim2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);


xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
ylim([-0.3 1])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');
title('Stimulus Align')
xlabel('Time (s)')
ylabel('Norm response')
legend('CorLarge','CorSmall','ErrLarge','ErrSmall')





subplot(6,3,13); hold on

plot(GrandPopRew2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopRew2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopRew2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopRew2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);

ylim([-0.8 1])

xlim([StartTime + (TimingVisualise(3,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(3,2)*sampleRate)]);

set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');
title('Outcome Align')
xlabel('Time (s)')



subplot(6,3,14); hold on

plot(GrandPopBeep2AwayDACorError(1,:)./ length(Animals),'g','LineWidth',2);

plot(GrandPopBeep2AwayDACorError(2,:)./ length(Animals),'--g','LineWidth',2);

plot(GrandPopBeep2AwayDACorError(3,:)./ length(Animals),'r','LineWidth',2);

plot(GrandPopBeep2AwayDACorError(4,:)./ length(Animals),'--r','LineWidth',2);

ylim([-0.3 0.7])

xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');

title('Beep Align')
xlabel('Time (s)')
ylabel('Norm response')

for c = 1:4
    
    
    subplot(6,3,15); hold on
    plot((GrandPopAbsStimResp(c,:) ./ length(Animals)),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Stimulus Align')


ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');


xlabel('Time (s)')
ylabel('Norm response')


for c = 1:4
    
    
    subplot(6,3,18); hold on
    plot((GrandPopAbsActionResp(c,:) ./ length(Animals)),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Action Align')


ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.2'},'TickDir','out','Box','off');

xlabel('Time (s)')
ylabel('Norm response')

for c = 2 % you could do the 1:4 and get all stimuli
    
    
    subplot(6,3,3); hold on
    plot(GrandPopAbsStimRespMiddleStimCorrect(c,:)./ length(Animals),'g','LineWidth',2);
    
    plot(GrandPopAbsStimRespMiddleStimError(c,:)./ length(Animals),'r','LineWidth',2);
    
end

subplot(6,3,3); hold on
title('middleSim,CorrError')

ylim([-0.3 0.5])


xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');

xlabel('Time (s)')
ylabel('Norm response')

for c = 2 % you could do the 1:4 and get all stimuli
    
    
    subplot(6,3,16); hold on
    plot(GrandPopAbsActionRespMiddleStimCorrect(c,:)./ length(Animals),'g','LineWidth',2);
    
    plot(GrandPopAbsActionRespMiddleStimError(c,:)./ length(Animals),'r','LineWidth',2);
    
end

subplot(6,3,16); hold on
title('middleSim,CorrError')

ylim([-0.3 0.5])


xlim([StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.2'},'TickDir','out','Box','off');

xlabel('Time (s)')
ylabel('Norm response')

subplot(6,3,6); hold on
title('middleSim,SmallLarge')

plot(GrandPopAbsStimRespMiddleStimCorLarge./ length(Animals),'g','LineWidth',2);
hold on
plot(GrandPopAbsStimRespMiddleStimCorSmall./ length(Animals),'--g','LineWidth',2);

xlim([3500 4900])
ylim([-0.3 0.5])


subplot(6,3,6); hold on
title('middleSim,SmallLarge')

plot(GrandPopAbsStimRespMiddleStimCorLarge./ length(Animals),'g','LineWidth',2);
hold on
plot(GrandPopAbsStimRespMiddleStimCorSmall./ length(Animals),'--g','LineWidth',2);

ylim([-0.3 0.5])

xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');

xlabel('Time (s)')
ylabel('Norm response')


subplot(6,3,17); hold on
title('middleSim,SmallLarge')

plot(GrandPopAbsActionRespMiddleStimCorLarge./ length(Animals),'g','LineWidth',2);
hold on
plot(GrandPopAbsActionRespMiddleStimCorSmall./ length(Animals),'--g','LineWidth',2);

ylim([-0.3 0.5])

xlim([StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.2'},'TickDir','out','Box','off');

xlabel('Time (s)')
ylabel('Norm response')

%%
figure

colorGray = [ 0.9 0.9 0.9
              0.8 0.8 0.8
              0.6 0.6 0.6
              0.5 0.5 0.5
              0.4 0.4 0.4
              0.2 0.2 0.2
               0 0 0  ];

for c = 1:7
        
    subplot(1,3,1); hold on
    plot(smooth((GrandPopStimLargeCorrect(c,:) ./ length(Animals)),100),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Stimulus Align')

ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');


xlabel('Time (s)')
ylabel('Norm response')

for c = 1:7
        
    subplot(1,3,2); hold on
    plot(smooth((GrandPopStimSmallCorrect(c,:) ./ length(Animals)),100),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Stimulus Align')

ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');


xlabel('Time (s)')
ylabel('Norm response')


for c = 1:7
        
    subplot(1,3,3); hold on
    plot(smooth((GrandPopStimLargeError(c,:) ./ length(Animals)),100),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Stimulus Align')

ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(1,1)*sampleRate) StartTime + (TimingVisualise(1,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(1,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');


xlabel('Time (s)')
ylabel('Norm response')


figure

TimingVisualise(2,:) = [-0.8 0.5]; 

for c = 1:7
        
    subplot(1,3,1); hold on
    plot(smooth((GrandPopActionLargeCorrect(c,:) ./ length(Animals)),100),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Action Align')

ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');


xlabel('Time (s)')
ylabel('Norm response')

for c = 1:7
        
    subplot(1,3,2); hold on
    plot(smooth((GrandPopActionSmallCorrect(c,:) ./ length(Animals)),100),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Action Align')

ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');


xlabel('Time (s)')
ylabel('Norm response')


for c = 1:7
        
    subplot(1,3,3); hold on
    plot(smooth((GrandPopActionLargeError(c,:) ./ length(Animals)),100),'color',colorGray(c,:),'LineWidth',2)
    
end


title('Action Align')

ylim([-0.3 1])


xlim([StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)])
set(gca, 'XTick', [StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)]);
set(gca, 'XTickLabel', {'0','0.8'},'TickDir','out','Box','off');


xlabel('Time (s)')
ylabel('Norm response')

