clear all
close all

%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 51
load('BehPhotoM_Exp23')

RTLimit = 10; % in s, excluding trials with RT longer than this

sample_rate = 12000
downSample = 1200
eventOnset = 3700


% ------------------ start stop times for task events ------------------

sstart = -0.2; %stimulus 
sstop = 0.8;

astart = -0.5; %action
astop = 0.5;

rstart = -0.2; %reward
rstop = 0.8;


%%


sessionz = 1:length(BehPhotoM(animal_ID).Session);


for iSession = sessionz
    
    % for each of these, the event is at column 3700
    
    TempBehData   = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
    TempBeepData  = BehPhotoM(animal_ID).Session(iSession).NeuronBeep; 
    TempStimData  = BehPhotoM(animal_ID).Session(iSession).NeuronStim;  
    TempRewardData= BehPhotoM(animal_ID).Session(iSession).NeuronReward;  
    TempStimz     = unique(TempBehData(:,2));
    

    % here you can use imagesc to make colourful plots
    % for each session lets make the following subplots
    
    f = figure('Position', [300 200 800 900]); hold on
    
    RT = TempBehData(:,10) - TempBehData (:,13); % compute choice reaction time from action-stim interval
   [i sortingIndex] = sort(RT);
    %------------------------------- plot psychometric curve------------------
   %% 
    subplot(2, 3, 1);
    
    c = 1;
    for istim = TempStimz'
    
        performance(c) = nanmean (TempBehData(TempBehData(:,2)==istim,3));
%         RT(c) = nanmean (ReactionTime(TrialTimingData(:,2)==istim));
    
        c=c+1; 
    end
    
    plot(TempStimz, performance,'k', 'LineWidth', 1.5)
    
    xticks([min(TempStimz) 0 max(TempStimz)])
    xlabel('Contrast')
    yticks([-1 0 1])
    yticklabels([0 0.5 1])
    ylabel('% Right')
    
    % -------------- plot title ------------------------------------
    
    text(TempStimz(end), 1.2, 'ALK071', 'FontWeight', 'bold', 'FontSize', 10);
    
    
%%
    
    % ---------- raster for all trials in session in order of RT ---------
%     figure; hold on; 
    subplot(2, 3, 4);
    
    imagesc((TempStimData(sortingIndex, (eventOnset-200):4500)))
    
    colormap('bluewhitered')
%     
%     for i = sortingIndex
%         line([200+RT(i)*downSample 200+RT(i)*downSample], [i i], 'Color', 'black', 'LineWidth', 1.5); % create tiny line for each session showing time of action 
%         
% %         if i == sortingIndex(1)
% %             text(200, 150, 
%     end
    
	line([200 200], [max(sortingIndex) min(sortingIndex)], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    text(100, -20, 'Stim', 'FontWeight', 'bold', 'FontSize', 8); %stim line label
    
    ylabel('Trials', 'FontWeight', 'bold')
    xticks([1:1200:4500])
%     xticklabels([0 1 2 3])
    xticklabels([TempStimData./1200])
    xlabel('Time (s)', 'FontWeight', 'bold')
    
%% STIM RESPONSE RASTERSSSS ---------------------------------------

    % ------------- 1. raster for stim reponse, -0.5 contrast
%     figure; hold on;
    subplot(4, 7, 4);

    imagesc((TempStimData((TempBehData(:,2)==0.5), (eventOnset-abs(sstart*downSample)):(eventOnset+abs(sstop*downSample)))))
    colormap('bluewhitered')

    xticklabels(['-.2' '0' '0.8'])
    xlabel('Time (s)', 'FontWeight', 'bold')     
    ylabel('+0.5 Contrast Trials', 'FontWeight', 'bold')
        
    
    
    % ------------- 2. raster for sitm response, 0 contrast 
    
    subplot(4,7, 11);
    
    imagesc((TempStimData((TempBehData(:,2)==0), (eventOnset-abs(sstart*downSample)):(eventOnset+abs(sstop*downSample)))))
    colormap('bluewhitered')
    
    xticklabels(['-.2' '0' '0.8'])
    xlabel('Time (s)', 'FontWeight', 'bold')     
    ylabel('0 Contrast Trials', 'FontWeight', 'bold')
    
    
    % ------------- 3. raster for stim response, 0.5 contrast 
    
    subplot(4,7, 18);
    
    imagesc((TempStimData((TempBehData(:,2)==-0.5), (eventOnset-abs(sstart*downSample)):(eventOnset+abs(sstop*downSample)))))
    colormap('bluewhitered')
    
    xticklabels(['-.2' '0' '0.8'])
    xlabel('Time (s)', 'FontWeight', 'bold')     
    ylabel('-0.5 Contrast Trials', 'FontWeight', 'bold')
    
    
    % ------------- 4. line plot; avg stim response over time by contrast
    
%     xlim([0 (sstop-sstart)* sample_rate]/downsampleScale)
%         xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
%     xticklabels([start:1:stop])
%xticklabels({'-2','-1','0','1','2'})

c=1;
for istim = StimzAbs
    h = plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim,:)),'color',colorGray4(c,:),'LineWidth',2)
    
    if length(StimzAbs)==4
        set(h, 'color', colorGray4(c,:))
    elseif length(StimzAbs)==3
        set(h, 'color', colorGray3(c,:))
    elseif length(StimzAbs)==2
        set(h, 'color', colorGray2(c,:))
    end
    
    c=c+1;
end
    
if length(StimzAbs)==4
    l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),num2str(StimzAbs(4)),'location','best')
elseif length(StimzAbs)==3
    l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),'location','best')
elseif length(StimzAbs)==2
    l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), 'location', 'best')
end
    


%% STIM RESPONSE RASTERS 2 (BROKEN BY CORRECT/ERROR) -------------------


    % ------------- 1. raster for stim reponse, correct trials
%     figure; hold on;
    subplot(2, 7, 5);

    imagesc((TempStimData((TempBehData(:,9)==1), (eventOnset-abs(sstart*downSample)):(eventOnset+abs(sstop*downSample)))))
    colormap('bluewhitered')

    xticklabels(['-.2' '0' '0.8'])
    xlabel('Time (s)', 'FontWeight', 'bold')     
    ylabel('Correct Trials', 'FontWeight', 'bold')
        

        % ------------- 2. raster for stim reponse, error trials
        
    subplot(2, 7, 12);

    imagesc((TempStimData((TempBehData(:,9)==0), (eventOnset-abs(sstart*downSample)):(eventOnset+abs(sstop*downSample)))))
    colormap('bluewhitered')

    xticklabels(['-.2' '0' '0.8'])
    xlabel('Time (s)', 'FontWeight', 'bold')     
    ylabel('Error Trials', 'FontWeight', 'bold')
    
    
    
    % --------------- 3. summary avg stim response 
%     
%     subplot(7,2,10); hold on
% title ('Stimulus aligned')
% c=1;
% 
% for istim = StimzAbs
%     
%     plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==1,:)),'color',colorGreen(c,:),'LineWidth',2)
%     plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==0,:)),'color',colorRed(c,:),'LineWidth',2)
%     
%     c=c+1;
%     
% end
% 

%%


    
    % aling to stimulus, sorted based on RT
    
    % subplots for different levels of stimuli (average L and R)
    
    % subplots (outcome aligned) for large reward, small reward and no
    % reward
    
    % so you will write imagesc
    
    % then call colormap(bluewhitered)
end