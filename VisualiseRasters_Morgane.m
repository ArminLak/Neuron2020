

clear all
close all

%[48, 50,51]  coresponding to ALK068, 70 and 71

animal_ID = 50
animal_name = animal_ID + 20;
load('BehPhotoM_Exp23')


concatenate = 'y' % or n (show all trials from single animal)


RTLimit = 10; % in s, excluding trials with RT longer than this

sample_rate = 12000
downSample = 1200
eventOnset = 3700

colorRange=[-10 15]

% ------------------ start stop times for task events ------------------

sStart = -0.2; %stimulus
sStop = 0.8;

aStart = -0.5; %action
aStop = 0.5;

rStart = -0.2; %reward
rStop = 0.8;

%---------------- colors for plotting-------------------------------------

colorGray4 = [0.8 0.8 0.8
    0.6 0.6 0.6
    0.3 0.3 0.3
    0 0 0];
colorGreen = [0 1 0
    0 0.8 0
    0 0.6 0
    0  0.3 0];
colorRed = [1 0 0
    0.8 0 0
    0.6 0  0
    0.3 0 0];


%%

    sessionz = 1:length(BehPhotoM(animal_ID).Session);
    
    
if concatenate == 'y'
    
    TempBehData=[];
    TempBeepData=[];
    TempStimData=[];
    TempActionData=[];
    TempRewardData=[];

        for iSession = sessionz

            % for each of these, the event is at column 3700

            TempBehData = [TempBehData; BehPhotoM(animal_ID).Session(iSession).TrialTimingData];
            TempBeepData   = [TempBeepData; BehPhotoM(animal_ID).Session(iSession).NeuronBeep];
            TempStimData   = [TempStimData; BehPhotoM(animal_ID).Session(iSession).NeuronStim];            
            TempActionData = [TempActionData; BehPhotoM(animal_ID).Session(iSession).NeuronAction];
            TempRewardData = [TempRewardData; BehPhotoM(animal_ID).Session(iSession).NeuronReward];
            
        end

        TempStimz      = unique(TempBehData(:,2));
        StimzAbs       = unique(abs(TempBehData(:,2)));

        sessionz = 1;

elseif concatenate == 'n'

    sessionz = 1:length(BehPhotoM(animal_ID).Session);

end

%% START PLOTTING!!! 

for iSession = sessionz
    
    % for each of these, the event is at column 3700
    if concatenate == 'n'
        
        TempBehData    = BehPhotoM(animal_ID).Session(iSession).TrialTimingData;
        TempBeepData   = BehPhotoM(animal_ID).Session(iSession).NeuronBeep;
        TempStimData   = BehPhotoM(animal_ID).Session(iSession).NeuronStim;
        TempActionData = BehPhotoM(animal_ID).Session(iSession).NeuronAction;
        TempRewardData = BehPhotoM(animal_ID).Session(iSession).NeuronReward;
        TempStimz      = unique(TempBehData(:,2));
        StimzAbs       = unique(abs(TempBehData(:,2)));
    
    end
    
    
    TempStimData   = smooth2a(TempStimData,0,100);
    TempActionData   = smooth2a(TempActionData,0,100);
    TempRewardData   = smooth2a(TempRewardData,0,100);
    
    % here you can use imagesc to make colourful plots
    % for each session lets make the following subplots
    
    
    
    RT = TempBehData(:,10) - TempBehData (:,13); % compute choice reaction time from action-stim interval
    [i sortingIndex] = sort(RT);
    
%     largeRew = sort(find(TempBehData(:,8)==1 & TempBehData(:,9)==1 & 
    
    largeRew = sort([(intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==-1), find(TempBehData(:,8)==1))); (intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==1), find(TempBehData(:,8)==2)))]);
    smallRew = setdiff([1:length(TempBehData)],largeRew)';
    
    %------------------------------- plot psychometric curve------------------
    %%
    
    f = figure('Position', [300 200 800 900]); hold on
    subplot(2, 4, 1);
    
    c = 1;
    
    if concatenate == 'n'      
        for istim = TempStimz'
            performance(c) = nanmean (TempBehData(TempBehData(:,2)==istim,3));
            c=c+1;
        end
        plot(TempStimz, performance,'k', 'LineWidth', 1.5) % plot one line if not concatenated
        
        
    elseif concatenate == 'y'
        
        largeRewTrials = [find(TempBehData(:,8)==1 & TempBehData(:,3)==-1); find((TempBehData(:,8)==2 & TempBehData(:,3)==1))];
        smallRewTrials = [find(TempBehData(:,8)==1 & TempBehData(:,3)==1); find((TempBehData(:,8)==2 & TempBehData(:,3)==-1))];


        for istim = TempStimz' 
            largePerformance(c) = nanmean(TempBehData(intersect(find(TempBehData(:,2)==istim), LTrials),3))
            smallPerformance(c) = nanmean(TempBehData(intersect(find(TempBehData(:,2)==istim), STrials),3))
            c=c+1;
        end
        
        plot(TempStimz, largePerformance,'k', 'LineWidth', 1.5)
        hold on;
        plot(TempStimz, smallPerformance,'color', [0.5 0.5 0.5], 'LineWidth', 1.5) %plot two lines if concatenated 
    
    end
    
    
    xticks([min(TempStimz) 0 max(TempStimz)])
    xlabel('Contrast')
    yticks([-1 0 1])
    yticklabels([0 0.5 1])
    ylabel('% Right', 'FontWeight', 'bold')
    
    % -------------- plot title ------------------------------------
    
    text(TempStimz(end), 1.2, ['ALK0' num2str(animal_name) '  Session ' num2str(iSession)], 'FontWeight', 'bold', 'FontSize', 10);
    
    
    %%
    
    % ---------- raster for all trials in session in order of RT ---------
    %     figure; hold on;
    subplot(2, 4, 5);
    
    imagesc((TempStimData(sortingIndex, (eventOnset-200):end)),colorRange)
    
    colormap('bluewhitered')
    
    %
    %     a = 1
    %     for i = sortingIndex
    %         line([200+RT(i) 200+RT(i)], [i i], 'Color', 'black', 'LineWidth', 1.5); % create tiny line for each session showing time of action
    %
    % %         if i == sortingIndex(1)
    % %             text(200, 150,
    %
    %     end
    %
    %
    %     line ( [ 200+RT(432)*downSample 200+RT(432)*downSample ], [sortingIndex(200) sortingIndex(400)], 'color', 'black', 'LineWidth', 1.5);
    line([200 200], [max(sortingIndex) min(sortingIndex)], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    text(100, -20, 'Stim', 'FontWeight', 'bold', 'FontSize', 8); %stim line label
    
    ylabel('Trials', 'FontWeight', 'bold')
    %     xticks([1:1200:4500])
    %     xticklabels([0 1 2 3])
    xticklabels([TempStimData./1200])
    xlabel('Time (s)', 'FontWeight', 'bold')
    
    %% STIM RESPONSE RASTERSSSS ---------------------------------------
    
    % ------------- 1. raster for stim reponse, max contrast (L+R)
    
    %     figure; hold on;
    
    
    if length(StimzAbs)>3
        subplot(5, 8, 3);

        imagesc(TempStimData(abs(TempBehData(:,2))==StimzAbs(length(StimzAbs)), (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        title('Stimulus')

        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(end))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line

        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end)) ' Contrast Trials'], 'FontWeight', 'bold')
    
    end
    
    % ------------- 2. raster for sitm response, 2nd max contrast (L+R)
    
    if length(StimzAbs)>1
        subplot(5, 8, 11);
        
        imagesc(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-1), (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        
        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end-1)) ' Contrast Trials'], 'FontWeight', 'bold')
        
    end
    
    
    % ------------- 3. raster for stim response, 3rd max contrast (L+R)
    
    if length(StimzAbs)>2
        subplot(5, 8, 19);
        
        imagesc(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-2), (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-2)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end-2)) ' Contrast Trials'], 'FontWeight', 'bold')
        
    end
    
    % ----------- 4. raster for stim response, 0 contrast
    
    if length(StimzAbs)>3
        subplot(5, 8, 27);
        
        imagesc(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-3), (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-3)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end-3)) ' Contrast Trials'], 'FontWeight', 'bold')
        
    end
    
    
    % ------------- 5. line plot; avg stim response over time by contrast
    
    subplot(5, 8, 35);
    
    c=1;
    for istim = StimzAbs'
        hold on
        h = plot(nanmean(TempStimData(abs(TempBehData(:,2))==istim, (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample)))), 'LineWidth', 2)
        
        if length(StimzAbs)==4
            set(h, 'color', colorGray4(c,:));
        elseif length(StimzAbs)==3
            set(h, 'color', colorGray3(c,:));
        elseif length(StimzAbs)==2
            set(h, 'color', colorGray2(c,:));
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
    
    xlim([0 abs(sStop*downSample)+abs(sStart*downSample)])
    xticks([1 abs(sStop*downSample)+abs(sStart*downSample)])
    xticklabels([sStart sStop])
    xlabel('Time (s)', 'FontWeight', 'bold')
    ylabel('{\Delta} F / F')
    
    
    %% ACTION RASTERS ------------------------------------------------------
    
    % ------------- 1. raster for action, 0.5 contrast (L+R)
    
%         figure; hold on;
    
    if length(StimzAbs)>3
        subplot(5, 8, 4);
        
        imagesc(TempActionData(abs(TempBehData(:,2))==StimzAbs(4), (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        title('Action')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(4))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
    
    end
    
    % --------- 2. raster for action, 0.25 contrast (L+R)
    
    if length(StimzAbs)>2
        subplot(5, 8, 12);
        
        imagesc(TempActionData(abs(TempBehData(:,2))==StimzAbs(3), (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(3))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
    
    end
    
    % --------- 3. raster for action, 0.12 contrast (L+R)
    
    if length(StimzAbs)>1
        subplot(5, 8, 20);
        
        imagesc(TempActionData(abs(TempBehData(:,2))==StimzAbs(2), (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(2))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
    
    end
    
    % --------- 4. raster for action, 0 contrast
    
    subplot(5, 8, 28);
    
    imagesc(TempActionData(abs(TempBehData(:,2))==StimzAbs(1), (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
    colormap('bluewhitered')
    
    line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(1))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
        xticklabels([ ])
        yticklabels([ ])
    
    
    % ------------- 5. line plot; avg stim response over time by contrast
    
    subplot(5, 8, 36);
    
    c=1;
    for istim = StimzAbs'
        hold on
        h = plot(nanmean(TempActionData(abs(TempBehData(:,2))==istim, (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample)))), 'LineWidth', 2)
        
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
    
    xlim([0 abs(aStop*downSample)+abs(aStart*downSample)])
    xticks([1 abs(aStop*downSample)+abs(aStart*downSample)])
    xticklabels([aStart aStop])
    xlabel('Time (s)', 'FontWeight', 'bold')
    
    %% REWARD RESPONSE RASTERS ---------------------------------------------
    
    % ------------- 1. raster for reward response, 0.5 contrast (L+R)
    
    %     figure; hold on;
    
    if length(StimzAbs)>3
        subplot(5, 8, 5);
        
        imagesc(TempRewardData(abs(TempBehData(:,2))==StimzAbs(4), (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
        colormap('bluewhitered')
        title('Reward')
        
        line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(4))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
    end
    
    % --------- 2. raster for reward response, 0.25 contrast (L+R)
    
    if length(StimzAbs)>2
        subplot(5, 8, 13);
        
        imagesc(TempRewardData(abs(TempBehData(:,2))==StimzAbs(3), (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(3))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
    end
    
    % --------- 3. raster for reward response , 0.12 contrast (L+R)
    
    if length(StimzAbs)>1
        subplot(5, 8, 21);
        
        imagesc(TempRewardData(abs(TempBehData(:,2))==StimzAbs(2), (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(2))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
    end
    
    % --------- 4. raster for reward response, 0 contrast
    
    subplot(5, 8, 29);
    
    imagesc(TempRewardData(abs(TempBehData(:,2))==StimzAbs(1), (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
    colormap('bluewhitered')
    
    line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(1))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
        xticklabels([ ])
        yticklabels([ ])
       
    % ------------- 5. line plot; avg rwd response over time by contrast
    
    subplot(5, 8, 37);
    
    c=1;
    for istim = StimzAbs'
        hold on
        h = plot(nanmean(TempRewardData(abs(TempBehData(:,2))==istim, (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample)))), 'LineWidth', 2)
        
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
    
    xlim([0 abs(rStop*downSample)+abs(rStart*downSample)])
    xticks([1 abs(rStop*downSample)+abs(rStart*downSample)])
    xticklabels([rStart rStop])
    xlabel('Time (s)', 'FontWeight', 'bold')
    
    
    %% STIM RESPONSE RASTERS 2 (BROKEN BY LARGE/SMALL/NO REWARD) -------------------
    
    
    % ------------- 1. raster for stim reponse, large reward trials
%         figure; hold on;
    subplot(4, 8, 7);
    
    imagesc(TempStimData(largeRew, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample)), colorRange)
    colormap('bluewhitered')
    
    line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
    title('Stimulus')
    xticklabels([ ])
    ylabel('Large Reward Trials', 'FontWeight', 'bold')
    
    % ------------- 2. raster for stim reponse, small reward trials
    
    subplot(4, 8, 15);
    
    imagesc(TempStimData(smallRew, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample)), colorRange)
    colormap('bluewhitered')
    
    line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
    xticklabels([ ])
    ylabel('Small Reward Trials', 'FontWeight', 'bold')
    

    % ------------- 3. raster for stim reponse, NO reward trials
    
    subplot(4, 8, 23);
    
    imagesc(TempStimData((TempBehData(:,9)==0), (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
    
    colormap('bluewhitered')
    
    line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
    xticklabels([ ])
    ylabel('No Reward Trials', 'FontWeight', 'bold')
    
    
    
    % -------------- 4. summary avg reward response large/small/error

        subplot(5, 8, 39);

        plot(nanmean(TempStimData(largeRew, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample))),'color',colorGreen(3,:),'LineWidth',2)
    hold on; 
        plot(nanmean(TempStimData(smallRew, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample))),'color',colorGreen(1,:),'LineWidth',2)
    hold on; 
        plot(nanmean(TempStimData(TempBehData(:,9)==0, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample))),'color',colorRed(2,:),'LineWidth',2)

    legend('Large reward','Small reward', 'No reward','location','best')
    
    xlim([0 abs(sStop*downSample)+abs(sStart*downSample)])
    xticks([1 abs(sStop*downSample)+abs(sStart*downSample)])
    xticklabels([sStart sStop])
    ylabel('{\Delta} F / F')
    xlabel('Time (s)', 'FontWeight', 'bold')
    
    
    %% REWARD RESPONSE RASTERS 2 (BROKEN BY LARGE/SMALL/NO REWARD) -------------------
    
    
    % ------------- 1. raster for stim reponse, large reward trials
%     figure; hold on;
    subplot(4, 8, 8);
 
    imagesc(TempRewardData(largeRew, eventOnset-abs(rStart*downSample):eventOnset+abs(rStop*downSample)), colorRange)
    colormap('bluewhitered')
    
    line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
    title('Reward')
    xticklabels([ ])
    yticklabels([ ])
    
    
    % ------------- 2. raster for stim reponse, small reward trials
    
    subplot(4, 8, 16);
    
    imagesc(TempRewardData(smallRew, eventOnset-abs(rStart*downSample):eventOnset+abs(rStop*downSample)), colorRange)
    colormap('bluewhitered')
    
    line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
    xticklabels([ ])
    yticklabels([ ])
    
    % ------------- 3. raster for stim reponse, NO reward trials
    
    subplot(4, 8, 24);
    
    imagesc(TempRewardData((TempBehData(:,9)==0), (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
    
    colormap('bluewhitered')
    
    line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line

    xticklabels([ ])
    yticklabels([ ])
    
    % ---------------- 4. summary avg reward response large/small/error

        subplot(5, 8, 40);
    %plot large rewards
    
        plot(nanmean(TempRewardData(largeRew, eventOnset-abs(rStart*downSample):eventOnset+abs(rStop*downSample))),'color',colorGreen(3,:),'LineWidth',2)
    hold on; 
        plot(nanmean(TempRewardData(smallRew, eventOnset-abs(rStart*downSample):eventOnset+abs(rStop*downSample))),'color',colorGreen(1,:),'LineWidth',2)
    hold on; 
    
    % plot no rewards:
    
    plot(nanmean(TempRewardData(TempBehData(:,9)==0, eventOnset-abs(rStart*downSample):eventOnset+abs(rStop*downSample))),'color',colorRed(2,:),'LineWidth',2)

        legend('Large reward','Small reward', 'No reward','location','best')
    
    xlim([0 abs(rStop*downSample)+abs(rStart*downSample)])
    xticks([1 abs(rStop*downSample)+abs(rStart*downSample)])
    xticklabels([rStart rStop])
    xlabel('Time (s)', 'FontWeight', 'bold')
    
    
    
    %%
    
    
    
    % aling to stimulus, sorted based on RT
    
    % subplots for different levels of stimuli (average L and R)
    
    % subplots (outcome aligned) for large reward, small reward and no
    % reward
    
    % so you will write imagesc
    
    % then call colormap(bluewhitered)
end


