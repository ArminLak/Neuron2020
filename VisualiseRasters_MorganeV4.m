
% This code plots trial by trial rasters of photometry signal
% Morgane Moss June 2018

clear all
close all

%[48, 50,51]  coresponding to ALK068, 70 and 71
% session 5 of ALK068 is chosen to be shown in paper figure

animal_ID = 48
animal_name = animal_ID + 20;
load('BehPhotoM_Exp23')  % Database with data summary

%%
concatenate = 'n' % or n (show all trials from single animal)

smoothFactor = 100;


RTLimit = 10; % in s, excluding trials with RT longer than this

sample_rate = 12000;
downSample = 1200;
eventOnset = 3700;


% color range for plotting imagesc data
colorRange=[-10 12];
colorRange=[-5 10];


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


%-- START PLOTTING!!!

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
    
    TempStimData   = smooth2a(TempStimData,0,smoothFactor);
    TempActionData   = smooth2a(TempActionData,0,smoothFactor);
    TempRewardData   = smooth2a(TempRewardData,0,smoothFactor);
    
    
    RT = TempBehData(:,10) - TempBehData (:,13); % compute choice reaction time from action-stim interval
    RT(RT < 0) = nan;   %few trials with error negative RTs
    
    % exclude trials with RTs longer than predefiend RT limits
    TempBehData(RT > RTLimit,: )=[];
    TempStimData(RT > RTLimit,: )=[];
    TempActionData(RT > RTLimit,: )=[];
    TempRewardData(RT > RTLimit,: )=[];
    RT(RT > RTLimit)=[];
    
    TempBehData(isnan(RT),: )=[];
    TempStimData(isnan(RT),: )=[];
    TempActionData(isnan(RT),: )=[];
    TempRewardData(isnan(RT),: )=[];
    RT(isnan(RT))=[];
    
    TempBehData(:,7)= RT;
    [i sortingIndex] = sort(RT);
    
    largeRew = sort([(intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==-1), find(TempBehData(:,8)==1))); (intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==1), find(TempBehData(:,8)==2)))]);
    
    smallRew = setdiff([1:length(TempBehData)],largeRew)';
    
    correct = find(TempBehData(:,9)==1);
    error = find(TempBehData(:,9)==0);
    
    
    
    %------------------------------- plot psychometric curve------------------
    %
    
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
        
        largeonLTrials = find(TempBehData(:,8)==1); % large on left
        largeonRTrials = find(TempBehData(:,8)==2); % large on right
        
        
        for istim = TempStimz'
            largeLPerformance(c) = nanmean(TempBehData(intersect(find(TempBehData(:,2)==istim), largeonLTrials),3));
            largeRPerformance(c) = nanmean(TempBehData(intersect(find(TempBehData(:,2)==istim), largeonRTrials),3));
            c=c+1;
        end
        
        plot(TempStimz, largeLPerformance,'k', 'LineWidth', 1.5)
        hold on;
        plot(TempStimz, largeRPerformance,'color', [0.5 0.5 0.5], 'LineWidth', 1.5) %plot two lines if concatenated
        
        
    end
    
    %  c=c+1;
    %     end
    
    xticks([min(TempStimz) 0 max(TempStimz)])
    xlabel('Contrast')
    yticks([-1 0 1])
    yticklabels([0 0.5 1])
    ylabel('% Right', 'FontWeight', 'bold')
    
    % -------------- plot title ------------------------------------
    
    text(TempStimz(end), 1.2, ['ALK0' num2str(animal_name) '  Session ' num2str(iSession)], 'FontWeight', 'bold', 'FontSize', 10);
    
    
    % ---------- raster for all trials in session in order of RT ---------
    %     figure; hold on;
    subplot(2, 4, 5);
    
    TempBehData2sort = TempBehData;
    TempBehData2sort(:,2)=abs(TempBehData2sort(:,2));
    TempBehData2sort(error,2)=-TempBehData2sort(error,2);
    
    [i j]= sortrows(TempBehData2sort,[9,2,7]);
    
    imagesc((TempStimData(j, 1:6700)),colorRange)
    
    colormap('bluewhitered')
    
  
       trace = 1;
        
        for ievent=RT(j)'
            
            H=line([downSample*ievent+eventOnset,downSample*ievent+eventOnset+20], [trace, trace]);
            set(H,'color',[0 0 0],'LineWidth',3)
            
            trace  = trace + 1;
        end
        
        trace = 1;
        
        StimOutcomeInt = TempBehData(:,14)-TempBehData(:,13); 
        for ievent=StimOutcomeInt(j)'
            
            H=line([downSample*ievent+eventOnset,downSample*ievent+eventOnset+20], [trace, trace]);
            set(H,'color',[0 1 0],'LineWidth',3)
            
            trace  = trace + 1;
        end
        
    line([eventOnset eventOnset], [max(sortingIndex) min(sortingIndex)], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    text(100, -20, 'Stim', 'FontWeight', 'bold', 'FontSize', 8); %stim line label
    
    ylabel('Trials', 'FontWeight', 'bold')
    xlim([3500 6700])
    xticks([1 RTLimit*downSample])
    xticklabels([0-(200/(sample_rate/downSample)) RTLimit*(downSample)])
    xlabel('Time (s)', 'FontWeight', 'bold')
    
    
    %-- STIM,ACTION and OUTCOME RESPONSE RASTERSSSS ---------------------------------------
    
    % ------------- 1.1 raster for stim reponse, max contrast (L+R)
    
    %     figure; hold on;
    
    
    if length(StimzAbs)>3
        subplot(5, 8, 3);
        
        tempIndex = intersect(find(abs(TempBehData(:,2))==StimzAbs(length(StimzAbs))),correct);
        
        imagesc(TempStimData(tempIndex, (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        title('Stimulus')
        
        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(end))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end)) ' Contrast Trials'], 'FontWeight', 'bold')
        
        % ------------- 1.2 raster for action, 0.5 contrast (L+R)
        subplot(5, 8, 4);
        
        imagesc(TempActionData(tempIndex, (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        title('Action')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(4))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        % ------------- 1.3 raster for outcome, 0.5 contrast (L+R)
        
        subplot(5, 8, 5);
        
        imagesc(TempRewardData(tempIndex, (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
        colormap('bluewhitered')
        title('Reward')
        
        line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(4))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        
    end
    
    % ------------- 2.1 raster for sitm response, 2nd max contrast (L+R)
    
    if length(StimzAbs)>1
        subplot(5, 8, 11);
        
        tempIndex = intersect(find(abs(TempBehData(:,2))==StimzAbs(end-1)),correct);
        
        imagesc(TempStimData(tempIndex, (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        
        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end-1)) ' Contrast Trials'], 'FontWeight', 'bold')
        
        % ------------- 2.2 raster for action, 0.25 contrast (L+R)
        
        subplot(5, 8, 12);
        
        imagesc(TempActionData(tempIndex, (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(3))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        % ------------- 2.3 raster for outcome, 0.25 contrast (L+R)
        
        subplot(5, 8, 13);
        
        imagesc(TempRewardData(tempIndex, (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(3))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        
    end
    
    
    % ------------- 3.1 raster for stim response, 3rd max contrast (L+R)
    
    if length(StimzAbs)>2
        subplot(5, 8, 19);
        
        tempIndex = intersect(find(abs(TempBehData(:,2))==StimzAbs(end-2)),correct);
        
        
        imagesc(TempStimData(tempIndex, (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-2)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end-2)) ' Contrast Trials'], 'FontWeight', 'bold')
        
        % ------------- 3.2 raster for action, 0.25 contrast (L+R)
        
        subplot(5, 8, 20);
        
        imagesc(TempActionData(tempIndex, (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(2))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        % ------------- 3.3 raster for action, 0.25 contrast (L+R)
        
        subplot(5, 8, 21);
        
        imagesc(TempRewardData(tempIndex, (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(2))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        
        
    end
    
    % ----------- 4.1 raster for stim response, 0 contrast
    
    if length(StimzAbs)>3
        subplot(5, 8, 27);
        
        tempIndex = intersect(find(abs(TempBehData(:,2))==StimzAbs(end-3)),correct);
        
        imagesc(TempStimData(tempIndex, (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData(abs(TempBehData(:,2))==StimzAbs(end-3)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticks([ ])
        xticklabels([ ])
        ylabel([num2str(StimzAbs(end-3)) ' Contrast Trials'], 'FontWeight', 'bold')
        
        % ------------- 4.2 raster for action, 0 contrast
        
        subplot(5, 8, 28);
        
        imagesc(TempActionData(tempIndex, (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(1))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        
        % --------- 4.3 raster for reward response, 0 contrast
        
        subplot(5, 8, 29);
        
        imagesc(TempRewardData(tempIndex, (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(1))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
    end
    
    
    % ------------- 5. line plot; avg stim response over time by contrast
    
    subplot(5, 8, 35);
    
    c=1;
    for istim = StimzAbs'
        
         tempIndex = intersect(find(abs(TempBehData(:,2))==istim),correct);

        hold on
        h = plot(nanmean(TempStimData(tempIndex, (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample)))), 'LineWidth', 2)
        
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
    
    
    % ------------- 5. line plot; avg stim response over time by contrast
    
    subplot(5, 8, 36);
    
    c=1;
    for istim = StimzAbs'
        hold on
         tempIndex = intersect(find(abs(TempBehData(:,2))==istim),correct);

        h = plot(nanmean(TempActionData(tempIndex, (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample)))), 'LineWidth', 2)
        
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
    
    
    % ------------- 5. line plot; avg rwd response over time by contrast
    
    subplot(5, 8, 37);
    
    c=1;
    for istim = StimzAbs'
        hold on
                 tempIndex = intersect(find(abs(TempBehData(:,2))==istim),correct);

        h = plot(nanmean(TempRewardData(tempIndex, (eventOnset-abs(rStart*downSample)):(eventOnset+abs(rStop*downSample)))), 'LineWidth', 2)
        
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
    
    
    
    %-- STIM RESPONSE RASTERS 2 (BROKEN BY LARGE/SMALL/NO REWARD) -------------------
    
    
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
    
    
    %-- REWARD RESPONSE RASTERS 2 (BROKEN BY LARGE/SMALL/NO REWARD) -------------------
    %
    
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
    
    line([abs(rStart*downSample) abs(rStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
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


