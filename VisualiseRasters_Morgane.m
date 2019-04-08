
% This code plots trial by trial rasters of photometry signal
% Morgane Moss June 2018

clear all
% close all

% things to do:
%1. tickdirection out
%2. render as painter
% in the sorting based on CrrErrContRT, plot a thick green line on the
% boarders between trials


%[48, 50,51]  coresponding to ALK068, 70 and 71
% session 5 of ALK068 is chosen to be shown in paper figure

animal_ID = 48
animal_name = animal_ID + 20;
load('BehPhotoM_Exp23_VTA')  % Database with data summary

%%
concatenate = 'y' % y or n (show all trials from single animal)

selectStimulus2Plot = 'y' %(this only works if concatenate = 'y')

sortDesign = 'CrrErrContRT'; % 'RT' or 'CrrErrContRT' sort based on RT or CorrectError Contrst and RTs

smoothFactor = 100;  % for smoothing the calcium data 

RTLimit = 3; % in s, excluding trials with RT longer than this

sample_rate = 12000;
downSample = 1200;
eventOnset = 3700;  % in the saved matrix, the ev


% color range for plotting imagesc data
colorRange=[-2 10];

if animal_ID ==48
colorRange=[-2 10]; end

if animal_ID == 51
    colorRange=[-5 20]; end

% ------------------ start stop times for task events in second ------------------

sStart = -0.2; %stimulus
sStop = 0.8;

aStart = -0.8; %action
aStop = 0.2;

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
rewardColor = [0 0.6 0
    0 1 0
    0.8 0 0];
% ------------------------------------------------------------------------

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
    
    
    if selectStimulus2Plot == 'y'  % making plots only for a subset of trials
        stim2plot = [-0.25 0.25];
        trialSubset=find(ismember(TempBehData(:,2),stim2plot));
        
        TempBehData = TempBehData(trialSubset,:);
        TempBeepData = TempBeepData (trialSubset,:);
        TempStimData = TempStimData (trialSubset,:);
        TempActionData = TempActionData (trialSubset,:);
        TempRewardData = TempRewardData (trialSubset,:);
        
    end
      
    TempStimz      = unique(TempBehData(:,2));
    StimzAbs       = unique(abs(TempBehData(:,2)));
    
    sessionz = 1;
    
elseif concatenate == 'n'
    
    sessionz = 1:length(BehPhotoM(animal_ID).Session);
    
end

% ------------------------------------------------------------------------
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
    
    TempBehData(:,7)= RT;  % put these RTs in the matrix
    [i sortingIndex] = sort(RT);
    
    largeRew = sort([(intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==-1), find(TempBehData(:,8)==1))); (intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==1), find(TempBehData(:,8)==2)))]);
    
    smallRew = sort([(intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==1), find(TempBehData(:,8)==1))); (intersect(find(TempBehData(:,9)==1 & TempBehData(:,3)==-1), find(TempBehData(:,8)==2)))]);
    
    TempBehData (:,16) = nan;
    TempBehData (largeRew,16) = 2;
    TempBehData (smallRew,16) = 1; % adding reward size to the data
    TempBehData (isnan(TempBehData(:,16)),16)= 0;
    
    correct = find(TempBehData(:,9)==1);
    error = find(TempBehData(:,9)==0);
      
    %------------------------------- plot psychometric curve------------------
    %
    
    f = figure('Position', [300 200 800 900]);
    set(gcf, 'Renderer', 'painters'); hold on
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
   
    
    xticks([min(TempStimz) 0 max(TempStimz)])
    xlabel('Contrast')
    xticks
    yticks([-1 0 1])
    yticklabels([0 0.5 1])
    ylabel('% Right', 'FontWeight', 'bold')
    set(gca, 'TickDir', 'out')
    
    % -------------- plot title ------------------------------------
    
    text(TempStimz(end), 1.2, ['ALK0' num2str(animal_name) '  Session ' num2str(iSession)], 'FontWeight', 'bold', 'FontSize', 10);
     
    % ---------- raster for all trials in session in order of RT ---------
    %     figure; hold on;
    subplot(2, 4, 5);
    
    % sort based on correct/error, contrast and then RT
    TempBehData2sort = TempBehData;
    TempBehData2sort(:,2)=abs(TempBehData2sort(:,2)); % sort based on abs of stimulus
    TempBehData2sort(TempBehData2sort(:,2)==0,2) =0.01; % to separate zero contraat rewarded and unrewarded
    
    TempBehData2sort(error,2)=-TempBehData2sort(error,2); % label error trials with negative so that they appear first
    
    [BehDatasorted j]= sortrows(TempBehData2sort,[2,16,7]); % final sorting (Correct/Error, abs contrast and RTs)
    
    if strcmp(sortDesign ,'CrrErrContRT')
        
        %%
        
        TempStimDataforR = TempStimData(j,:);
        tempj = j;
        tempRT = RT(j,:);
        StimOutcomeInt = TempBehData(:,14)-TempBehData(:,13);
        StimOutcomeInt = StimOutcomeInt(j,:);
        
        tempStimOutcomeInt = StimOutcomeInt(j,:);
        
        indexConditionChange = (find(diff(BehDatasorted(:,2)))'); % after this trial, the condition (stim, etc) changes
        
        if concatenate == 'y'
            row2add = 10;
        elseif concatenate == 'n'
            row2add = 3;
            
        end
         
      % put a band of nan between different stimulus contrasts (and error/
      % correct) but not reward sizes
        for iIndex = fliplr(indexConditionChange)
            
            TempStimDataforR = insertrows(TempStimDataforR, nan(row2add,size(TempStimData,2)), iIndex);
            tempj = insertrows(tempj, nan(row2add,1), iIndex);
            tempRT = insertrows(tempRT, nan(row2add,1), iIndex);
            tempStimOutcomeInt = insertrows(tempStimOutcomeInt, nan(row2add,1), iIndex);
            
            BehDatasorted = insertrows(BehDatasorted, nan(row2add,size(BehDatasorted,2)), iIndex);
        
        end
         
        %then plot all these rows
        imagesc((TempStimDataforR(:, 1:7400)),colorRange)
        hold on;
        
        trace = 1;
        
        for c = 1:length(TempStimDataforR)
            
            for ievent=tempRT'
                
                if trace <= size(TempStimDataforR,1) & isnan(TempStimDataforR(trace,1))==0
                    
                    H=line([downSample*ievent+eventOnset,downSample*ievent+eventOnset+20], [trace, trace]);
                    set(H,'color',[0 0 0],'LineWidth',3)
                    
                end
                
                trace  = trace + 1;
                
            end
            
        end
        
         trace = 1;
        
        
        for ievent=tempStimOutcomeInt'
            H=line([downSample*ievent+eventOnset,downSample*ievent+eventOnset+20], [trace, trace]);
            
            if BehDatasorted(trace,16) ==2
                set(H,'color',[0 0.8 0],'LineWidth',3)
            elseif BehDatasorted(trace,16) ==1
                set(H,'color',[0 1 0],'LineWidth',3)
            else
                set(H,'color',[1 0.2 0.8],'LineWidth',3)
            end
            
            trace  = trace + 1;
        end
        
        line([eventOnset eventOnset], [length(tempRT) 1], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        text(100, -20, 'Stim', 'FontWeight', 'bold', 'FontSize', 8); %stim line label
        
        %create thick green line every time
          
    elseif     strcmp( sortDesign , 'RT') % sort only based on RTs
        
        
        imagesc((TempStimData(sortingIndex, 1:7400)),colorRange)
        
        trace = 1;
        
        for ievent=RT(sortingIndex)'
            
            H=line([downSample*ievent+eventOnset,downSample*ievent+eventOnset+20], [trace, trace]);
            set(H,'color',[0 0 0],'LineWidth',3)
            
            trace  = trace + 1;
        end
     end
    
    
    colormap('bluewhitered')
    
    
    ylabel('Trials', 'FontWeight', 'bold')
    xlim([3500 7400])
    xticks([3700 4900 6100,7400])
    xticklabels({'0','1','2','3'})
    
    set(gca, 'TickDir', 'out')
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
        set(gca, 'TickDir', 'out')
        
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
        set(gca, 'TickDir', 'out')
        
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
        set(gca, 'TickDir', 'out')
        
        % ------------- 3.2 raster for action, 0.25 contrast (L+R)
        
        subplot(5, 8, 20);
        
        imagesc(TempActionData(tempIndex, (eventOnset-abs(aStart*downSample)):(eventOnset+abs(aStop*downSample))),colorRange)
        colormap('bluewhitered')
        
        line([abs(aStart*downSample) abs(aStart*downSample)], [0.5 0.5+length(TempBehData(:,2)==StimzAbs(2))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
        
        xticklabels([ ])
        yticklabels([ ])
        
        % ------------- 3.3 raster for reward, 0.25 contrast (L+R)
        
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
        set(gca, 'TickDir', 'out')
        
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
    
    %     if length(StimzAbs)==4
    %         l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),num2str(StimzAbs(4)),'location','best')
    %     elseif length(StimzAbs)==3
    %         l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),'location','best')
    %     elseif length(StimzAbs)==2
    %         l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), 'location', 'best')
    %     end
    
    xlim([0 abs(sStop*downSample)+abs(sStart*downSample)])
    xticks([1 abs(sStop*downSample)+abs(sStart*downSample)])
    xticklabels([sStart sStop])
    xlabel('Time (s)', 'FontWeight', 'bold')
    ylabel('{\Delta} F / F')
    set(gca, 'TickDir', 'out')
    
    
    % ------------- 5. line plot; avg action response over time by contrast
    
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
    
    %     if length(StimzAbs)==4
    %         l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),num2str(StimzAbs(4)),'location','best')
    %     elseif length(StimzAbs)==3
    %         l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),'location','best')
    %     elseif length(StimzAbs)==2
    %         l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), 'location', 'best')
    %     end
    
    xlim([0 abs(aStop*downSample)+abs(aStart*downSample)])
    xticks([1 abs(aStop*downSample)+abs(aStart*downSample)])
    xticklabels([aStart aStop])
    xlabel('Time (s)', 'FontWeight', 'bold')
    set(gca, 'TickDir', 'out')
    
    
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
    set(gca, 'TickDir', 'out')
     
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
    set(gca, 'TickDir', 'out')
    
    % ------------- 2. raster for stim reponse, small reward trials
    
    subplot(4, 8, 15);
    
    imagesc(TempStimData(smallRew, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample)), colorRange)
    colormap('bluewhitered')
    
    line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
    xticklabels([ ])
    ylabel('Small Reward Trials', 'FontWeight', 'bold')
    set(gca, 'TickDir', 'out')
    
    
    % ------------- 3. raster for stim reponse, NO reward trials
    
    subplot(4, 8, 23);
    
    imagesc(TempStimData((TempBehData(:,9)==0), (eventOnset-abs(sStart*downSample)):(eventOnset+abs(sStop*downSample))),colorRange)
    
    colormap('bluewhitered')
    
    line([abs(sStart*downSample) abs(sStart*downSample)], [0.5 0.5+length(TempStimData((TempBehData(:,9)==1)))], 'Color', 'black', 'LineWidth', 1.5); % stim onset line
    
    xticklabels([ ])
    ylabel('No Reward Trials', 'FontWeight', 'bold')
    set(gca, 'TickDir', 'out')
    
     
    % -------------- 4. summary avg reward response large/small/error
    
    subplot(5, 8, 39);
    
    plot(nanmean(TempStimData(largeRew, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample))),'color',colorGreen(3,:),'LineWidth',2)
    hold on;
    plot(nanmean(TempStimData(smallRew, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample))),'color',colorGreen(1,:),'LineWidth',2)
    hold on;
    plot(nanmean(TempStimData(TempBehData(:,9)==0, eventOnset-abs(sStart*downSample):eventOnset+abs(sStop*downSample))),'color',colorRed(2,:),'LineWidth',2)
    
    %     legend('Large reward','Small reward', 'No reward','location','best')
    
    xlim([0 abs(sStop*downSample)+abs(sStart*downSample)])
    xticks([1 abs(sStop*downSample)+abs(sStart*downSample)])
    xticklabels([sStart sStop])
    ylabel('{\Delta} F / F')
    xlabel('Time (s)', 'FontWeight', 'bold')
    set(gca, 'TickDir', 'out')
    
    
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
    set(gca, 'TickDir', 'out')
    
    
    % this is something temporary: 
    %%

%     figure
%     
%     for ireward = [0 1 2]
%      indexes=mintersect(find(BehDatasorted(:,16)==ireward), find(ismember(TempBehData(:,2),[-0.25 0.25])));
%        
%      hold on 
%      plot(nanmean(TempStimDataforR(indexes,:)))
%     end
%      title('stimulus aligned for 0.25 contrast')
%      xlim([3500 7400])

    figure
    
    color = 3
    for ireward = [0 1 2]
     indexes=mintersect(find(BehDatasorted(:,16)==ireward), find(ismember(TempBehData(:,2),[-0.25 0.25])));
       
     hold on 
     plot(nanmean(TempStimDataforR(indexes,:)), 'color', rewardColor(color,:),'LineWidth',2)
     color = color-1;
    end
    
     legend('No reward','Small reward', 'Large reward','location','best')
     title('stimulus aligned for 0.25 contrast')
     xlim([3500 7400])
     xticks([3700 4900 6100,7400])
     xticklabels({'0','1','2','3'})
     ylabel('Norm {\Delta} F / F', 'FontWeight', 'bold')

    
     set(gca, 'TickDir', 'out')
     xlabel('Time (s)', 'FontWeight', 'bold')
  
    
    % then call colormap(bluewhitered)
end


