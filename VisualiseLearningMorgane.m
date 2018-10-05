% Visualise avg responses to stimulus and reward in specified sessions

% Morgane August 2018


clear all
close all

animal_name = 'ALK070'

z = 'n' % option to z-score data 


if animal_name == 'ALK068'
    SessionList = [1:11];
    ylimrwd = [-2 2];
    ylimstim = [-1 2.5];
    load('BehPhotoM_Exp7_VTA')                                   % load beh data databseB
    
elseif animal_name == 'ALK070'
    SessionList = [1:9]; %weaker signal in sessions 10 and 11 
    ylimrwd = [-2 3];
    ylimstim = [-1 1.5];
    load('BehPhotoM_Exp7_VTA')                                   % load beh data databse
    
elseif animal_name == 'ALK071'
    SessionList = [1:9];
    ylimrwd = [-8 15];
    ylimstim = [-3 6];
    load('BehPhotoM_Exp7_VTA')                                   % load beh data databse

elseif animal_name == 'MMM001'
    SessionList = [1, 2, 3, 5, 7]; % to review: error with session 6
    ylimrwd = [-3 6];
    ylimstim = [-3 2];
    load('BehPhotoM_Exp23_NAc')                                   % load beh data databse

elseif animal_name == 'MMM002'
    SessionList = [1, 2, 3, 4, 5, 6, 7];
    ylimrwd = [-4 7];
    ylimstim = [-3 3];
    load('BehPhotoM_Exp23_NAc')                                   % load beh data databse

end

if z == 'y'
    ylimstim = [-3 3];
    ylimrwd = [-3 3];
end


% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = -1 % s this should be -1 or less
stop = 2    % s
event_time = 3; % this is the time in the summary matrix where the event took place

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

start2stop = (event_time+start)*sample_rate/downsampleScale:(event_time+stop)*sample_rate/downsampleScale; %window of interest


% event is at 3 seconds i.e. point 3600

preAlignStim = (event_time-0.4)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignStim = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;
preAlignRwd = (event_time-0.2)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignRwd = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;


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

%-------------------------------find path, add path and load data----------------------------------
% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

% -------- prep for summary figures --------

StimResponseSummary = nan(max(SessionList), 5); % 1 to 5 are 0, 0.12, 0.25, 0.5, 1.0
RwdResponseSummary = nan(max(SessionList), 5);

% ---- start figure ----
figure; 
plotnum = 1;


% loop through each session in the list

for iSession = SessionList
    
    TargetSession = iSession
 
    % load Beh data and photometry data
    TrialTimingData = BehPhotoM(animal_ID).Session(TargetSession).TrialTimingData;
    TrialTimingDataCor = TrialTimingData(TrialTimingData(:,9)==1, :);
    
    NeuronStim = BehPhotoM(animal_ID).Session(TargetSession).NeuronStim;
    NeuronStimCor = NeuronStim(TrialTimingData(:,9)==1, :);
    
    NeuronReward = BehPhotoM(animal_ID).Session(TargetSession).NeuronReward;
    NeuronRewardCor = NeuronReward(TrialTimingData(:,9)==1, :);
        
    %StimzAbs=unique(abs(TrialTimingData(:,2)))';
    StimzAbs = [0 0.12 0.25 0.5 1.0];

    % ------------------- stimulus plot ------------------------------

    subplot(length(SessionList), 2, plotnum ); hold on
    
    
  %  c=4;
    
    
    for stimcount = 1:length(StimzAbs)
        stimcolumn = 6 - stimcount; % % column 1 is 0 contrast and column 5 is 1.0 contrast
        istim = StimzAbs(stimcount);
        
        if z == 'y'
            h = plot(zscore(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==istim, start2stop))),'color',colorGray4(stimcount,:),'LineWidth',2)
        elseif z == 'n'
            h = plot(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==istim, start2stop)),'color',colorGray4(stimcount,:),'LineWidth',2)
        end
        
        StimResponseSummary(TargetSession,stimcolumn)  = nanmean(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), postAlignStim)))...
            - nanmean(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), preAlignStim))); %difference between before and after; indicator of relative signal change 
        RwdResponseSummary(TargetSession,stimcolumn)  = nanmean(nanmean(NeuronRewardCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), postAlignRwd))) ...
            - nanmean(nanmean(NeuronRewardCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), preAlignRwd))); %difference between before and after; indicator of relative signal change 

        
        
        if length(StimzAbs)==4
            set(h, 'color', colorGray4(stimcount,:));
        elseif length(StimzAbs)==3
            set(h, 'color', colorGray3(stimcount,:));
        elseif length(StimzAbs)==2
            set(h, 'color', colorGray2(stimcount,:));
        elseif length(StimzAbs)==1
            set(h, 'color', colorGray2(2,:));
        end

   %     c=c-1;
    end
    

    
    
    if plotnum == 1
        title('Stimulus response')
        if length(StimzAbs)==4
            l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),num2str(StimzAbs(4)),'location','northeast');
        elseif length(StimzAbs)==3
            l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),'location','northeast');
        elseif length(StimzAbs)==2
            l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), 'location', 'northeast');
        end
    end
    
    ylim(ylimstim)
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])
    ylabel('{\Delta} F / F')
    
    if plotnum == length(SessionList)*2 - 1
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
    end

    plotnum = plotnum + 1;
    
    
    
    
    
    % ------------------- reward plot --------------------------------
    
    subplot(length(SessionList), 2, plotnum ); hold on
    

    
   if numel(StimzAbs) > 2
        correcteasy = abs(TrialTimingData(:,2)==[1, 0.5]);
        correcthard = abs(TrialTimingData(:,2)==[0.25, 0.12, 0]);
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop)), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop)), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop))), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop))), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        end
        
   elseif numel(StimzAbs) == 2
        correcteasy = abs(TrialTimingData(:,2)==[1]);
        correcthard = abs(TrialTimingData(:,2)==[0.5, 0.25, 0.12, 0]);
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop)), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop)), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop))), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop))), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        end
        
   elseif numel(StimzAbs) == 1
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1, start2stop)), 'color', colorGreen(3,:), 'LineWidth', 2)
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1, start2stop))), 'color', colorGreen(3,:), 'LineWidth', 2)
        end
   end
    
    hold on; 
    plot(nanmean(NeuronReward(TrialTimingData(:,9)==0, start2stop)), 'red', 'LineWidth', 2) % plot error trials 

    ylim(ylimrwd)
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])
    
    if plotnum == 2
        title('Reward response')
        legend('Easy & correct', 'Hard & correct', 'Error', 'location', 'northeast')
    end

    if plotnum == length(SessionList)*2
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
    end
    
    plotnum = plotnum + 1;


end
%%
% ---- plot summaries ----

figure; hold on 

subplot(1, 2, 1) % stim first 
c = 1;
for i = 1:size(StimResponseSummary, 2)
    
    plot(StimResponseSummary(:, i), 'o-', 'color', colorGray4(c,:), 'markerFaceColor', colorGray4(c,:), 'MarkerSize', 3.5, 'LineWidth', 1)
    
    c = c+1;
    
    hold on
    
end
c = 1;

subplot(1, 2, 2) % reward second 

for i = 1:size(StimResponseSummary, 2)
    
    plot(RwdResponseSummary(:, i), 'o-', 'color', colorGray4(c,:), 'markerFaceColor', colorGray4(c,:), 'MarkerSize', 3.5, 'LineWidth', 1)
    
    c = c+1;
    
    hold on
    
end


