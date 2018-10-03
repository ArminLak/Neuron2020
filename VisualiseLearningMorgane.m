% Visualise avg responses to stimulus and reward in specified sessions

% Morgane August 2018


clear all
close all

animal_name = 'ALK068'

z = 'y' % option to z-score data 


if animal_name == 'ALK068'
    SessionList = [1, 2, 3, 4, 5, 6];
    ylimrwd = [-2 2];
    ylimstim = [-1 2.5];
    load('BehPhotoM_Exp23')                                   % load beh data databseB
    
elseif animal_name == 'ALK070'
    SessionList = [1, 2, 3, 4, 5, 6];
    ylimrwd = [-2 3];
    ylimstim = [-1 1.5];
    load('BehPhotoM_Exp23')                                   % load beh data databse
    
elseif animal_name == 'ALK071'
    SessionList = [1, 2, 3, 4, 5, 6];
    ylimrwd = [-8 15];
    ylimstim = [-3 6];
    load('BehPhotoM_Exp23')                                   % load beh data databse

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



%--------------- useful information --------------------------------------
% task event
% 10: action time
% 12: beep 
% 13: stimulus 
% 14: reward
% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = 0 % s this should be -1 or less
stop = 2    % s

 

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses


%define time zones for normalisation after alignment

preAlign = ((sample_rate/downsampleScale)*abs(start)-800):((sample_rate/downsampleScale)*abs(start)-70);
postAlign = ((sample_rate/downsampleScale)*abs(start)+200):((sample_rate/downsampleScale)*(abs(start))+700);


% plot colours
colorGray4 = [0.8 0.8 0.8
    0.6 0.6 0.6
    0.3 0.3 0.3
    0 0 0];
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


% ---- start figure ----
figure; 
plotnum = 1;


% loop through each session in the list

for iSession = SessionList
    
    TargetSession = iSession
 
    % load Beh data and photometry data
    TrialTimingData = BehPhotoM(animal_ID).Session(TargetSession).TrialTimingData;
    
    NeuronStim = BehPhotoM(animal_ID).Session(TargetSession).NeuronStim;
    NeuronStim = NeuronStim(:, (3+start)*sample_rate/downsampleScale:(3+stop)*sample_rate/downsampleScale);
    
    NeuronReward = BehPhotoM(animal_ID).Session(TargetSession).NeuronReward;
    NeuronReward = NeuronReward(:, (3+start)*sample_rate/downsampleScale:(3+stop)*sample_rate/downsampleScale);
        
    StimzAbs=unique(abs(TrialTimingData(:,2)))';

    % ------------------- stimulus plot ------------------------------

    subplot(length(SessionList), 2, plotnum ); hold on
    
    
    c=1;
    for istim = StimzAbs
        if z == 'y'
            h = plot(zscore(nanmean(NeuronStim(abs(TrialTimingData(:,2))==istim,:))),'color',colorGray4(c,:),'LineWidth',2)
        elseif z == 'n'
            h = plot(nanmean(NeuronStim(abs(TrialTimingData(:,2))==istim,:)),'color',colorGray4(c,:),'LineWidth',2)
        end
        
        
        if length(StimzAbs)==4
            set(h, 'color', colorGray4(c,:));
        elseif length(StimzAbs)==3
            set(h, 'color', colorGray3(c,:));
        elseif length(StimzAbs)==2
            set(h, 'color', colorGray2(c,:));
        elseif length(StimzAbs)==1
            set(h, 'color', colorGray2(2,:));
        end

        c=c+1;
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
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1,:)), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1,:)), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1,:))), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1,:))), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        end
        
   elseif numel(StimzAbs) == 2
        correcteasy = abs(TrialTimingData(:,2)==[1]);
        correcthard = abs(TrialTimingData(:,2)==[0.5, 0.25, 0.12, 0]);
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1,:)), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1,:)), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1,:))), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim 
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1,:))), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim 
        end
        
   elseif numel(StimzAbs) == 1
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1,:)), 'color', colorGreen(3,:), 'LineWidth', 2)
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1,:))), 'color', colorGreen(3,:), 'LineWidth', 2)
        end
   end
    
    hold on; 
    plot(nanmean(NeuronReward(TrialTimingData(:,9)==0,:)), 'red', 'LineWidth', 2)

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


