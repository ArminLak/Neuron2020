% Visualise avg responses to stimulus and reward in specified sessions

% Morgane August 2018


clear all
close all

animal_name = 'ALK068'




if animal_name == 'MMM001'
    SessionList = [1, 2, 3, 5, 7]; % to review: error with session 6
%     exp_dates   = [{'2018-07-09', '2018-07-10', '2018-07-11', '2018-07-12', '2018-07-16'}];
%     exp_series =[{'1', '1', '1', '1', '3'}];
    ylimrwd = [-3 6];
    ylimstim = [-3 2];
elseif animal_name == 'ALK071'
    SessionList = [1, 3, 4, 5];
%     exp_dates   = [{'2018-02-23', '2018-02-27', '2018-02-28', '2018-03-01', '2018-03-02'}];
%     exp_series =[{'1', '1', '1', '1', '3'}];
    ylimrwd = [-1 15];
    ylimstim = [-2 8];
elseif animal_name == 'ALK070'
    SessionList = [1, 2, 3, 4, 5];
%     exp_dates   = [{'2018-01-24', '2018-01-25', '2018-01-29', '2018-01-30'}];
%     exp_series =[{'1', '3', '2', '2'}];
    ylimrwd = [-1 2];
    ylimstim = [-1 2];
elseif animal_name == 'MMM002'
    SessionList = [1, 2, 3, 4, 5, 6, 7];
%     exp_dates   = [{'2018-08-01', '2018-08-02', '2018-08-06', '2018-08-07', '2018-08-08', '2018-08-10'}];
%     exp_series =[{'1', '1', '1', '1', '1', '2'}];
    ylimrwd = [-4 7];
    ylimstim = [-3 3];
elseif animal_name == 'ALK068'
    SessionList = [1, 3, 5, 6];
    ylimrwd = [-4 4];
    ylimstim = [-1 3];
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


load('MiceExpInfoPhotoM')                                   % load beh data databse
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
    TrialTimingData = MiceExpInfo.mice(animal_ID).session(TargetSession).TrialTimingData;
    
    StimzAbs=unique(abs(TrialTimingData(:,2)))';
    
    % delay for wheel movement
    FileAlignDelay = MiceExpInfo.mice(animal_ID).session(TargetSession).AlignDelay;

    % load photoM data
    photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile(1:end-4);
    load(photoMFileName);
    
    DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
    TimeStamps=photoMdata.Time_s_;
    
    
    % ------------------- stimulus plot ------------------------------
    
    event_times = TrialTimingData(:,13); % stimulus onset
    [Raster_MatrixStim]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);

    subplot(length(SessionList), 2, plotnum ); hold on
    
    
    c=1;
    for istim = StimzAbs
        h = plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim,:)),'color',colorGray4(c,:),'LineWidth',2)

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
    
    if length(StimzAbs)==4
        l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),num2str(StimzAbs(4)),'location','best');
    elseif length(StimzAbs)==3
        l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),'location','best');
    elseif length(StimzAbs)==2
        l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), 'location', 'best');
    end
    
    
    if plotnum == 1
        title('Stimulus response')
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
    
    event_times = TrialTimingData(:,14); %reward onset
    [Raster_MatrixReward]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);
    
    subplot(length(SessionList), 2, plotnum ); hold on
    
    plot(nanmean(Raster_MatrixReward(TrialTimingData(:,9)==1,:)), 'green', 'LineWidth', 2)
    plot(nanmean(Raster_MatrixReward(TrialTimingData(:,9)==0,:)), 'red', 'LineWidth', 2)
    
    if plotnum == 2
        title('Reward response')
    end
    
    ylim(ylimrwd)
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])

    if plotnum == length(SessionList)*2
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
    end
    
    plotnum = plotnum + 1;


end


