clear all
close all

animal_name = 'MMM001'

if animal_name == 'MMM001'
    exp_dates   = [{'2018-07-09', '2018-07-10', '2018-07-11', '2018-07-12', '2018-07-16'}];
    exp_series =[{'1', '1', '1', '1', '3'}];
    ylimrwd = [-1.5 2];
elseif animal_name == 'ALK

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

for iSession = 1:numel(exp_dates)

    temp_date = exp_dates{iSession};
    temp_series = exp_series{iSession};

    % find path to beh and photometry data  
    path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',temp_date,'\photoM'];
    path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',temp_date,'\',temp_series];

    addpath (genpath(path2Beh))
    addpath (genpath(path2photoM))

    
    % load wheel data
    load([temp_date,'_',temp_series,'_',animal_name,'_Block.mat']);
    
    
    % find neuronal recording info about the session
    TargetSessionFound = 0;
    isession = 1;
    while  TargetSessionFound == 0

        TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[temp_date,'_',temp_series,'_',animal_name,'_Block.mat']);

        isession = isession + 1;
    end
    TargetSession = isession - 1;

    
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

    subplot(numel(exp_dates), 2, plotnum ); hold on
    
    
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
    
    
    if plotnum == 1
        title('Stimulus response')
    end
    
    ylim([-1 2])
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])
    ylabel('{\Delta} F / F')
    if plotnum == numel(exp_dates)*2 - 1
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
    end

    plotnum = plotnum + 1;
    % ------------------- reward plot --------------------------------
    
    event_times = TrialTimingData(:,14); %reward onset
    [Raster_MatrixReward]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);
    
    subplot(numel(exp_dates), 2, plotnum ); hold on
    
    plot(nanmean(Raster_MatrixReward), 'k', 'LineWidth', 2)
    
    
    if plotnum == 2
        title('Reward response')
    end
    
    ylim([-1.5 2])
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])

    if plotnum == numel(exp_dates)*2
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
    end
    
    plotnum = plotnum + 1;


end


