% Visualise individual trials in early sessions 3D 

% Morgane Sept 2018


clear all
% close all

animal_name = 'MMM002'


if animal_name == 'MMM001'
    SessionList = [1, 2, 3, 4, 5, 7]; % to review: error with session 6
%     exp_dates   = [{'2018-07-09', '2018-07-10', '2018-07-11', '2018-07-12', '2018-07-16'}];
%     exp_series =[{'1', '1', '1', '1', '3'}];
elseif animal_name == 'ALK071'
    SessionList = [1, 3, 4, 5];
elseif animal_name == 'ALK070'
    SessionList = [2, 3, 4, 5];
elseif animal_name == 'MMM002'
    SessionList = [1, 2, 3, 4, 5, 6, 7];
elseif animal_name == 'ALK068'
    SessionList = [1, 3, 5, 6];
end

start = 0 % s this should be -1 or less
stop = 2    % s

SmoothFactor = 10

load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

Zstim = [];
Zrew = [];

% ------------ first get the Raster_MatrixStim ------------------------

for TargetSession = SessionList 
    
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
    
    event_times = TrialTimingData(:,13); % stimulus onset
    [Raster_MatrixStim]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);
    
    event_times = TrialTimingData(:,14); %reward onset
    [Raster_MatrixReward]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);
    
%     tempZ = tempRaster_MatrixStim'; % get trial # on x axis, time on y axis
    
    
%     tempZ = smooth2a(tempZ, SmoothFactor);
%     tempZ = flipud(tempZ); % get y axis to read left to right
    
    Zstim = [Zstim; Raster_MatrixStim]; % 
    Zrew = [Zrew; Raster_MatrixReward];

end

Zstim = flipud(Zstim'); % get trial # on x axis, time on y axis and y axis to read left to right 
Zstim = smooth2a(Zstim, SmoothFactor);

Zrew = flipud(Zrew');
Zrew = smooth2a(Zrew, SmoothFactor);



% ------------ plot stim resp---------------------------------------------


figure; 
% [X, Y] = meshgrid(1:size(Z, 2), 1:size(Z, 1));
h = surf(Zstim)
hold on;
set(h, 'LineStyle', 'none')
ylim([0 (stop-start)* sample_rate]/downsampleScale)
yticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
yticklabels(fliplr([start: 1 : stop]))
ylabel('Time (s)')
xlabel('Trial #')
zlabel('Response')
stimTitle = [animal_name,' Stim responses'];
title(stimTitle)
colorbar

% ---------- plot reward resp -----------------------------------------

figure; 

h = surf(Zrew)
hold on;
set(h, 'LineStyle', 'none')
ylim([0 (stop-start)* sample_rate]/downsampleScale)
yticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
yticklabels(fliplr([start: 1 : stop]))
ylabel('Time (s)')
xlabel('Trial #')
zlabel('Response')
rewTitle = [animal_name,' Reward responses'];
title(rewTitle)
colorbar

