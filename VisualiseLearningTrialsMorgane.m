% Visualise individual trials in early sessions 3D 

% Morgane Sept 2018


clear all
% close all

animal_name = 'MMM001'


if animal_name == 'MMM001'
    TargetSession = 1
%     SessionList = [1, 2, 3, 5, 7]; % to review: error with session 6
%     exp_dates   = [{'2018-07-09', '2018-07-10', '2018-07-11', '2018-07-12', '2018-07-16'}];
%     exp_series =[{'1', '1', '1', '1', '3'}];
%     ylimrwd = [-3 6];
%     ylimstim = [-3 2];
% elseif animal_name == 'ALK071'
%     SessionList = [1, 3, 4, 5];
%     ylimrwd = [-1 15];
%     ylimstim = [-2 8];
% elseif animal_name == 'ALK070'
%     SessionList = [1, 2, 3, 4, 5];
%     ylimrwd = [-1 2];
%     ylimstim = [-1 2];
% elseif animal_name == 'MMM002'
%     SessionList = [1, 2, 3, 4, 5, 6, 7];
%     ylimrwd = [-4 7];
%     ylimstim = [-3 3];
% elseif animal_name == 'ALK068'
%     SessionList = [1, 3, 5, 6];
%     ylimrwd = [-4 4];
%     ylimstim = [-1 3];
end

start = 0 % s this should be -1 or less
stop = 2    % s

SmoothFactor = 5

load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

% ------------ first get the Raster_MatrixStim ------------------------

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
    
    Z = Raster_MatrixStim'; % get trial # on x axis, time on y axis
    
   
    Z = smooth2a(Z, SmoothFactor);
    Z = flipud(Z); % get y axis to read left to right


% ------------ plot ---------------------------------------------


figure; 
[X, Y] = meshgrid(1:174, 1:2300);
h = surf(X, Y, Z)
hold on;
set(h, 'LineStyle', 'none')
    ylim([0 (stop-start)* sample_rate]/downsampleScale)
    yticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    yticklabels(fliplr([start: 1 : stop]))


