clear all
close all

% specify anial name and list of sessions of interest
% it calcuate trial-by-trial data from each session and stores them in a
% large strucutre (i.e. BehPhotoM_Exp23.mat)

% if you are happy with it, then save it at our Github folder.


animal_name = 'ALK068'

% give a list of sessions

SessionList = [14, 15, 16, 17, 18, 19, 20, 22, 23, 24]; % ALK068 Exp23
%SessionList = [13, 14, 15, 16,17, 18, 19, 20, 21, 22, 23, 24]; % ALK070 Exp23
%SessionList = [10, 11, 12,13, 14]; % ALK071 Exp23


% This is the structure that will hold the data
load('BehPhotoM_Exp23.mat')

%--------------- useful information --------------------------------------
% task event
% 10: action time
% 12: beep
% 13: stimulus
% 14: reward
% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = -3 % s
stop=3     % s

load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses


% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

SessionC = 1;
for iSession =  SessionList

  
%--------------------------------------------------------------------------
% load Beh data and photometry data
TrialTimingData = MiceExpInfo.mice(animal_ID).session(iSession).TrialTimingData;
%TrialTimingData(TrialTimingData(:,3)==-1,3)=0; % define left choice as 0 (rather than -1)

ReactionTime = TrialTimingData(:,10) -TrialTimingData(:,13);

Stimz=unique(TrialTimingData(:,2))';
StimzAbs=unique(abs(TrialTimingData(:,2)))';

% delay for wheel movement
FileAlignDelay = MiceExpInfo.mice(animal_ID).session(iSession).AlignDelay;

% load photoM data
photoMFileName=MiceExpInfo.mice(animal_ID).session(iSession).Neuronfile(1:end-4);
%photoMdata = readtable([path2photoM,'\',photoMFileName]);
load(photoMFileName);
DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
TimeStamps=photoMdata.Time_s_;


%------------------------define event time for event-alinged responses--------------------------

event_times = TrialTimingData(:,12); % initial beep onset

[Raster_MatrixBeep]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);

event_times = TrialTimingData(:,13); % vis stimulus onset

[Raster_MatrixStim]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);

event_times = TrialTimingData(:,10); % action onset

event_times(find(isnan(TrialTimingData(:,10))))= TrialTimingData(find(isnan(TrialTimingData(:,10))),13)+0.2;

[Raster_MatrixAction]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);

event_times = TrialTimingData(:,14); %reward onset

[Raster_MatrixReward]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);

BehPhotoM(animal_ID).Session(SessionC).TrialTimingData =  TrialTimingData;
BehPhotoM(animal_ID).Session(SessionC).NeuronBeep  = Raster_MatrixBeep;
BehPhotoM(animal_ID).Session(SessionC).NeuronStim  = Raster_MatrixStim;
BehPhotoM(animal_ID).Session(SessionC).NeuronAction  = Raster_MatrixAction;
BehPhotoM(animal_ID).Session(SessionC).NeuronReward  = Raster_MatrixReward;

 SessionC = SessionC + 1;
end
 



