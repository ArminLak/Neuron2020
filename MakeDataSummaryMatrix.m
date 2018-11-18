clear all
close all

% specify anial name and list of sessions of interest

% it calcuate trial-by-trial data from each session and stores them in a
% large strucutre (i.e. BehPhotoM_Exp23.mat)

% if you are happy with it, then save it

% Armin Feb 2018
% Armin July 2018 added the possbiltiy of saving 2 channels per recoring ( L and R
% hem)


animal_name = 'ALK083'

%Unilatral or bilateral ('Un' or 'Bi')
Implant = 'Un'

if strcmp(Implant,'Un')
    ChanNum =1;
elseif strcmp(Implant,'Bi')
    ChanNum =[1 2];
end;

% ------ give a list of sessions (don't modify these numbers)

% VTA animals, Exp 23
SessionList = [14, 15, 16, 17, 18, 19, 20, 22, 23, 24];           % ALK068 Exp23
SessionList = [13, 14, 15, 16,17, 18, 19, 20, 21, 22, 23, 24];    % ALK070 Exp23
SessionList = [10, 11, 12,13, 14,19];                             % ALK071 Exp23 % 15-18 are bad

% DMS animals Exp 23
SessionList = [21,22,23,24,25,26,27]; % ALK074, exp 23
SessionList = [15, 16,17,18,19];      % ALK075, exp 23

% NAc animals Exp 23

 SessionList = [13,14,15,16, 18, 19, 20, 21, 22, 23, 24, 25]; % MMM001, exp 23
 SessionList = [23:32];  % ALK078, 
 SessionList = [15,16,17,18,20,21,22,23]; % MMM002, exp 23, 

 %% Learning Exp
 % Naive learning VTA animals
SessionList = [1,3,5:13]; % ALk068, Exp 7
SessionList = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12]; % ALk070, Exp 7
SessionList = [1:9]; %ALK071 Exp 7
SessionList = [1:3]; %ALK083 Exp 7

%  % Naive learning DMS animals
% SessionList = [1:20]; %ALK074 Exp 7
% SessionList = [1:14]; %ALK075 Exp 7
% 
%  % Naive learning NAc animals
% SessionList = [1:12]; %MMM001 Exp 7
% SessionList = [1:10]; %MMM002 Exp 7

%%
% ------------------------------------------------------------------------

% This is the structure that will hold the data

%select database in case it is saved, u can load it to add
%more animals

%load('BehPhotoM_Exp23')

%load('BehPhotoM_Exp23_NAc')

%load('BehPhotoM_Exp23_DMS')

%load('BehPhotoM_Exp7_VTA')

load('BehPhotoM_Exp7_DMS')

% load('BehPhotoM_Exp7_NAc')

%--------------- useful information --------------------------------------
% task event
% 10: action time
% 12: beep
% 13: stimulus
% 14: reward
% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = -3 % s
stop=8     % s

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
    
    for iChan = ChanNum
        
        if iChan == 1
            DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
        elseif iChan == 2
            DeltaFoverF = photoMdata.AnalogIn_4_dF_F0;
        end
        
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
        
%         if length(ChanNum) == 1
%             HemChan = [1];
%         elseif length(ChanNum) == 2
%             if MiceExpInfo.mice(animal_ID).session(iSession).RChan == 2
%                 HemChan = [2 1];
%             elseif MiceExpInfo.mice(animal_ID).session(iSession).RChan == 4
%                 HemChan = [1 2];
%                 
%             end
%         end
        
        BehPhotoM(animal_ID).Session(SessionC).TrialTimingData               =  TrialTimingData;
        
        
        if length(ChanNum) == 1 % only one channel
            
            BehPhotoM(animal_ID).Session(SessionC).NeuronBeep    = Raster_MatrixBeep;
            BehPhotoM(animal_ID).Session(SessionC).NeuronStim    = Raster_MatrixStim;
            BehPhotoM(animal_ID).Session(SessionC).NeuronAction  = Raster_MatrixAction;
            BehPhotoM(animal_ID).Session(SessionC).NeuronReward  = Raster_MatrixReward;
            
        elseif length(ChanNum) == 2
            
            if iChan ==1 && MiceExpInfo.mice(animal_ID).session(iSession).RChan == 2 % Analog2 looks at Right Hem
                
                BehPhotoM(animal_ID).Session(SessionC).NeuronBeepR    = Raster_MatrixBeep;
                BehPhotoM(animal_ID).Session(SessionC).NeuronStimR    = Raster_MatrixStim;
                BehPhotoM(animal_ID).Session(SessionC).NeuronActionR  = Raster_MatrixAction;
                BehPhotoM(animal_ID).Session(SessionC).NeuronRewardR  = Raster_MatrixReward;
            end
            if iChan ==2 && MiceExpInfo.mice(animal_ID).session(iSession).RChan == 2 % Analog2 looks at Right Hem
                
                BehPhotoM(animal_ID).Session(SessionC).NeuronBeep    = Raster_MatrixBeep;
                BehPhotoM(animal_ID).Session(SessionC).NeuronStim    = Raster_MatrixStim;
                BehPhotoM(animal_ID).Session(SessionC).NeuronAction  = Raster_MatrixAction;
                BehPhotoM(animal_ID).Session(SessionC).NeuronReward  = Raster_MatrixReward;
                
            end
            
            if iChan ==2 && MiceExpInfo.mice(animal_ID).session(iSession).RChan == 4 % Analog2 looks at left Hem
                
                BehPhotoM(animal_ID).Session(SessionC).NeuronBeepR    = Raster_MatrixBeep;
                BehPhotoM(animal_ID).Session(SessionC).NeuronStimR    = Raster_MatrixStim;
                BehPhotoM(animal_ID).Session(SessionC).NeuronActionR  = Raster_MatrixAction;
                BehPhotoM(animal_ID).Session(SessionC).NeuronRewardR  = Raster_MatrixReward;
            end
            if iChan ==1 && MiceExpInfo.mice(animal_ID).session(iSession).RChan == 4 % Analog2 looks at left Hem
                
                BehPhotoM(animal_ID).Session(SessionC).NeuronBeep    = Raster_MatrixBeep;
                BehPhotoM(animal_ID).Session(SessionC).NeuronStim    = Raster_MatrixStim;
                BehPhotoM(animal_ID).Session(SessionC).NeuronAction  = Raster_MatrixAction;
                BehPhotoM(animal_ID).Session(SessionC).NeuronReward  = Raster_MatrixReward;
                
            end
            
        end
        

    end
    
    SessionC = SessionC + 1;
end




