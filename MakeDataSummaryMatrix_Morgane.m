clear all
close all

% specify anial name, experiment ID, and list of sessions of interest

% it calcuate trial-by-trial data from each session and stores them in a
% large strucutre (i.e. BehPhotoM_Exp23.mat)

% be careful because it SAVES without you having to. 

% Armin Feb 2018
% Armin July 2018 added the possbiltiy of saving 2 channels per recoring ( L and R hem)
% Morgane December 2018 modified so separates into left/right and vta/dms/nac 
        % BUT requires new system of Chan2 / Chan4 in database


animal_name = 'ALK084'
Implant = 'Bi' %Unilatral or bilateral ('Un' or 'Bi')

exp_ID = '7';


% ---------- 

if strcmp(Implant,'Un')
    ChanNum =1;
elseif strcmp(Implant,'Bi')
    ChanNum =[1 2];
end;

% ------ give a list of sessions (don't modify these numbers)

% Exp 7 session lists
SessionList = [1,3,5:13]; % ALk068 
SessionList = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12]; % ALk070
SessionList = [1:9]; %ALK071
SessionList = [1:20]; %ALK074
SessionList = [1:14]; %ALK075
SessionList = [1:3]; %ALK083
SessionList = [1:9]; %ALK084
SessionList = [1:12]; %MMM001
SessionList = [1:10]; %MMM002
SessionList = [1:9]; %MMM005


% Exp 23
SessionList = [14, 15, 16, 17, 18, 19, 20, 22, 23, 24];           % ALK068
SessionList = [13, 14, 15, 16,17, 18, 19, 20, 21, 22, 23, 24];    % ALK070
SessionList = [10, 11, 12,13, 14,19];                             % ALK071 : 15-18 are bad
SessionList = [21,22,23,24,25,26,27]; % ALK074
SessionList = [15, 16,17,18,19];      % ALK075
SessionList = [23:32];  % ALK078
SessionList = [13:21];  % ALK083
SessionList = [11:29];  % ALK084
SessionList = [13,14,15,16, 18, 19, 20, 21, 22, 23, 24, 25]; % MMM001
SessionList = [15,16,17,18,20,21,22,23]; % MMM002


%%

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

if Implant == 'Un'
    animal_chanz = cellstr([MiceExpInfo.mice(animal_ID).session(SessionList(1)).Chan2]);
elseif Implant == 'Bi'
    animal_chanz = [MiceExpInfo.mice(animal_ID).session(SessionList(1)).Chan2];
    animal_chanz = cellstr([animal_chanz; MiceExpInfo.mice(animal_ID).session(SessionList(1)).Chan4]);

    SessionList1 = find({MiceExpInfo.mice(animal_ID).session(SessionList).Chan2}==string(animal_chanz(1))); % sessions where animal_chanz(1) is chan2
    SessionList2 = find({MiceExpInfo.mice(animal_ID).session(SessionList).Chan2}==string(animal_chanz(2)));
end



for iChan = ChanNum

    load(['BehPhotoM_Exp', exp_ID, '_', char(animal_chanz(iChan))]);

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
    
        
        if iChan ==1 && ismember(iSession, SessionList1) || iChan ==2 && ismember(iSession, SessionList2)
            DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
        elseif iChan ==1 && ismember(iSession, SessionList2) || iChan ==2 && ismember (iSession, SessionList1)
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

        
        BehPhotoM(animal_ID).Session(SessionC).TrialTimingData               =  TrialTimingData;
            
            BehPhotoM(animal_ID).Session(SessionC).NeuronBeep    = Raster_MatrixBeep;
            BehPhotoM(animal_ID).Session(SessionC).NeuronStim    = Raster_MatrixStim;
            BehPhotoM(animal_ID).Session(SessionC).NeuronAction  = Raster_MatrixAction;
            BehPhotoM(animal_ID).Session(SessionC).NeuronReward  = Raster_MatrixReward;
            

    SessionC = SessionC + 1;
end

save(['BehPhotoM_Exp', exp_ID, '_', char(animal_chanz(iChan))], 'BehPhotoM')

end



