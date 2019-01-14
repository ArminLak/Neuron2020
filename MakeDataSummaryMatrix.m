
% specify anial name, experiment ID, and list of sessions of interest

% it calcuate trial-by-trial data from each session and stores them in a
% large strucutre (i.e. BehPhotoM_Exp23.mat)

% be careful because it SAVES without you having to. 

% Armin Feb 2018
% Armin July 2018 added the possbiltiy of saving 2 channels per recoring ( L and R hem)
% Morgane December 2018 modified so separates into L/R and vta/dms/nac 
        % BUT requires new system of Chan2 / Chan4 in database
        
clear all
close all


animal_name = 'ALK068'
Implant = 'Un' %Unilatral or bilateral ('Un' or 'Bi')

exp_ID = '7';


% ---------- 
if strcmp(Implant,'Un')
    ChanNum =1;
elseif strcmp(Implant,'Bi')
    ChanNum =[1 2];
end;

path2data = ['\\zubjects.cortexlab.net\Subjects\',animal_name];
addpath(genpath(path2data))

[SessionList] = getSessionList_photoM(animal_name, exp_ID);


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
%addpath('C:\Users\morga\Dropbox\Morgane Project\Code')
load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

% get brain regions and hems 
if Implant == 'Un'
    animal_chanz = cellstr([MiceExpInfo.mice(animal_ID).session(SessionList(1)).Chan2]);
    
    r1 = char(animal_chanz{1}); r1 = r1(end-2:end); % brain Region
    h1 = char(animal_chanz{1}); h1 = h1(1);         % brain Hemi
    
    SessionList1 = SessionList(find({MiceExpInfo.mice(animal_ID).session(SessionList).Chan2}==string(animal_chanz(1))));

elseif Implant == 'Bi'
    animal_chanz = [MiceExpInfo.mice(animal_ID).session(SessionList(1)).Chan2];
    animal_chanz = cellstr([animal_chanz; MiceExpInfo.mice(animal_ID).session(SessionList(1)).Chan4]);
    
    r1 = char(animal_chanz{1}); r1 = r1(end-2:end); r2 = char(animal_chanz{2}); r2 = r2(end-2:end);
    h1 = char(animal_chanz{1}); h1 = h1(1);         h2 = char(animal_chanz{2});         h2 = h2(1);

        SessionList_a = SessionList(find({MiceExpInfo.mice(animal_ID).session(SessionList).Chan2}==string(animal_chanz(1))));
        SessionList_b = SessionList(find({MiceExpInfo.mice(animal_ID).session(SessionList).Chan4}==string(animal_chanz(2))));
    SessionList1 = unique([SessionList_a; SessionList_b]); % sessions where animal_chanz(1) is chan2
    SessionList2 = setdiff(SessionList, SessionList1); % sessions where animal_chanz(1) is chan4

end


%%

for iChan = ChanNum
    
	if iChan == 1  
        hem = h1; 
        load(['BehPhotoM_Exp', exp_ID, '_', r1]);
        
    elseif iChan ==2
        hem = h2;
        
    	if r1 ~= r2 % if the second channel is looking into a new brain region
            load(['BehPhotoM_Exp', exp_ID, '_', r2]);
        end
    end
    
       
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
        if isempty(MiceExpInfo.mice(animal_ID).session(iSession).Chan2) 
            continue                                                         %skip this iteration if channel 2 was not recorded from in this session
        end
        
        DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
    elseif iChan ==1 && ismember(iSession, SessionList2) || iChan ==2 && ismember (iSession, SessionList1)
        if isempty(MiceExpInfo.mice(animal_ID).session(iSession).Chan4)
            continue                                                         %skip this iteration if channel 4 was not recorded from in this session
        end
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

        if iChan == 1 || string(r1) ~= string(r2)
            BehPhotoM(animal_ID).Session(SessionC).TrialTimingData   =  TrialTimingData;
        end
        
        if hem == 'L'
            
            BehPhotoM(animal_ID).Session(SessionC).NeuronBeepL    = Raster_MatrixBeep;
            BehPhotoM(animal_ID).Session(SessionC).NeuronStimL   = Raster_MatrixStim;
            BehPhotoM(animal_ID).Session(SessionC).NeuronActionL  = Raster_MatrixAction;
            BehPhotoM(animal_ID).Session(SessionC).NeuronRewardL  = Raster_MatrixReward;
            
        elseif hem == 'R'
            BehPhotoM(animal_ID).Session(SessionC).NeuronBeepR    = Raster_MatrixBeep;
            BehPhotoM(animal_ID).Session(SessionC).NeuronStimR   = Raster_MatrixStim;
            BehPhotoM(animal_ID).Session(SessionC).NeuronActionR  = Raster_MatrixAction;
            BehPhotoM(animal_ID).Session(SessionC).NeuronRewardR  = Raster_MatrixReward;
            
        end
        
    SessionC = SessionC + 1;
end

if strcmpi(getComputerName,'zopamine2')
    cd ('C:\Users\Armin\Dropbox\Work\UCL\Science\Analysis Code\PhotoM')
elseif strcmpi(getComputerName, 'zebrafish')
    cd ('C:\Users\morga\Documents\MATLAB')
end

        if length(ChanNum) ==1 || (length(ChanNum)==2 && iChan ==1 && string(r1) ~= string(r2)) 
            save(['BehPhotoM_Exp', exp_ID, '_', r1], 'BehPhotoM', '-v7.3');
        elseif iChan ==2
            save(['BehPhotoM_Exp', exp_ID, '_', r2], 'BehPhotoM', '-v7.3');
        end
        
    end
    


