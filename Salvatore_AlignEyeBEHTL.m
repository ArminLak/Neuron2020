

function Salvatore_AlignEyeBEHTL(Tag,animal_name,exp_date,exp_series)
% 
% %example
%  Salvatore_AlignEyeBEHTL('photoM','ALK071','2018-03-20','1')

% Armin Lak April 2018 London (actually in the train)

if strcmp(Tag,'photoM')
    
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];
    
else
    disp('for now the code is only working for photoM data')
    return
end

% load data and micePhotoM database
load([exp_date,'_',exp_series,'_',animal_name,'_eye_proc.mat']);
load([exp_date,'_',exp_series,'_',animal_name,'_eye.mat']);
load([exp_date,'_',exp_series,'_',animal_name,'_Timeline']);
load('MiceExpInfoPhotoM.mat')


% Align the timing of eyeLog and TimeLine and compute pupilArea according to TL timing 
SR = 1000;   % sampling rate of timeline

[animal_ID, ~] =Salvatore_Get_chan_order(animal_name);

% camera strobe on timeline
[flipTimesStrobe, flipsUpStrobe, flipsDownStrobe]=schmittTimes(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,4) , [1 4]);

% find the time of the first and last CameraStrobe in the TL
startPointEyeCamera = flipTimesStrobe(1) * SR; % first strobe on TL
endPointEyeCamera = flipTimesStrobe(end) * SR;

% timing of eyeLog data and pupil area
ft = [eyeLog.TriggerData.Time];
pupilArea = proc.data.pupil.area;

% check to see if the timing of eyeLog and TL are comparable
if abs((ft(end)-ft(1)) - (endPointEyeCamera-startPointEyeCamera)/SR) > 1
  disp('CANNOT HANDLE THIS')
  return
end

% iterpolate pupilArea data. This gives pupil data aligned on the Timelines
% timestamps
pupilArea = interp1(ft-ft(1), pupilArea,0:1/30:(endPointEyeCamera-startPointEyeCamera)/SR);

assert(sum(isnan(pupilArea)) < 10, 'interpolation failure') % a fewNaNs in the very end are acceptable
pupilArea(isnan(pupilArea)) = 0;



% ---------------------------------------
% Align photoM or ephy synch data with Timeline

% water data from TL
[flipTimesWater, flipsUpWater, flipsDownWater]=schmittTimes(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,8) , [1 4]);


% water data from ephys or photoM

if strcmp(Tag,'photoM')
    
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];
    

addpath (genpath(path2Beh))
addpath (genpath(path2photoM))


% find neuronal recording info about the session
TargetSessionFound = 0;
isession = 1;
while  TargetSessionFound == 0
  
TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[exp_date,'_',exp_series,'_',animal_name,'_Block.mat']);
    
    isession = isession + 1;
end

TargetSession = isession - 1;

photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile;

load([path2photoM,'\',photoMFileName(1:end-4),'.mat']);

WaterDig = photoMdata.TTL_1;
TimeStamps=photoMdata.Time_s_;


WaterUp = find(diff(WaterDig)>0.5);

WaterUpTime = TimeStamps(WaterUp); % watertime data from photometry data file

end


s=regstats(WaterUpTime,flipsUpWater);               
lag = s.beta(1)                 % this is the duration that photoM or ephys data was earlier than TL
% ---------------------------------------------

%ShorterTrace = min(length(flipsUpStrobe),length(pupilArea));

TimeIntervalOnTLwithMeasuredPupil =startPointEyeCamera/SR:1/30: endPointEyeCamera/SR;

%pupilArea = pupilArea(1:ShorterTrace);
eyeTimeAlingedtoNeurData = TimeIntervalOnTLwithMeasuredPupil + lag ;  % this is the time aligned to neuronal responses
%eyeTimeAlingedtoNeurData=eyeTimeAlingedtoNeurData(1:ShorterTrace);

% save
cd(['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series])


proc.data.pupil.pupilArea=pupilArea;
proc.data.pupil.AlignedTime=eyeTimeAlingedtoNeurData';


save([exp_date,'_',exp_series,'_',animal_name,'_eye_proc.mat'],'proc')


figure

plot(eyeTimeAlingedtoNeurData,pupilArea)
xlabel('alignedTime(s)')
title('PupilSize')
