clear all

% for now only works for ALK084. 
% must have done eye processing and alignment

animal_name = 'ALK084'
exp_date = '2018-11-14'
exp_series = '2'

load('MiceExpInfoPhotoM')

sample_rate = 12000;                                        % photoM recording sampling rate
video_fps = 30;                                             % video frames per second 
downsampleScale = sample_rate/video_fps;                    % factor downsampling the Ca responses to match video frame rate


% ------ read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);


% ------ get paths to photoM, beh, and video data
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];
% path2eye = ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];

addpath (genpath(path2Beh))
addpath (genpath(path2photoM))


% ------ find neuronal recording info about the session
TargetSessionFound = 0;
isession = 1;
while  TargetSessionFound == 0
    
    TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[exp_date,'_',exp_series,'_',animal_name,'_Block.mat']);
    
    isession = isession + 1;
end
TargetSession = isession - 1;


% ------ load Beh data and video data
TrialTimingData = MiceExpInfo.mice(animal_ID).session(TargetSession).TrialTimingData;

load([exp_date,'_',exp_series,'_',animal_name,'_eye_proc.mat']);
VideoFrameTimes = proc.data.pupil.AlignedTime; 


% ------ load photoM data
photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile(1:end-4);
load(photoMFileName);
TimeStamps=photoMdata.Time_s_;
if strcmp(MiceExpInfo.mice(animal_ID).session(TargetSession).Chan2(end-2:end), 'VTA')
    DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
elseif strcmp(MiceExpInfo.mice(animal_ID).session(TargetSession).Chan4(end-2:end), 'VTA')
    DeltaFoverF = photoMdata.AnalogIn_4_dF_F0;
end


% ------- cut photoM data to fit length of video (by times)


DeltaFoverF = DeltaFoverF(round(sample_rate*VideoFrameTimes(1)):round(sample_rate*VideoFrameTimes(end)));
DeltaFoverF = downsample(DeltaFoverF, downsampleScale);



% ----- create a wave 

amp = 10;
Fs = 1000;                                                  %sampling frequency must be less than min freq in audio

Audio = rescale(DeltaFoverF, 1500, 6000);






% to do: 
% align photoM data to video data 
% filter out so u can add any animal and it looks for specific brain region