clear all

% for now only works for ALK084. 
% must have done eye processing and alignment

animal_name = 'ALK084'
exp_date = '2018-11-14'
exp_series = '2'

start_time_s = 1; %interval of video
stop_time_s = 10;


sample_rate = 12000;                                        % photoM recording sampling rate
video_fps = 30;                                             % video frames per second 
downsampleScale = sample_rate/video_fps;                    % factor downsampling the Ca responses to match video frame rate

load('MiceExpInfoPhotoM')

% ------ read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);


% ------ get paths to photoM, beh, and video data
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];
% path2eye = ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];
filename = [exp_date,'_',exp_series,'_',animal_name,'_eye.mj2'];

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
Video = VideoReader([exp_date,'_',exp_series,'_',animal_name,'_eye.mj2']);

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

DeltaFoverF = DeltaFoverF(round(sample_rate*VideoFrameTimes(1)):round(sample_rate*VideoFrameTimes(end))); %align video to photom
TimeStamps = TimeStamps(round(sample_rate*VideoFrameTimes(1)):round(sample_rate*VideoFrameTimes(end)));

% ----- create sound and frame stack

startframe = round((VideoFrameTimes(1)+start_time_s)*video_fps); %start frame
stopframe = round((VideoFrameTimes(1)+stop_time_s)*video_fps); %stop frame 

% DeltaFoverF = downsample(DeltaFoverF, sample_rate/video_fps);
%%
duration = range(VideoFrameTimes);
t = 0:1/sample_rate:duration;
Fc = 200;               % baseline frequency
FDev = 200;             % frequency deviation 
Fs = sample_rate/video_fps;    % sampling frequency set so that audio samples per second matches video frames per second
soundwave = fmmod(smooth(DeltaFoverF,10), Fc, Fs, FDev);
% sound(soundwave, 1000)
[frame_stack] = getframes(path2Beh, filename, startframe, stopframe); % create framestack

% ------ write video file

% nRep = floor(length(sound)/length(frame_stack)); % audio frames per video frame
% nDiff = length(sound) - nRep*length(frame_stack); %offset between audio and video 
% 
% if nDiff
%         % if length(frame_stack) does not evenly divide nsoundwave, then subsample audio to match nRep*length(frame_stack)
%         selector = round(linspace(1, length(sound), nRep*length(frames)));
%         subsoundwave = sound(selector, :);
% end

VideoFWriter = vision.VideoFileWriter(fullfile(path2Beh, [exp_date,'_',exp_series,'_',animal_name,'_sonifiedDA.avi']), ...
    'FileFormat', 'Indexed AVI', 'FrameRate', video_fps, 'AudioInputPort', true);
VideoFWriter.Colormap = parula(256)

FS_dim = size(frame_stack);
audio_frame_length = length(soundwave)/FS_dim(3)
audio_c = 0

for iFrame = 1:FS_dim(3)
    videoFrame = frame_stack(:,:,iFrame);
    soundFrame = soundwave(audio_c+1:(audio_c+audio_frame_length));
    
    step(VideoFWriter, videoFrame, soundFrame);
    
	audio_c = audio_c + audio_frame_length;
end



% to do:
% fix Fs for creating soundwave. be sure total sampling frequency makes sense 


function [frame_stack] = getframes(path, filename, startframe, stopframe)

    frame_stack = [];

    addpath(path);
    vidobj = VideoReader(fullfile(path,filename));
    
    framec = 1;
    for iFrame = startframe:stopframe
        frame_stack(:,:,framec)=read(vidobj,iFrame);
        framec = framec+1;
    end
    
    
end
