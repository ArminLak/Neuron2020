clear all
close all

animal_name = 'ALK071'

exp_date   = '2018-03-16'
exp_series ='1'

%--------------- useful information --------------------------------------
% task event
% 10: action time
% 12: beep
% 13: stimulus
% 14: reward
% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = -1 % s
stop=1     % s

load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

%---------------- colors for plotting-------------------------------------

colorGray = [0 0 0
    0.3 0.3 0.3
    0.6 0.6 0.6
    0.8 0.8 0.8];
colorGreen = [0 1 0
    0 0.8 0
    0 0.6 0
    0  0.3 0];
colorRed = [1 0 0
    0.8 0 0
    0.6 0  0
    0.3 0 0];
%-------------------------------find path, add path and load data----------------------------------
% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

% find path to beh and photometry data
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];

addpath (genpath(path2Beh))
addpath (genpath(path2photoM))

% load wheel data
load([exp_date,'_',exp_series,'_',animal_name,'_Block.mat']);
%--------------------------------------------------------------------------
% find neuronal recording info about the session
TargetSessionFound = 0;
isession = 1;
while  TargetSessionFound == 0
    
    TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[exp_date,'_',exp_series,'_',animal_name,'_Block.mat']);
    
    isession = isession + 1;
end

TargetSession = isession - 1;
%--------------------------------------------------------------------------
% load Beh data and photometry data
TrialTimingData = MiceExpInfo.mice(animal_ID).session(TargetSession).TrialTimingData;
TrialTimingData(TrialTimingData(:,3)==-1,3)=0; % define left choice as 0 (rather than -1)

ReactionTime = TrialTimingData(:,10) -TrialTimingData(:,13);

Stimz=unique(TrialTimingData(:,2))';
StimzAbs=unique(abs(TrialTimingData(:,2)))';

% delay for wheel movement
FileAlignDelay = MiceExpInfo.mice(animal_ID).session(TargetSession).AlignDelay;

% load photoM data
photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile(1:end-4);
%photoMdata = readtable([path2photoM,'\',photoMFileName]);
load(photoMFileName);
DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
TimeStamps=photoMdata.Time_s_;

    
%------------------------------- plot psychoimetric curve------------------
% psych curve:
c = 1;
for istim = Stimz
    
    performance(c) = nanmean (TrialTimingData(TrialTimingData(:,2)==istim,3));
    RT(c) = nanmean (ReactionTime(TrialTimingData(:,2)==istim));
    
    c=c+1; end
%------------------------define event time for event-alinged responses--------------------------

event_times = TrialTimingData(:,13); % vis stimulus onset

[Raster_MatrixStim]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);

event_times = TrialTimingData(:,14); %reward onset

[Raster_MatrixReward]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);

% --------------- wheel movement ----------------

rawPos = block.inputSensorPositions;
rawPos = wheel.correctCounterDiscont(rawPos); % correction because sometimes negative wheel positions wrap around
rawTimes = block.inputSensorPositionTimes;

% interpolate it to be regularly sampled

Fs = 1000;
t = rawTimes(1):1/Fs:rawTimes(end);
pos = interp1(rawTimes, rawPos, t, 'linear');

wheelRadius = 31; % mm (burgess wheel) (measured by Chris)
wheelRadius = 150; % mm (running wheel) (guess!)

rotaryEncoderResolution = 360*4; % number of ticks for one revolution (factor of 4 is according to CB)

pos = pos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm

nSampToPlot = min(700000, length(pos));

% make the position signal be relative to a particular time(s)
tr = [block.trial];
tr = tr(1:block.numCompletedTrials);
eventTimes = [tr.stimulusCueStartedTime];

posRel = wheel.resetAtEvent(t, pos, eventTimes);
%--------------------------------------------------------------------------


% --------- make plots --------------
%%
figure; hold on
subplot(6,2,1); % psych curve

plot(Stimz, performance,'k')

subplot(6,2,2); % psych curve

plot(Stimz, RT,'k')


subplot(4,1,2); % snapshot of CA trace
plot(downsample(TimeStamps,10),smooth(downsample(DeltaFoverF,10)),'k')

ax = gca;
for ievent = 1: size(TrialTimingData,1)
    
    h=rectangle(ax, 'Position',[TrialTimingData(ievent,13) -2 TrialTimingData(ievent,14)-TrialTimingData(ievent,13) 10],'EdgeColor',[1 1 1]);
    
    if TrialTimingData(ievent,9)==1
        h.FaceColor = [0 .5 .5 0.3];
    else
        h.FaceColor = [0.5 .1 0 0.3];
    end
end

a = downsample(TimeStamps,10);
xlim([0 a(end)])
xlabel ('Time (s)')
ylabel ('Response')
title ('Ca responses')
xlim([10 210])
ylim([-4 8])

subplot(4,1,3);
plot(t(1:end-1)+FileAlignDelay, smooth(diff(posRel)./diff(t)),'k')

xlabel ('Time (s)')
title ('Wheel acceleration')
xlim([10 210])
ylim([-200 200])


subplot(6,2,9); hold on
title ('Stimulus aligned ')
xlim([0 (stop-start)* sample_rate]/downsampleScale)

xticks([0 1200 2400 ])
xticklabels({'-1','0','1'})

c=1;
for istim = StimzAbs
    plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim,:)),'color',colorGray(c,:))
    c=c+1;
    
end
if length(Stimz)==4
    legend(num2str(Stimz(1)), num2str(Stimz(2)), num2str(Stimz(3)),num2str(Stimz(4)),'location','best')
elseif length(Stimz)==3
    legend(num2str(Stimz(1)), num2str(Stimz(2)), num2str(Stimz(3)),'location','best')
end


subplot(6,2,10); hold on
title ('Stimulus aligned')
c=1;

for istim = StimzAbs
    
    plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==1,:)),'color',colorGreen(c,:))
    plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==0,:)),'color',colorRed(c,:))
    
    c=c+1;
    
end


xlim([0 (stop-start)* sample_rate]/downsampleScale)
xticks([0 1200 2400 ])
xticklabels({'-1','0','1'})




subplot(6,2,11); hold on

xlim([0 (stop-start)* sample_rate]/downsampleScale)
xticks([0 1200 2400 ])
xticklabels({'-1','0','1'})
title ('Outcome aligned ')

c=1;

for istim = StimzAbs
    
    plot(nanmean(Raster_MatrixReward(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==1,:)),'color',colorGreen(c,:))
    plot(nanmean(Raster_MatrixReward(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==0,:)),'color',colorRed(c,:))
    
    c=c+1;
    
end


legend('no reward','reward','location','best')

xlabel ('Time (s)')

%%
figure(2); hold on
c=1;

for istim = Stimz
    AvRespCue(c) = nanmean(nanmean(Raster_MatrixStim((TrialTimingData(:,2))==istim,1800:end)));
    c=c+1;
    
end



subplot(6,2,1); hold on

plot(Stimz,AvRespCue)

c=1;

AvRespCueError = nan(1,length(Stimz));
for istim = Stimz
    AvRespCueCorrect(c) = nanmean(nanmean(Raster_MatrixStim((TrialTimingData(:,9))==1 & (TrialTimingData(:,2))==istim,1600:end)));
    AvRespCueError(c) = nanmean(nanmean(Raster_MatrixStim((TrialTimingData(:,9))==0 & (TrialTimingData(:,2))==istim,1600:end)));
    
    c=c+1;
    
end



subplot(6,2,2); hold on

plot(Stimz,AvRespCueCorrect,'g')

plot(Stimz,AvRespCueError,'r')

c=1;
AvNormRewRespError = nan(1,length(Stimz));
for istim = Stimz
    IndexCorrect = (TrialTimingData(:,9))==1 & (TrialTimingData(:,2))==istim;
    IndexError = (TrialTimingData(:,9))==0 & (TrialTimingData(:,2))==istim;
    
    AvNormRewRespCorrect(c) = nanmean(...
        nanmean(Raster_MatrixReward(IndexCorrect,1400:end),2)- nanmean( Raster_MatrixReward(IndexCorrect,200:1200),2));
    
     AvNormRewRespError(c) = nanmean(...
        nanmean(Raster_MatrixReward(IndexError,1400:end),2)- nanmean( Raster_MatrixReward(IndexError,200:1200),2));
       
       
    c=c+1;
    
end


subplot(6,2,3); hold on
title('reward response')
plot(Stimz,AvNormRewRespCorrect,'g')
plot(Stimz,AvNormRewRespError,'r')


