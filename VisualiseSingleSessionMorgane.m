clear all
%close all

animal_name = 'ALK071'
exp_date   = '2018-03-06'
exp_series ='6'

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

% original data is sammple 12K per s. We downsample 10 times making it 1.2K
% per s. We define 2 s before and 2 s after the event (4 * 1.2k) we stay
% conservative by keeping 4700 of 4800 (so event is at 2320  sample)

load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses


%define time zones for normalisation after alignment

preAlign = (0.1*(sample_rate/downsampleScale)*abs(start)):((sample_rate/downsampleScale)*abs(start)-70);
postAlign = ((sample_rate/downsampleScale)*abs(start)+200):(0.9*(sample_rate/downsampleScale)*(abs(start)+stop));


%---------------- colors for plotting-------------------------------------

colorGray4 = [0.8 0.8 0.8
    0.6 0.6 0.6
    0.3 0.3 0.3
    0 0 0];
colorGray3 = [0.8 0.8 0.8
    0.4 0.4 0.4
    0 0 0];
colorGray2 = [0.7 0.7 0.7
    0 0 0];
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

Stimz=unique(TrialTimingData(:,2))';
StimzAbs=unique(abs(TrialTimingData(:,2)))';

% delay for wheel movement
FileAlignDelay = MiceExpInfo.mice(animal_ID).session(TargetSession).AlignDelay;

% load photoM data
photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile(1:end-4);
%photoMdata = readtable([path2photoM,'\',photoMFileName]);
load(photoMFileName);DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;TimeStamps=photoMdata.Time_s_;

    
%------------------------------- plot psychoimetric curve------------------
% psych curve:
Noutof = histc(TrialTimingData(:,2), Stimz)';
Nhit = nan(1, size(Noutof, 2));
Ncorrect = nan(1, size(Noutof, 2));
c = 1;
for istim = Stimz
    Nhit(c)=sum(TrialTimingData(TrialTimingData(:,2)==istim,3));
    Ncorrect(c)=sum(TrialTimingData(TrialTimingData(:,2)==istim,9));
    performance(c) = nanmean (TrialTimingData(TrialTimingData(:,2)==istim,3));
    
    c=c+1; end

% compute errorbars for psych function
pR = Nhit ./ Noutof;    
psuccess=Ncorrect ./ Noutof;
[M,V]=binostat(Noutof, psuccess)
std = sqrt(V)./100


%------------------------define event time for event-alinged responses--------------------------

event_times = TrialTimingData(:,12); % trial initiation/ onset beep

[Raster_MatrixOnsetBeep]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale);


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
%f = figure; hold on
f = figure('Position', [300 200 800 900]); hold on


% ---------- psych curve ------------------------------------------------
subplot(7,2,1); % psych curve
errorbar(unique(TrialTimingData(:,2)'), performance, std, 'o', 'MarkerSize', 1)
plot(unique(TrialTimingData(:,2)'), performance, 'color', [74/255 127/255 189/255],'LineWidth',2,'Marker','o','MarkerFaceColor', [74/255 127/255 189/255],'MarkerSize',3)
xlabel('Contrast')
xticks([min(Stimz) 0 max(Stimz)])
title('Psych Curve')
ylabel('% Right-ward')
yticklabels({0 50 100})



% --------- Ca response trace ---------------------------
subplot(4,1,2); % snapshot of CA trace

% 12: onset beep(solid line)
% 13: stimulus (dotted line)
% 14: reward (green = reward, red = no reward)


plot(downsample(TimeStamps,10),smooth(downsample(DeltaFoverF,10)),'color', [0.25 0.25 0.25])
hold on
% line([90 90], [-2 4])

ax = gca;
Visstart = 30; % visualise trace from this trial 
Visstop  = 36;  % visualise trace up to this trial
ymin = -13
ymax = 21
ylim([ymin ymax])

for ievent = Visstart: Visstop
    
    h=rectangle(ax, 'Position',[TrialTimingData(ievent,13) ymin TrialTimingData(ievent,14)-TrialTimingData(ievent,13) ymax-ymin],'EdgeColor',[1 1 1], 'FaceColor', [144/255 186/255 212/255 0.2]);
    
    line([TrialTimingData(ievent, 12) TrialTimingData(ievent, 12)], [min(smooth(downsample(DeltaFoverF, 10))) max(smooth(downsample(DeltaFoverF, 10)))], 'Color', [74/255 127/255 189/255], 'LineStyle', '--', 'LineWidth', 1.5);
    line([TrialTimingData(ievent, 13) TrialTimingData(ievent, 13)], [min(smooth(downsample(DeltaFoverF, 10))) max(smooth(downsample(DeltaFoverF, 10)))], 'color', [74/255 127/255 189/255] , 'LineWidth', 1.5);
    rl = line([TrialTimingData(ievent, 14) TrialTimingData(ievent, 14)], [min(smooth(downsample(DeltaFoverF, 10))) max(smooth(downsample(DeltaFoverF, 10)))], 'LineWidth', 1.5);
    
    
    if TrialTimingData(ievent,9)==1
        rl.Color = [47/255 189/255 28/255];
        
    else
        rl.Color = [189/255 89/255 28/255];
        
    end
end

 

a = downsample(TimeStamps,10);
xlim([0 a(end)])
xlabel ('Time (s)')
ylabel ('Response, {\Delta} F / F')
title ('Ca responses')
xlim([TrialTimingData(Visstart, 12) - 1  TrialTimingData(Visstop, 14) + 2])

% find the max and min of signal: something like this:    
%Max4Yaxis = max ( DeltaFoverF (TimeStamps (((TrialTimingData(3, 12)-1):(TrialTimingData(ievent, 14) + 2)), 1)))



% subplot(4,1,3);
% plot(t(1:end-1)+FileAlignDelay-0.2, smooth(diff(posRel)./diff(t)),'k')
% 
% xlabel ('Time (s)')
% title ('Wheel acceleration')
% xlim([10 110])
% ylim([-200 200])

%----------------------------- subplot for stumulus aligned normalised
subplot(7,2,2); hold on
title('Normalised post-stim response')
%xlim([ 0 (stop-start)* sample_rate]/downsampleScale)

% normalise trial-by-trial stim responses by the activivty before stimulus 

prepostmeandiff = nan(size(Raster_MatrixStim,1), 1); %pre-allocate

for ievent = 1:size(Raster_MatrixStim,1)
    
    prepostmeandiff(ievent) = mean(Raster_MatrixStim(ievent, postAlign) - mean(Raster_MatrixStim(ievent, preAlign))); %difference between before and after; indicator of relative signal change 
    
end
% next for each istim

normalised = nan(size(unique(TrialTimingData(:,2)),1),1);
c=1;
for istim = unique(TrialTimingData(:,2))'
    
    normalised(c)= nanmean(prepostmeandiff(TrialTimingData(:,2)==istim, 1));
    c=c+1;
end



ylim([min(normalised)-min(normalised)*.2 max(normalised)*1.2])
xticks([Stimz])
xticklabels([Stimz])
plot (Stimz, normalised', 'LineWidth',2,'Marker','o','MarkerFaceColor', [74/255 127/255 189/255],'MarkerSize',3, 'Color', [74/255 127/255 189/255])
xlabel('Contrast')
ylabel('Rel')

% 
% -------------------------------------------------subplot for stimulus aligned, not normalised
subplot(7,2,9); hold on
title ('Stimulus aligned ')
xlim([0 (stop-start)* sample_rate]/downsampleScale)
% xticks([0 1200 2400 3600 4800])
xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
xticklabels([start:1:stop])
%xticklabels({'-2','-1','0','1','2'})

c=1;
for istim = StimzAbs
    h = plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim,:)),'color',colorGray4(c,:),'LineWidth',2)
    
    if length(StimzAbs)==4
        set(h, 'color', colorGray4(c,:))
    elseif length(StimzAbs)==3
        set(h, 'color', colorGray3(c,:))
    elseif length(StimzAbs)==2
        set(h, 'color', colorGray2(c,:))
    end
    
    c=c+1;
end
    
if length(StimzAbs)==4
    l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),num2str(StimzAbs(4)),'location','best')
elseif length(StimzAbs)==3
    l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),'location','best')
elseif length(StimzAbs)==2
    l = legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), 'location', 'best')
end

title(l, 'Contrast')

ylabel('{\Delta} F / F')
xlabel ('Time (s)')


subplot(7,2,10); hold on
title ('Stimulus aligned')
c=1;

for istim = StimzAbs
    
    plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==1,:)),'color',colorGreen(c,:),'LineWidth',2)
    plot(nanmean(Raster_MatrixStim(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==0,:)),'color',colorRed(c,:),'LineWidth',2)
    
    c=c+1;
    
end


xlim([0 (stop-start)* sample_rate]/downsampleScale)
xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
% 0:((stop*start)*sample_rate)/downsampleScale):(sample_rate/downsampleScale)
xticklabels([start:1:stop])
ylabel('{\Delta} F / F')
xlabel('Time (s)')



subplot(7,2,11); hold on

xlim([0 (stop-start)* sample_rate]/downsampleScale)
xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
xticklabels([start:1:stop])
title ('Outcome aligned ')

c=1;

for istim = StimzAbs
    
    plot(nanmean(Raster_MatrixReward(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==1,:)),'color',colorGreen(c,:),'LineWidth',2)
    plot(nanmean(Raster_MatrixReward(abs(TrialTimingData(:,2))==istim & TrialTimingData(:,9)==0,:)),'color',colorRed(c,:),'LineWidth',2)
    
    c=c+1;
    
end


legend('Reward','No reward','location','best')
ylabel('{\Delta} F / F')
xlabel ('Time (s)')


c=1;

for istim = Stimz
    AvResp(c) = nanmean(nanmean(Raster_MatrixStim((TrialTimingData(:,2))==istim,:)));
    c=c+1;
    
end

subplot(7,2,12); hold on
xlim([0 (stop-start)* sample_rate]/downsampleScale)
xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
xticklabels([start:1:stop])
% ylim([-1 1])
yticks([-1 0 1])
title('Onset Beep Aligned')
ylabel('{\Delta} F / F')
xlabel('Time (s)')

plot(nanmean(Raster_MatrixOnsetBeep), 'color', colorGray4(1,:),'LineWidth',2)


%------------------------------------------ post stim response by contrast for errors and correct trials-----
subplot(7, 2, 13); hold on
title('Normalised post stim response') %separate line for correct and error trials

% pre-allocate; col 1 = correct, col 2 = error
poststim = nan(length(normalised), 2);

c = 1;

for istim = unique(TrialTimingData(:,2))'
    
    poststim(c, 1)= nanmean(prepostmeandiff(TrialTimingData(:,2)==istim & TrialTimingData(:,9)==1, 1));
    poststim(c, 2)= nanmean(prepostmeandiff(TrialTimingData(:,2)==istim & TrialTimingData(:,9)==0, 1));
    c=c+1;
end

plot(Stimz, poststim(:,1), 'color', [47/255 189/255 28/255], 'LineWidth', 2)
hold on;
plot(Stimz, poststim(:,2), 'color', [189/255 89/255 28/255], 'LineWidth', 2)

ylabel('{\Delta} F / F')
xlabel('Contrast')
xticks([Stimz])
legend('Reward', 'No reward')

%------------------------------- normalised post outcome response -------------

subplot(7, 2, 14); hold on
title('Normalised post outcome response')

outcomediff = nan(size(Raster_MatrixReward,1), 1); %pre-allocate

for ievent = 1:size(Raster_MatrixReward,1)
    
    outcomediff(ievent) = mean(Raster_MatrixReward(ievent, postAlign) - mean(Raster_MatrixReward(ievent, preAlign))); %difference between before and after; indicator of relative signal change 
    
end

postoutcome = nan(length(normalised), 2); %col 1 = correct, col 2 = error
c=1;

for istim = unique(TrialTimingData(:,2))'
    
    postoutcome(c, 1)= nanmean(outcomediff(TrialTimingData(:,2)==istim & TrialTimingData(:,9)==1, 1));
    postoutcome(c, 2)= nanmean(outcomediff(TrialTimingData(:,2)==istim & TrialTimingData(:,9)==0, 1));
    c=c+1;
end

plot(Stimz, postoutcome(:,1), 'color', [47/255 189/255 28/255], 'LineWidth', 2)
hold on;
plot(Stimz, postoutcome(:,2), 'color', [189/255 89/255 28/255], 'LineWidth', 2)

ylabel('{\Delta} F / F')
xlabel('Contrast')
xticks([Stimz])
legend('Reward', 'No reward')


% c=1;
% for istim = StimzAbs
%    plot(nanmean(Raster_MatrixOnsetBeep(abs(TrialTimingData(:,2))==istim,:)),'color',colorGray(c,:));
%    c=c+1;
    
% end
