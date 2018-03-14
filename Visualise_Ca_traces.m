clear all
% close all

animal_name = 'ALK071'
exp_date   = '2018-03-13'
exp_series ='1'

trial_from = 1          % visualise trace from this trial
trial_to = 50          % visualise trace up to this trial


% NB trial_to-trial_from = 50


%--------------- useful information --------------------------------------
% task event
% 10: action time
% 12: beep 
% 13: stimulus 
% 14: reward
% ------------------------------------------------------------------------

load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses


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

%---------------------------------------------------------------------------------------------------
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
load(photoMFileName);
DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
TimeStamps=photoMdata.Time_s_;


%% 

% --------- make plots --------------

% trial_from = 1          % visualise trace from this trial
% trial_to = 51          % visualise trace up to this trial
% 
% 
% % NB trial_to-trial_from = 50

% 12: onset beep(solid line)
% 13: stimulus (dotted line)
% 14: reward (green = reward, red = no reward)


% --------- Ca response trace ---------------------------
f = figure('Position', [300 200 1200 900]); hold on

formatSpec1 = '%s %s Trials %d to %d';
formatSpec2 = 'Trials %d to %d';

i = 1;

for i = 1 : round((trial_to-trial_from)./10)
    
    subplot(round((trial_to-trial_from)./10), 1, i); % snapshot of CA trace

    plot(downsample(TimeStamps,10),smooth(downsample(DeltaFoverF,10)),'color', [0.25 0.25 0.25])
    hold on
    % line([90 90], [-2 4])

    ax = gca;
    ymin = -15
    ymax = 27
    ylim([ymin ymax])
    xlim ([TrialTimingData(min(10*(i - 1) + 1:i*10), 12) - 1  TrialTimingData(max(10*(i - 1) + 1:i*10), 14) + 2])

            for ievent = min(10*(i - 1) + 1:i*10):max(10*(i - 1) + 1:i*10)

                h=rectangle(ax, 'Position',[TrialTimingData(ievent,13) ymin TrialTimingData(ievent,14)-TrialTimingData(ievent,13) ymax+abs(ymin)],'EdgeColor',[1 1 1], 'FaceColor', [144/255 186/255 212/255 0.2]);
                text(TrialTimingData(ievent, 13), 30, num2str(TrialTimingData(ievent, 2)), 'FontWeight', 'bold')
                line([TrialTimingData(ievent, 12) TrialTimingData(ievent, 12)], [min(smooth(downsample(DeltaFoverF, 10))) max(smooth(downsample(DeltaFoverF, 10)))], 'Color', [74/255 127/255 189/255], 'LineStyle', '--', 'LineWidth', 1.5);
                line([TrialTimingData(ievent, 13) TrialTimingData(ievent, 13)], [min(smooth(downsample(DeltaFoverF, 10))) max(smooth(downsample(DeltaFoverF, 10)))], 'color', [74/255 127/255 189/255] , 'LineWidth', 1.5);
                rl = line([TrialTimingData(ievent, 14) TrialTimingData(ievent, 14)], [min(smooth(downsample(DeltaFoverF, 10))) max(smooth(downsample(DeltaFoverF, 10)))], 'LineWidth', 1.5);


                if TrialTimingData(ievent,9)==1
                    rl.Color = [47/255 189/255 28/255];

                else
                    rl.Color = [189/255 89/255 28/255];

                end
            end

    
    
            
            if i ==1 
                text(TrialTimingData(trial_from,13), 35, sprintf(formatSpec1,animal_name, exp_date, trial_from, ievent), 'FontWeight', 'bold', 'FontSize', 10);
            else
                text(TrialTimingData(ievent-9,13), 35, sprintf(formatSpec2, ievent-9, ievent), 'FontWeight', 'bold', 'FontSize', 10);
            end
            
    
%     text(TrialTimingData(14,13), 35, 'ALK070 2018-03-06', 'FontWeight', 'bold', 'FontSize', 10);
    a = downsample(TimeStamps,10);
    xlim([0 a(end)])
    xlabel ('Time (s)')
    ylabel ('{\Delta} F / F')
    %title ('ALK070 2018-03-06')
    if trial_from ==1  && i==1 
     
        xlim ([0  TrialTimingData(max(10*(i - 1) + 1:i*10), 14) + 2])
    else 
        xlim ([TrialTimingData(min(10*(i - 1) + 1:i*10), 12) - 1  TrialTimingData(max(10*(i - 1) + 1:i*10), 14) + 2])
            end
    
    
        
    %xlim is 1 : 10 + (i-1)*10

    
    i = i+1;
    
end

















