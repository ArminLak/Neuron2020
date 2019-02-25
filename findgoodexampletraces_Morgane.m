% get trial #s where that trial and the subsequent 7 trials are varied and
% potentially good for example Ca traces =)

clear all

% --------- enter reqs -----------------------------

animal_name = 'ALK074'
exp_date   = '2018-05-02'
exp_series ='1'


desired_contrasts = [-0.5, -0.12, 0.12, 0.25, 0.5];
trial_seq_n = 6; %must be equal to or less tha nnumber of desired contrasts 

% --------------------------------------------------

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

load('MiceExpInfoPhotoM')                                   % load beh data databse
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name); %get animal ID

path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series]; addpath (genpath(path2Beh)) %add path to beh
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM']; addpath (genpath(path2photoM)) %add path to photom

TargetSessionFound = 0;
isession = 1;
while  TargetSessionFound == 0
    TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[exp_date,'_',exp_series,'_',animal_name,'_Block.mat']); %find target session
    isession = isession + 1;
end
TargetSession = isession - 1;

TrialTimingData = MiceExpInfo.mice(animal_ID).session(TargetSession).TrialTimingData; %load trial timing data 
Stimz=unique(TrialTimingData(:,2))';                % stims for that session 

photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile(1:end-4);% load photoM data
load(photoMFileName);
DeltaFoverF = photoMdata.AnalogIn_2_dF_F0; % get DF/F trace data
TimeStamps=photoMdata.Time_s_; % get time stamps 
ds_TimeStamps = downsample(TimeStamps,10);
sm_ds_DeltaFoverF = downsample(DeltaFoverF,10);


trials2return = [];

for i = 1: length(TrialTimingData)-trial_seq_n+1
    match = 0;
    tempTrialz = unique(TrialTimingData(i:i+(trial_seq_n-1), 2)); %trials to test 
    
    for contrast = 1: length(tempTrialz)        
        if max(ismember(desired_contrasts, tempTrialz(contrast))) % is this trial's contrast in the desired contrasts list ?
            match = match + 1;
        end
    end

    if match == length(desired_contrasts)
        trials2return = [trials2return;i];
    end

end

trials2return

%% figs each have 4 example traces 
itrial = 1;

for i = 1:ceil(length(trials2return)/4)     % loop thru figs
    
    figure 
    
    for iplot = 1:4                         % loop thru subplots
        subplot(4, 1, iplot)
        
        %get time of trials2return(itrial)
        starttime = TrialTimingData(trials2return(itrial), 12) - 1;
        endtime = TrialTimingData(trials2return(itrial+trial_seq_n-1), 14) + 1;
        
        plot(ds_TimeStamps(starttime*(sample_rate/downsampleScale):endtime*(sample_rate/downsampleScale)), ...
            sm_ds_DeltaFoverF(starttime*(sample_rate/downsampleScale):endtime*(sample_rate/downsampleScale)),'color', [0.25 0.25 0.25]) %plot df/f between start -> end time 
        
        ax = gca;
        ymin = ax.YLim(1); ymax = ax.YLim(2);

        for ievent = TrialTimingData(trials2return(itrial):trials2return(itrial)+trial_seq_n-1,1)
            
            h=rectangle(ax, 'Position',[TrialTimingData(ievent,13) ymin TrialTimingData(ievent,14)-TrialTimingData(ievent,13) ymax+abs(ymin)],'EdgeColor',[1 1 1], 'FaceColor', [144/255 186/255 212/255 0.2]);
            text(TrialTimingData(ievent, 13), 2, num2str(TrialTimingData(ievent, 2)), 'FontWeight', 'bold')
            line([TrialTimingData(ievent, 12) TrialTimingData(ievent, 12)], [min(sm_ds_DeltaFoverF) max(sm_ds_DeltaFoverF)], 'Color', [74/255 127/255 189/255], 'LineStyle', '--', 'LineWidth', 1.5);
            line([TrialTimingData(ievent, 13) TrialTimingData(ievent, 13)], [min(sm_ds_DeltaFoverF) max(sm_ds_DeltaFoverF)], 'color', [74/255 127/255 189/255] , 'LineWidth', 1.5);
            rl = line([TrialTimingData(ievent, 14) TrialTimingData(ievent, 14)], [min(sm_ds_DeltaFoverF) max(sm_ds_DeltaFoverF)], 'LineWidth', 1.5);
            
            if TrialTimingData(ievent,9)==1
                rl.Color = [47/255 189/255 28/255];
                
            else
                rl.Color = [189/255 89/255 28/255];
            end
        end

        hold on

    end

end






