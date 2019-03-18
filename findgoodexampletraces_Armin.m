% get trial #s where that trial and the subsequent 7 trials are varied and
% potentially good for example Ca traces =)

clear all
close all

% --------- enter reqs -----------------------------

% animal_name = 'MMM009'
% exp_date   = '2019-03-12'
% exp_series ='2';

animal_name = 'ALK068'
exp_ID   = '23'

desired_contrasts = [-0.5, -0.25, 0, 0.12, 0.25, 0.5];
trial_seq_n = 7; %must be equal to or more tha nnumber of desired contrasts

% --------------------------------------------------

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses
scale = sample_rate/downsampleScale;                        % for use in plotting

load('MiceExpInfoPhotoM')                                   % load beh data databse
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name); %get animal ID


%path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series]; addpath (genpath(path2Beh)) %add path to beh
%path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM']; addpath (genpath(path2photoM)) %add path to photom


path2data = ['\\zubjects.cortexlab.net\Subjects\',animal_name];
addpath(genpath(path2data))

[SessionList] = getSessionList_photoM(animal_name, exp_ID);


for iSession = SessionList

 

if isfield (MiceExpInfo.mice(animal_ID).session(iSession), 'Chan4') % if mouse is bi hem
    if MiceExpInfo.mice(animal_ID).session(iSession).Chan2(1) == 'L'
        chan_ori = 'LR';
    elseif MiceExpInfo.mice(animal_ID).session(iSession).Chan2(1) == 'R'
        chan_ori = 'RL';
    end
else
    chan_ori = 1; %if mouse is uni hem
end

exp_date   = MiceExpInfo.mice(animal_ID).session(iSession).Blockname(1:10);
exp_series = MiceExpInfo.mice(animal_ID).session(iSession).Blockname(12);

% extract aligned wheel trace : 
[WheelTime,WheelMove] = Salvatore_AlignWheelPhotoM(animal_name, exp_date, exp_series);


TrialTimingData = MiceExpInfo.mice(animal_ID).session(iSession).TrialTimingData; %load trial timing data

%desired_contrasts =unique(TrialTimingData(:,2))';
Stimz=unique(TrialTimingData(:,2))';                % stims for that session

photoMFileName=MiceExpInfo.mice(animal_ID).session(iSession).Neuronfile(1:end-4);% load photoM data
load(photoMFileName);
% TimeStamps=photoMdata.Time_s_; % get time stamps

TimeStamps = photoMdata.Time_s_;
ds_TimeStamps = downsample(TimeStamps, 10);

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

trials2return %#ok<NOPTS>

trials2return(trials2return==1)=[]; 

%% figs each have 4 example traces

for ihem = 1:numel(chan_ori)
    
    itrial = 0;
    
    if ihem == 1
        DeltaFoverF       = photoMdata.AnalogIn_2_dF_F0; % get DF/F trace data
        sm_ds_DeltaFoverF = smooth(downsample(DeltaFoverF,10),20);
    elseif ihem == 2
        DeltaFoverF       = [];
        DeltaFoverF       = photoMdata.AnalogIn_4_dF_F0; % get DF/F trace data
        sm_ds_DeltaFoverF = smooth(downsample(DeltaFoverF,10),20);
    end
    
    for ifig = 1:ceil(length(trials2return)/4)     % loop thru figs
        
        figure
        
        for iplot = 1:4 % loop thru subplots
            
            itrial = itrial +1;
            if itrial > length(trials2return)
                continue
            end
            
            subplot(4, 1, iplot)
            
            start_s = floor(TrialTimingData(trials2return(itrial), 12)) - 1;
            end_s   = ceil(TrialTimingData(trials2return(itrial)+trial_seq_n-1, 14)) + 1;
            
%             plot(ds_TimeStamps,sm_ds_DeltaFoverF,'color', [0.25 0.25 0.25]) %plot df/f between start -> end time
            
            ymax = 1.9 * max(sm_ds_DeltaFoverF(scale*start_s:scale*end_s));
            ymin = 0.4 * min(sm_ds_DeltaFoverF(scale*start_s:scale*end_s));
            
            hold on
            
            WheelStart= find(floor(100*WheelTime)/100==start_s,1);
            WheelEnd  = find(floor(100*WheelTime)/100==end_s,1);
            
            %   WheelMove(WheelMove<1)=0;
            
            
%             plot(WheelTime(WheelStart:WheelEnd),smooth(WheelMove(WheelStart:WheelEnd),8)-1,'color',[0.8 0 0.8])
%             %     plot(WheelTime,WheelMove,'b')
            
            
            %             ax = gca;
            %             ymin = ax.YLim(1); ymax = ax.YLim(2);
            yyaxis left
            ymin=ymin-3;
            ylim([ymin ymax])
            xlim([start_s  end_s])
            xticks([start_s:5:end_s])
            xticklabels([start_s:5:end_s])
            ax = gca;
            
            if iplot == 4 || itrial == length(trials2return)
                xlabel ('Time (s)')
                ylabel ('{\Delta} F / F')
            end
            
            if numel(chan_ori)==2 && iplot==1
                text(TrialTimingData(trials2return(itrial),12)-4, 1.4*ymax, [chan_ori(ihem), ' HEM', animal_name], 'FontWeight', 'bold')
            end
            
            for ievent = trials2return(itrial):trials2return(itrial)+trial_seq_n-1
                
                h=rectangle(ax, 'Position',[TrialTimingData(ievent,13) ymin (TrialTimingData(ievent,14)-TrialTimingData(ievent,13)) ymax+abs(ymin)], ...
                    'EdgeColor',[144/255 186/255 212/255 0.2], 'FaceColor', [144/255 186/255 212/255 0.2]); % trial in progress light blue rectangle
                
                text(TrialTimingData(ievent,13), 1.2*ymax, num2str(TrialTimingData(ievent, 2)), 'FontWeight', 'bold')
                line([TrialTimingData(ievent,12) TrialTimingData(ievent,12)], [ymin ymax], 'Color', [74/255 127/255 189/255], 'LineStyle', '--', 'LineWidth', 1.5);
                line([TrialTimingData(ievent,13) TrialTimingData(ievent,13)], [ymin ymax], 'Color', [74/255 127/255 189/255], 'LineWidth', 1.5); % stim line
                rl = line([TrialTimingData(ievent,14) TrialTimingData(ievent,14)], [ymin ymax], 'LineWidth', 1.5); %rwd line
                
                if TrialTimingData(ievent,9)==1
                    rl.Color = [47/255 189/255 28/255];
                    
                else
                    rl.Color = [189/255 89/255 28/255];
                end
            end
            
            plot(ds_TimeStamps,sm_ds_DeltaFoverF,'LineWidth', 1.2, 'color', [0.25 0.25 0.25]) %plot df/f
            
            yyaxis right 
            plot(WheelTime(WheelStart:WheelEnd),smooth(WheelMove(WheelStart:WheelEnd),8)-1, 'LineWidth', 1.2, 'color',[0.8 0 0.8])
            ylabel('Wheel position')
            %     plot(WheelTime,WheelMove,'b')
            
            hold on
            
        end
        
        cd('\\zserver.cortexlab.net\Lab\Share\Lak\Morgane\DMS Traces\autosave')
        
        figname=num2str(ifig);
        
        if exist([animal_name,'_',exp_date,'_',exp_series,'_fig',figname,'.fig'], 'file') > 1
           figname = num2str(  exist([animal_name,'_',exp_date,'_',exp_series,'_fig',figname,'.fig'], 'file')+1);
        end

        savefig([animal_name,'_',exp_date,'_',exp_series,'_fig',num2str(figname)])
    end
end

end
