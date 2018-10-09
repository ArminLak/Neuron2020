% Visualise avg responses to stimulus and reward in specified sessions

% Morgane August 2018


% to do instead of z-scoring normalise to max response 

clear all
close all

<<<<<<< HEAD
animal_name = 'MMM002'

z = 'n' % option to z-score data


if animal_name == 'ALK068'
    SessionList = [1:11];
    ylimrwd = [-2 2];
    ylimstim = [-1 2.5];
    load('BehPhotoM_Exp7_VTA')                                   % load beh data databseB
    
elseif animal_name == 'ALK070'
    SessionList = [1:9]; %weaker signal in sessions 10 and 11
    ylimrwd = [-2 3];
    ylimstim = [-1 1.5];
    load('BehPhotoM_Exp7_VTA')                                   % load beh data databse
    
elseif animal_name == 'ALK071'
    SessionList = [1:9];
    ylimrwd = [-8 15];
    ylimstim = [-3 6];
    load('BehPhotoM_Exp7_VTA')                                   % load beh data databse
    
elseif animal_name == 'ALK074'
    SessionList = [1:11];
    ylimrwd = [-5 15];
    ylimstim = [-5 15];
    load('BehPhotoM_Exp7_DMS')
     
elseif animal_name == 'ALK075'
    SessionList = [1:11];
    ylimrwd = [-0.5 1];
    ylimstim = [-.5 1];
    load('BehPhotoM_Exp7_DMS')
    
elseif animal_name == 'MMM001'
    SessionList = [1:11];
    ylimrwd = [-1 2];
    ylimstim = [-1 2];
    load('BehPhotoM_Exp7_NAc')
    
elseif animal_name == 'MMM002'
    SessionList = [1:11];
    ylimrwd = [-2 5];
    ylimstim = [-2 5];
    load('BehPhotoM_Exp7_NAc')
end
=======
% animal_name= ['ALK068'];
animal_list= [{'ALK068', 'ALK070', 'ALK071'}];
% animal_ID_list = [48, 50, 51];
SessionList = [1:9];
load('BehPhotoM_Exp7_VTA')

z = 'y' %#ok<NOPTS> % option to z-score data / must be y if averaging animals 


% if strcmp(animal_name,'ALK068')
%     SessionList = [1:11];
%     ylimrwd = [-3 3];
%     ylimstim = [-1 2.5];
%     load('BehPhotoM_Exp7_VTA')                                   % load beh data databseB
%     
% elseif strcmp(animal_name,'ALK070')
%     SessionList = [1:9]; %weaker signal in sessions 10 and 11 
%     ylimrwd = [-2 3];
%     ylimstim = [-1 1.5];
%     load('BehPhotoM_Exp7_VTA')                                   % load beh data databse
%     
% elseif strcmp(animal_name, 'ALK071')
%     SessionList = [1:9];
%     ylimrwd = [-8 25];
%     ylimstim = [-4 11];
%     load('BehPhotoM_Exp7_VTA')                                   % load beh data databse
% 
% elseif strcmp(animal_name,'ALK074')
%     SessionList = [1:11];
%     ylimrwd = [-5 15];
%     ylimstim = [-5 15];
%     load('BehPhotoM_Exp7_DMS')
% 
% elseif strcmp(animal_name, 'ALK075')
%     SessionList = [1:11];
%     ylimrwd = [-0.5 1];
%     ylimstim = [-.5 1];
%     load('BehPhotoM_Exp7_DMS')
% 
% elseif strcmp(animal_name, 'MMM001')
%     SessionList = [1:11];
%     ylimrwd = [-1 2];
%     ylimstim = [-1 2];
%     load('BehPhotoM_Exp7_NAc')
% end
>>>>>>> fde2c5f41f3f69eb704c36aa00875984951200ad



if z == 'y'
    ylimstim = [-2 3];
    ylimrwd = [-3 3];
end


% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = -1 % s this should be -1 or less
stop = 2    % s
event_time = 3; % this is the time in the summary matrix where the event took place

sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses

start2stop = (event_time+start)*sample_rate/downsampleScale:(event_time+stop)*sample_rate/downsampleScale; %window of interest



% event is at 3 seconds i.e. point 3600

preAlignStim = (event_time-0.4)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignStim = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;
preAlignRwd = (event_time-0.2)*sample_rate/downsampleScale : (event_time-0)*sample_rate/downsampleScale;
postAlignRwd = (event_time+0.2)*sample_rate/downsampleScale : (event_time+0.8)*sample_rate/downsampleScale;


% plot colours
colorGray4 = [0.8 0.8 0.8 %lightest
    0.6 0.6 0.6
    0.4 0.4 0.4
    0.2 0.2 0.2
    0 0 0]; % black
colorGray3 = [0.8 0.8 0.8
    0.4 0.4 0.4
    0 0 0];
colorGray2 = [0.7 0.7 0.7
    0 0 0];
colorGreen = [0 1 0
    0 0.8 0
    0 0.6 0
    0  0.3 0];


% -------- prep for summary figures --------

StimzAbs = [1.0 0.5 0.25 0.12 0];
StimResponseAvg = nan(length(StimzAbs)*length(animal_list)*length(SessionList), length(start2stop)) ;
RwdResponseAvg = nan(3*length(animal_list)*length(SessionList), length(start2stop)) ;


% ---- grand summary matrix for avg across animals figure ----
% e.g. 3 animals, 5 stim levels =>
    % every 5 rows is one animal, one session
    % every 15 rows is all 3 animals, one session 
    % every 5 rows correspond to the stim levels 
    
stimStartRow = 1:length(StimzAbs):length(StimzAbs)*length(animal_list);
rwdStartRow = 1: 3 : 3*length(animal_list);

for animalcount = 1:length(animal_list)
    animal_name = animal_list(animalcount);
    [animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

    
%     SessionList = 1:length(BehPhotoM(animal_ID).Session);
    stimRow = stimStartRow(animalcount);
    rwdRow = rwdStartRow(animalcount);
    for iSession = SessionList

        TargetSession = iSession;

        TrialTimingData = BehPhotoM(animal_ID).Session(TargetSession).TrialTimingData;
        TrialTimingDataCor = TrialTimingData(TrialTimingData(:,9)==1, :);

        NeuronStim = BehPhotoM(animal_ID).Session(TargetSession).NeuronStim;
        NeuronStimCor = NeuronStim(TrialTimingData(:,9)==1, :);

        NeuronReward = BehPhotoM(animal_ID).Session(TargetSession).NeuronReward;
        NeuronRewardCor = NeuronReward(TrialTimingData(:,9)==1, :);

        for stimcount = 1:length(StimzAbs) % fat stim matrix
            istim = StimzAbs(stimcount);
            StimResponseAvg(stimRow,:) = zscore(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==istim, start2stop)));

            stimRow = stimRow + 1;
        end
        stimRow = stimRow + 10;
            
        % rwd matrix : correct easy/hard (2 separate rows)
        if numel(unique(abs(TrialTimingData(:, 2)))) > 2
            correcteasy = abs(TrialTimingData(:,2)==[1, 0.5, 0.25]);
            correcthard = abs(TrialTimingData(:,2)==[ 0.12, 0]);
            RwdResponseAvg(rwdRow,:) = zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop)));
            rwdRow = rwdRow+1;
            RwdResponseAvg(rwdRow,:) = zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop)));
        elseif numel(unique(abs(TrialTimingData(:, 2)))) == 2
            correcteasy = abs(TrialTimingData(:,2)==[1]);
            correcthard = abs(TrialTimingData(:,2)==[0.5, 0.25, 0.12, 0]);
            RwdResponseAvg(rwdRow,:) = zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop)));
            rwdRow = rwdRow + 1;
            RwdResponseAvg(rwdRow,:) = zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop)));
        else 
            RwdResponseAvg(rwdRow,:) = zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1, start2stop)));
            rwdRow = rwdRow+1;
        end
        rwdRow = rwdRow + 1;
        % error trials : separate row 
        RwdResponseAvg(rwdRow,:) = zscore(nanmean(NeuronReward(TrialTimingData(:,9)==0, start2stop)));
           
            rwdRow = rwdRow + 7;
        
    end
        
end


%%
% ------ figure for avg across animals -------------------

<<<<<<< HEAD
% ---- start figure ----
figure;
=======
figure; hold on
>>>>>>> fde2c5f41f3f69eb704c36aa00875984951200ad
plotnum = 1;
title('average across animals')

stimSess = 1 : length(StimzAbs) : length(StimzAbs)*numel(animal_list);
rwdSess = 1 : 3 : 3*numel(animal_list);


for isession = 1:length(SessionList)

    c = 5;
    
    % first stim on LHS 
    subplot(length(SessionList), 2, plotnum)
    for stimcount = 1:length(StimzAbs) % stim resp averages plots 
        hold on;
        plot(nanmean(StimResponseAvg(stimSess, :)),'color',colorGray4(c,:),'LineWidth',2)
        stimSess = stimSess + 1;
        c = c-1;
        
        ylim(ylimstim)
        xlim([0 (stop-start)* sample_rate]/downsampleScale)
        xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
        xticklabels([''])          
    
        if plotnum == length(SessionList)*2 - 1
            ylabel('Z({\Delta} F / F)')
            xlabel ('Time (s)')
            xticklabels([start:1:stop])
        end
        
    end
    
    plotnum = plotnum + 1; % now reward RHS -----------

    subplot(length(SessionList), 2, plotnum)
    hold on
        plot(nanmean(RwdResponseAvg(rwdSess, :)), 'color', colorGreen(3,:), 'LineWidth', 2)
        plot(nanmean(RwdResponseAvg(rwdSess+1, :)), 'color', colorGreen(2,:), 'LineWidth', 2)
        plot(nanmean(RwdResponseAvg(rwdSess+2, :)), 'red', 'LineWidth', 2)
    ylim(ylimrwd)
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])          

    if plotnum == length(SessionList)*2
        ylabel('Z({\Delta} F / F)')
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
        
    end % ----------- end of reward RHS
    
    plotnum = plotnum + 1;
    
stimSess = stimSess + 10;
rwdSess = rwdSess + 9;

end



%% for one animal only 
% loop through each session in the list

for iSession = SessionList
    
<<<<<<< HEAD
    TargetSession = iSession
    
=======
    TargetSession = iSession;
 
>>>>>>> fde2c5f41f3f69eb704c36aa00875984951200ad
    % load Beh data and photometry data
    TrialTimingData = BehPhotoM(animal_ID).Session(TargetSession).TrialTimingData;
    TrialTimingDataCor = TrialTimingData(TrialTimingData(:,9)==1, :);
    
    NeuronStim = BehPhotoM(animal_ID).Session(TargetSession).NeuronStim;
    NeuronStimCor = NeuronStim(TrialTimingData(:,9)==1, :);
    
    NeuronReward = BehPhotoM(animal_ID).Session(TargetSession).NeuronReward;
    NeuronRewardCor = NeuronReward(TrialTimingData(:,9)==1, :);
    
    %StimzAbs=unique(abs(TrialTimingData(:,2)))';
    StimzAbs = [0 0.12 0.25 0.5 1.0];
    
    % ------------------- Grans summary ------------------------------
    
    subplot(length(SessionList), 2, plotnum ); hold on
    
    
    %  c=4;
    
    
    for stimcount = 1:length(StimzAbs)
        stimcolumn = 6 - stimcount; % % column 1 is 0 contrast and column 5 is 1.0 contrast
        istim = StimzAbs(stimcount);
        
        if z == 'y'
            h = plot(zscore(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==istim, start2stop))),'color',colorGray4(stimcount,:),'LineWidth',2)
        elseif z == 'n'
            h = plot(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==istim, start2stop)),'color',colorGray4(stimcount,:),'LineWidth',2)
        end
        
        StimResponseSummary(TargetSession,stimcolumn)  = nanmean(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), postAlignStim)))...
            - nanmean(nanmean(NeuronStimCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), preAlignStim))); %difference between before and after; indicator of relative signal change
        RwdResponseSummary(TargetSession,stimcolumn)  = nanmean(nanmean(NeuronRewardCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), postAlignRwd))) ...
            - nanmean(nanmean(NeuronRewardCor(abs(TrialTimingDataCor(:,2))==StimzAbs(stimcolumn), preAlignRwd))); %difference between before and after; indicator of relative signal change
        
        
        
        if length(StimzAbs)==4
            set(h, 'color', colorGray4(stimcount,:));
        elseif length(StimzAbs)==3
            set(h, 'color', colorGray3(stimcount,:));
        elseif length(StimzAbs)==2
            set(h, 'color', colorGray2(stimcount,:));
        elseif length(StimzAbs)==1
            set(h, 'color', colorGray2(2,:));
        end
        
        %     c=c-1;
    end
    
    
    
    
    if plotnum == 1
        title('Stimulus response')
    end
    
    if plotnum == 2*length(SessionList)-1
        if length(StimzAbs)==4
            legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),num2str(StimzAbs(4)),'location','northeast');
        elseif length(StimzAbs)==3
            legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), num2str(StimzAbs(3)),'location','northeast');
        elseif length(StimzAbs)==2
            legend(num2str(StimzAbs(1)), num2str(StimzAbs(2)), 'location', 'northeast');
        end
    end
    
    ylim(ylimstim)
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])
    ylabel('{\Delta} F / F')
    
    if plotnum == length(SessionList)*2 - 1
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
    end
    
    plotnum = plotnum + 1;
    
    
    
    
    
    % ------------------- reward plot --------------------------------
    
    subplot(length(SessionList), 2, plotnum ); hold on
    
    
<<<<<<< HEAD
    
    if numel(StimzAbs) > 2
        correcteasy = abs(TrialTimingData(:,2)==[1, 0.5]);
        correcthard = abs(TrialTimingData(:,2)==[0.25, 0.12, 0]);
=======
   if numel(unique(abs(TrialTimingData(:, 2)))) > 2
        correcteasy = abs(TrialTimingData(:,2)==[1, 0.5, 0.25]);
        correcthard = abs(TrialTimingData(:,2)==[ 0.12, 0]);
>>>>>>> fde2c5f41f3f69eb704c36aa00875984951200ad
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop)), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop)), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop))), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop))), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim
        end
        
<<<<<<< HEAD
    elseif numel(StimzAbs) == 2
        correcteasy = abs(TrialTimingData(:,2)==[1]);
        correcthard = abs(TrialTimingData(:,2)==[0.5, 0.25, 0.12, 0]);
=======
   elseif numel(unique(abs(TrialTimingData(:, 2)))) == 2
        correcteasy = abs(TrialTimingData(:,2)==[1, 0.5, 0.25]);
        correcthard = abs(TrialTimingData(:,2)==[ 0.12, 0]);
>>>>>>> fde2c5f41f3f69eb704c36aa00875984951200ad
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop)), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop)), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcteasy,2)==1, start2stop))), 'color', colorGreen(3,:), 'LineWidth', 2)  % plot for stim = max stim
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1 & sum(correcthard,2)==1, start2stop))), 'color', colorGreen(2,:), 'LineWidth', 2)  % plot for stim = min stim
        end
        
<<<<<<< HEAD
    elseif numel(StimzAbs) == 1
=======
   elseif numel(unique(abs(TrialTimingData(:, 2)))) == 1
>>>>>>> fde2c5f41f3f69eb704c36aa00875984951200ad
        if z == 'n'
            plot(nanmean(NeuronReward(TrialTimingData(:,9)==1, start2stop)), 'color', colorGreen(3,:), 'LineWidth', 2)
        elseif z == 'y'
            plot(zscore(nanmean(NeuronReward(TrialTimingData(:,9)==1, start2stop))), 'color', colorGreen(3,:), 'LineWidth', 2)
        end
    end
    
    hold on;
    plot(nanmean(NeuronReward(TrialTimingData(:,9)==0, start2stop)), 'red', 'LineWidth', 2) % plot error trials
    
    ylim(ylimrwd)
    xlim([0 (stop-start)* sample_rate]/downsampleScale)
    xticks([0:(sample_rate/downsampleScale):((stop-start)* sample_rate/downsampleScale)])
    xticklabels([''])
    
    if plotnum == 2
        title('Reward response')
        legend('Easy & correct', 'Hard & correct', 'Error', 'location', 'northeast')
    end
    
    if plotnum == length(SessionList)*2
        xlabel ('Time (s)')
        xticklabels([start:1:stop])
    end
    
    plotnum = plotnum + 1;
    
    
end
%%
% ---- plot summaries ----

figure; hold on

subplot(1, 2, 1) % stim first
c = 1;
for i = 1:size(StimResponseSummary, 2)
    
    plot(StimResponseSummary(:, i), 'o-', 'color', colorGray4(c,:), 'markerFaceColor', colorGray4(c,:), 'MarkerSize', 3.5, 'LineWidth', 1)
    
    c = c+1;
    
    hold on
    
end
c = 1;

subplot(1, 2, 2) % reward second

for i = 1:size(StimResponseSummary, 2)
    
    plot(RwdResponseSummary(:, i), 'o-', 'color', colorGray4(c,:), 'markerFaceColor', colorGray4(c,:), 'MarkerSize', 3.5, 'LineWidth', 1)
    
    c = c+1;
    
    hold on
    
end


