% clear all
% close all

animal_name = 'MMM001'
exp_dates   = {'2018-07-06', '2018-07-09', '2018-07-10_1', '2018-07-11'};
exp_series ={'1', '1', '1', '1'};


%--------------- useful information --------------------------------------
% task event
% 10: action time
% 12: beep 
% 13: stimulus 
% 14: reward
% ------------------------------------------------------------------------
% start and stop of time axis for plot (in second before and after the event)
start = -1 % s this should be -1 or less
stop = 1     % s


load('MiceExpInfoPhotoM')                                   % load beh data databse
sample_rate = 12000;                                        % photoM recording sampling rate
downsampleScale = 10;                                       % factor downsampling the Ca responses


%define time zones for normalisation after alignment

preAlign = ((sample_rate/downsampleScale)*abs(start)-800):((sample_rate/downsampleScale)*abs(start)-70);
postAlign = ((sample_rate/downsampleScale)*abs(start)+200):((sample_rate/downsampleScale)*(abs(start))+700);


%-------------------------------find path, add path and load data----------------------------------
% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);


% loop through each session in the list

for iSession = 1:numel(exp_dates)

    temp_date = expDate(iSession);
    temp_series = exp_series(iSession);

    % find path to beh and photometry data  
    path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',temp_date,'\photoM'];
    path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',temp_date,'\',temp_series];

    addpath (genpath(path2Beh))
    addpath (genpath(path2photoM))

    % load wheel data
    load([exp_date,'_',temp_series,'_',animal_name,'_Block.mat']);

end


