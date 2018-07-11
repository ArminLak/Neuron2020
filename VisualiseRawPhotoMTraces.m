
% This code plots a section of both raw photoM traces and the DeltaF/F
% Morgane Moss July 2018

clear all
close all


%[48, 50,51]  coresponding to ALK068, 70 and 71
% session 5 of ALK068 is chosen to be shown in paper figure

% animal_ID = 51

animal_name = 'ALK068'
exp_date   = '2018-01-25'
exp_series ='8'

load('MiceExpInfoPhotoM')                                   % load beh data databse

%%

start = 320 % trace time section in seconds
stop = 390
sample_rate = 12000

%-------------------------------find path, add path and load data----------------------------------
% read animals' ID
[animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name);

% find path to beh and photometry data  
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];

addpath (genpath(path2Beh))
addpath (genpath(path2photoM))

% load photoM data
load([animal_name,'_',exp_date,'_',exp_series,'_F.mat']);

AnalogIn1 = photoMdata.AnalogIn_1;
AnalogIn2 = photoMdata.AnalogIn_2;
DeltaFoverF = photoMdata.AnalogIn_2_dF_F0;
TTL = photoMdata.TTL_1;




figure; hold on

subplot(4, 1, 2);

plot(AnalogIn2(sample_rate*start : sample_rate*stop));
xlim([0 (sample_rate*stop - sample_rate*start)])
xticklabels([]);
xlabel('');
ylabel('Calcium dependent');
ax = ylim;


subplot(4, 1, 1);

plot(AnalogIn1(sample_rate*start:sample_rate*stop));
xlim([0 (sample_rate*stop - sample_rate*start)])
xticklabels([]);
xlabel('');
ylabel('Calcium independent');
ylim(ax);



subplot(4, 1, 3);

plot(DeltaFoverF(sample_rate*start : sample_rate*stop));
xlim([0 (sample_rate*stop - sample_rate*start)])
xticklabels([]);
xlabel('');
ylabel('{\Delta} F / F');


subplot(4, 1, 4);

plot(TTL(sample_rate*start : sample_rate*stop));
xlim([0 (sample_rate*stop - sample_rate*start)])
xticks(0 : 5*sample_rate : (sample_rate*stop - sample_rate*start));
yticks(0 : 1 : 1);
xticklabels(start : 5 : stop);
xlabel('Time (s)');
ylabel('Reward');














