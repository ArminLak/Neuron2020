


clear all
% close all

animals = ['MMM009', 'MMM010'] % ['MMM009', 'MMM010']

% -------------------------------------------------------------


vis_from    = -1; 
vis_to      = 1.5;

exp_date = '2019-04-04';

if length(animals) == 6
    nanimals = 1;
elseif length(animals) == 12
    nanimals = 2;
end

load('MiceExpInfoPhotoM')
samplingRate = 12000;
downsampleScale = 10;

airPuffMatrix = [];
noiseOnlyMatrix = [];

for ianimal = 1:nanimals
    if nanimals == 1
        animal_name = animals;
    elseif nanimals == 2 && ianimal == 1
        animal_name = animals(1:6);
    elseif nanimals == 2 && ianimal == 2
        animal_name = animals(7:12);
    end
    
    if strcmp(animal_name, 'MMM009')
        Session = 37;
        animal_ID = 71;
        exp_series = '3';
    elseif strcmp(animal_name, 'MMM010')
        Session = 26;
        animal_ID = 72;
        exp_series = '2';
    end
    
    % find path to beh and photometry data
    path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
    path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];
    
    addpath (genpath(path2Beh))
    addpath (genpath(path2photoM))
    
    % load photoM data
    photoMFileName=MiceExpInfo.mice(animal_ID).session(Session).Neuronfile(1:end-4);
    %photoMdata = readtable([path2photoM,'\',photoMFileName]);
    load(photoMFileName);
    DeltaFoverF     = [];
    DeltaFoverF     = photoMdata.AnalogIn_2_dF_F0;
    TTL             = [];
    TTL             = photoMdata.TTL_1;
    
    
    % get times of air puff and control noise
    airUpDown     = [];
    airUpDown     = diff(TTL);
    airPuff       = [];
    airPuff       = find(airUpDown == 1);
    noiseOnly     = [];
    noiseOnly     = airPuff + (8*samplingRate); % 8 is number of seconds between air puff and noise events always
    
    
    % create raster matrices
    
    for i = 1:25 % length(airPuff)
        airPuffMatrix      = [airPuffMatrix; (smooth(downsample(DeltaFoverF(airPuff(i)   - 5*samplingRate : airPuff(i)   + 5*samplingRate), downsampleScale), 20))']; %get 5s either side
        noiseOnlyMatrix    = [noiseOnlyMatrix; (smooth(downsample(DeltaFoverF(noiseOnly(i) - 5*samplingRate : noiseOnly(i) + 5*samplingRate), downsampleScale), 20))'];
    end
    
end

%%
% plot

figure; hold on



plot(nanmean(airPuffMatrix), 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5)
xlim ([(5+vis_from)*samplingRate/downsampleScale (5+vis_to)*samplingRate/downsampleScale])
xticks([(5+vis_from) 5 (5+vis_to)]*samplingRate/downsampleScale)

hold on;
% 
% subplot(2,1,2) % noise only
% 
plot(nanmean(noiseOnlyMatrix), 'Color', [1 0 0], 'LineWidth', 1.5)
xlim([(5+vis_from)*samplingRate/downsampleScale (5+vis_to)*samplingRate/downsampleScale])
xticks([(5+vis_from) 5 (5+vis_to)]*samplingRate/downsampleScale)
xticklabels([vis_from 0 vis_to])
ylim([-0.5 2])
ylabel ('{\Delta} F / F')
xlabel ('Time (s)')

