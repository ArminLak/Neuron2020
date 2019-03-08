function [WheelTime,WheelMove] = Salvatore_AlignWheelPhotoM(animal_name, exp_date, exp_series)

% Armin Mar 2019


[animal_ID, ~] =Salvatore_Get_chan_order(animal_name);


load('MiceExpInfoPhotoM')                                   % load beh data databse

% find path to beh and photometry data
path2Data= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];

addpath (genpath(path2photoM))


AlignDelay = MiceExpInfo.mice(animal_ID).session.AlignDelay;

WheelTime= block.inputSensorPositionTimes(2:end) + AlignDelay;

WheelMove = abs(diff(block.inputSensorPositions));

