function [WheelTime,WheelMove] = Salvatore_AlignWheelPhotoM(animal_name, exp_date, exp_series)

% Armin Mar 2019


[animal_ID, ~] =Salvatore_Get_chan_order(animal_name);


load('MiceExpInfoPhotoM')                                   % load beh data databse

% find path to beh and photometry data
path2Data= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];

addpath (genpath(path2Data))


filename_block = [exp_date,'_',num2str(exp_series),'_',animal_name,'_Block','.mat']
%filename_Timeline = [exp_date,'_',num2str(exp_series),'_',animal_name,'_Timeline','.mat']




TargetSessionFound = 0;
isession = 1;
while  TargetSessionFound == 0
    TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[exp_date,'_',exp_series,'_',animal_name,'_Block.mat']); %find target session
    isession = isession + 1;
end
TargetSession = isession - 1;

    
%TrialTimingData = MiceExpInfo.mice(animal_ID).session(TargetSession).TrialTimingData; %load trial timing data 

load(filename_block)
AlignDelay = MiceExpInfo.mice(animal_ID).session(TargetSession).AlignDelay;

WheelTime= block.inputSensorPositionTimes(2:end) + AlignDelay;

WheelMove = abs(diff(block.inputSensorPositions));

