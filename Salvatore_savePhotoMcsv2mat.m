function Salvatore_savePhotoMcsv2mat(animal_name, exp_date, exp_series)

% Example : Salvatore_savePhotoMcsv2mat('ALK068', '2018-01-31', '1')
% load csv photoM file, conevert them to .mat structre, and save it in the same folder of the csv file
% Armin Lak 2018-02-01

[animal_ID, ~] =Salvatore_Get_chan_order(animal_name);


load('MiceExpInfoPhotoM')                                   % load beh data databse

% find path to beh and photometry data
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];

addpath (genpath(path2photoM))

TargetSessionFound = 0;
isession = 1;
while  TargetSessionFound == 0
    
    TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[exp_date,'_',exp_series,'_',animal_name,'_Block.mat']);
    
    isession = isession + 1;
end

TargetSession = isession - 1;

photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile;
photoM = readtable([path2photoM,'\',photoMFileName]);


photoMdata.AnalogIn_1=photoM.AnalogIn_1;
photoMdata.AnalogIn_2=photoM.AnalogIn_2;
photoMdata.AnalogIn_2_dF_F0=photoM.AnalogIn_2_dF_F0;
photoMdata.Time_s_=photoM.Time_s_;
photoMdata.TTL_1=photoM.TTL_1;

cd(path2photoM)

filename = [photoMFileName(1:end-4),'.mat']

save(filename, 'photoMdata')










