function TrialTimingData = Salvatore_PrepareTrialTimingDataPhotoM (animal_name, exp_date, exp_series)

% Example : TrialTimingData = Salvatore_PrepareTrialTimingDataPhotoM ('ALK068', '2018-01-29', '8');
% Armin Lak , London Dec 2017


load('MiceExpInfoPhotoM.mat')                                                               % load beh data databse
sample_rate = 12000;                                                              % Neuronal recording sampling rate


[animal_ID, ~] =Salvatore_Get_chan_order(animal_name);


%path2photoM= ['\\zserver.cortexlab.net\Data2\Subjects\',animal_name,'\',exp_date,'\photoM'];
%path2Beh= ['\\zserver.cortexlab.net\Data2\Subjects\',animal_name,'\',exp_date,'\',exp_series];
    
path2photoM= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\photoM'];
path2Beh= ['\\zubjects.cortexlab.net\Subjects\',animal_name,'\',exp_date,'\',exp_series];
    

addpath (genpath(path2Beh))
addpath (genpath(path2photoM))


% find neuronal recording info about the session
TargetSessionFound = 0;
isession = 1;
while  TargetSessionFound == 0
  
TargetSessionFound = strcmp(MiceExpInfo.mice(animal_ID).session(isession).Blockname,[exp_date,'_',exp_series,'_',animal_name,'_Block.mat']);
    
    isession = isession + 1;
end

TargetSession = isession - 1;


TrialTimingData= Salvatore_PrepareNeuronBehData(animal_name,exp_date,exp_series);             % construct big behavioural matrix of the session

photoMFileName=MiceExpInfo.mice(animal_ID).session(TargetSession).Neuronfile;

load([path2photoM,'\',photoMFileName(1:end-4),'.mat']);

WaterDig = photoMdata.TTL_1;
TimeStamps=photoMdata.Time_s_;


WaterUp = find(diff(WaterDig)>0.5);


WaterUpTime = TimeStamps(WaterUp);%  / sample_rate;


[TrialTimingData,lag] = Salvatore_AlignBehPhotoM(TrialTimingData,WaterUpTime,sample_rate);   % Water Time alignment
    


% remove correction trials
TrialTimingData (TrialTimingData(:,15)>1,:)=[];
% this is a bit stupid but best we could do. 16 is go cue time and we put
% at as 15. otherwise other codes will be affected.
TrialTimingData (:,15) = TrialTimingData (:,16);
TrialTimingData (:,16)=[];



% save list into the databse
MiceExpInfo.mice(animal_ID).session(TargetSession).TrialTimingData = TrialTimingData;
MiceExpInfo.mice(animal_ID).session(TargetSession).AlignDelay = lag - 0.1;


% MachineName=getComputerName;
% 
% if strcmp(MachineName,'zopamine')
%     
%     cd 'G:\Dropbox\Work\UCL\Science\Analysis Code'
%     
% else
%     cd '/Users/Armin/Dropbox/Work/UCL/Science/Analysis Code'
% end


path2save =which ('Salvatore_savePhotoMcsv2mat');
cd(path2save(1:end-30));

save('MiceExpInfoPhotoM.mat','MiceExpInfo')
save('MiceExpInfoPhotoM_BackUP.mat','MiceExpInfo')


