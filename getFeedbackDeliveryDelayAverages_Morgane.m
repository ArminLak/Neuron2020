% get avg feedback delivery delay for each animal and brain region, Exp 23 
% go to line ~35 to change computer =) 

% Morgane September 2019 


clear all 
close all 

% ---- ENTER REQS ---------------------------------------------------------

exp_ID = '23';

% VTA
animal_names = ["ALK068", "ALK070", "ALK071", "ALK084"]; %place comma between animal names 
brain_region = 'VTA';

% NAC
% animal_names = ["ALK078", "MMM001", "MMM002", "MMM005"] ;
% brain_region = 'NAc';

% DMS
% animal_names = ["ALK074", "MMM003", "ALK083", "MMM009", "MMM010"];
%  brain_region = 'DMS';

% -------------------------------------------------------------------------

load('MiceExpInfoPhotoM')

FeedbackDeliveryDelay = zeros(1, numel(animal_names));
totalSessions = 0;

c = 1;
for animal_name = animal_names(c);
    
    [animal_ID, chan_order] = Salvatore_Get_chan_order(animal_name);
%     path2data = ['D:\Raw data for striatal DA paper\',convertStringsToChars(animal_name)]; % morgane's oxford computer
    path2data = ['D\\zubjects.cortexlab.net\Subjects\',convertStringsToChars(animal_name)]; % morgane's london computer 
    addpath(genpath(path2data))
    [SessionList] = getSessionList_photoM(animal_name, exp_ID);
        
    for iSession = SessionList
        
        totalSessions = length(SessionList) + totalSessions;
    
        ParamFileName = MiceExpInfo.mice(animal_ID).session(iSession).Paramname;
        load(ParamFileName);
        
        FeedbackDeliveryDelay(c) = parameters.feedbackDeliveryDelay + FeedbackDeliveryDelay(c);
    
    end
    
    FeedbackDeliveryDelay(c) = FeedbackDeliveryDelay(c)./length(SessionList);

    c = c + 1;
end

PopFeedbackDeliveryDelay = mean(FeedbackDeliveryDelay); 

    
    


