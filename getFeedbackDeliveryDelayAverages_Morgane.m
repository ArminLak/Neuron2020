% get avg feedback delivery delay for each animal and brain region, Exp 23 
% Morgane September 2019 


clear all 
close all 

% ---- ENTER REQS ---------------------------------------------------------

exp_ID = '23';

% VTA
% Animals = [48 50 51 64]
animal_names = ["ALK068"]; %place comma between animal names 
brain_region = 'VTA';

% NAC
% Animals = [56 57 59 66]
% animal_name = ;
% brain_region = 'NAc';

% DMS
%  Animals = [53, 62, 63, 71,72]
% animal_name = ;
%  brain_region = 'DMS';

% -------------------------------------------------------------------------

load('MiceExpInfoPhotoM')

FeedbackDeliveryDelay = zeros(numel(animal_names));

c = 1;
for animal_name = animal_names(c);
    
    [animal_ID, chan_order] = Salvatore_Get_chan_order(animal_name);
    path2data = ['D:\Raw data for striatal DA paper\',convertStringsToChars(animal_name)];
    addpath(genpath(path2data))
    [SessionList] = getSessionList_photoM(animal_name, exp_ID);
        
    for iSession = SessionList
    
        ParamFileName = MiceExpInfo.mice(animal_ID).session(iSession).Paramname;
        load(ParamFileName);
        
        FeedbackDeliveryDelay(c) = parameters.feedbackDeliveryDelay + FeedbackDeliveryDelay;
    
    end
    
    FeedbackDeliveryDelay(c) = FeedbackDeliveryDelay(c)./length(SessionList);

    c = c + 1;
end

PopFeedbackDeliveryDelay = mean(FeedbackDeliveryDelay); 

    
    


