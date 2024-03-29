function [SessionList] = getSessionList_photoM(animal_name, exp_ID)

if strcmp(animal_name,'ALK068')
    if exp_ID == '7'
        SessionList = [1:13]; % 1,3,5:13
    elseif exp_ID == '23'
        SessionList = [14, 15, 16, 17, 18, 19, 20, 22, 23, 24];
    end
    
elseif strcmp(animal_name, 'ALK070')
    if exp_ID == '7'
        SessionList = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12];
    elseif exp_ID == '23'
        SessionList = [13, 14, 15, 16,17, 18, 19, 20, 21, 22, 23, 24];
    end
    
elseif strcmp(animal_name, 'ALK071')
    if exp_ID == '7'
        SessionList = [1:9];
    elseif exp_ID == '23'
        SessionList = [10, 11, 12,13, 14,19];
    end
    
elseif strcmp(animal_name, 'ALK074')
    if exp_ID == '7'
        SessionList = [1:20];
    elseif exp_ID == '23'
        SessionList = [21,22,23,24,25,26,27];
    end
    
elseif strcmp(animal_name, 'ALK075')
    if exp_ID == '7'
        SessionList = [1:14];
    elseif exp_ID == '23'
        SessionList = [15, 16,17,18,19];
    end
    
elseif strcmp(animal_name, 'ALK078')
    if exp_ID == '7'
        SessionList = [1:22];
    elseif exp_ID == '23'
        SessionList = [23:32];
    end
    
elseif strcmp(animal_name, 'ALK083')
    if exp_ID == '7'
        SessionList = [1:3]; %check this
    elseif exp_ID == '23'
        SessionList = [13:21];
    elseif exp_ID == '38'
        SessionList = [23 24];
    end
    
elseif strcmp(animal_name, 'ALK084')
    if exp_ID == '7'
        SessionList = [1:7];
        
    elseif exp_ID == '23'
        SessionList = [12:29];
    elseif exp_ID == '38'
        SessionList = [33:38];
        
    end
    
elseif strcmp(animal_name, 'MMM001')
    if exp_ID == '7'
        SessionList = [1:12];
    elseif exp_ID == '23'
        SessionList = [13,14,15,16, 18, 19, 20, 21, 22, 23, 24, 25];
    end
    
elseif strcmp(animal_name, 'MMM002')
    if exp_ID == '7'
        SessionList = [1:10];
    elseif exp_ID == '23'
        SessionList = [15,16,17,18,20,21,22,23];
    end
    
elseif strcmp(animal_name, 'MMM003')
    if exp_ID == '7'
        SessionList = [1:3, 5:12];
    elseif exp_ID == '23'
        SessionList = [13:16];
    end
    
elseif strcmp(animal_name, 'MMM005')
    if exp_ID == '7'
        SessionList = [1:7];
    elseif exp_ID == '23'
        SessionList = [10:15, 17, 18, 19,20];
    elseif exp_ID == '38'
        SessionList = [22:27];
        
    end
    
elseif strcmp(animal_name, 'ALK085')
    if exp_ID == '7'
        SessionList = [1:9, 12];
    elseif exp_ID == '23'
        SessionList = [15:23];
        
    elseif exp_ID == '38'
        SessionList = [24:27];
    end
    
elseif strcmp(animal_name, 'MMM006')
    if exp_ID == '7'
        SessionList = [7, 9:11]; % had a look - unusual responses to reward and stim
    elseif exp_ID == '23'
        SessionList = [14:21]; 
        
    end
       
elseif strcmp(animal_name, 'MMM008')
    if exp_ID == '7'
        SessionList = [5, 7]; % 5 good. 6 bad. 7 ok. 8 bad. 
    elseif exp_ID == '23'
        SessionList = []; 
        
    end
    
elseif strcmp(animal_name, 'MMM009')
    if exp_ID == '7'
        SessionList = [1:4, 6:9]; % could not align 5
    elseif exp_ID == '23'
        SessionList = [10:19]; 
    elseif exp_ID =='38'
        SessionList = [30 34 35 36];
    end
    
elseif strcmp(animal_name, 'MMM010')
    if exp_ID == '7'
        SessionList = [1:5];
    elseif exp_ID == '23'
        SessionList = [6:14]; 
    elseif exp_ID == '38'
        SessionList = [20 21];
    end
    
end
end

