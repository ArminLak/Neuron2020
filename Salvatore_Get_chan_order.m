function [animal_ID, chan_order] =Salvatore_Get_chan_order(animal_name)

if strcmp(animal_name,'EJ008DA')   % 16 chan, 4 * 4 tetrode
    
    chan_order = [1 3 13 15
        2 4 14 16
        5 7 9  11
        6 8 10 12];
    
    animal_ID  = 7;
    
    
elseif strcmp(animal_name,'M150602_MW')   % 16 chan, 4 * 4 tetrode
    
    chan_order = [1 3 13 15
        2 4 14 16
        5 7 9  11
        6 8 10 12];
    
    animal_ID  = 10;
    
elseif strcmp(animal_name,'ALK008')   % 32 chan, 4 * 8 nontetrode, some broken
    
    chan_order = [21 22 31 30 nan nan nan nan
        24 25 28 26 23 27 29 nan
        6  5  3  4   8  2 12 nan
        11 10 13 9  14  0 1 nan];
    
    animal_ID  = 11;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK011')   % 32 chan buzsaki, 4 * 8 nontetrode,
    
    chan_order = [0:7
        8:15
        16:23
        24:31];
    animal_ID  = 14;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK017')   % 32 chan Holtzman, 2 * 16 linear parallel,
    
    chan_order = [0:15
        nan nan nan nan nan 21 22 nan nan 25 26 27 nan 29 30 31];
    animal_ID  = 19;
    
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK028')   % 32 chan Holtzman, 2 * 16 linear parallel,
    
    chan_order = [0:15
        16:31];
    animal_ID  = 26;
    
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK045')   % 64 chan Holtzman, 4 * 16 linear parallel,
    % not finalsied 2017-06-09:
    chan_order = [0:15
        16:31
        32:47
        48: 63];
    animal_ID  = 36;
    
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK052')   % 374 neuropixel,
    chan_order = [];
    animal_ID  = 42;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK068')   % photometry,
    chan_order = [];
    animal_ID  = 48;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK070')   % photometry,
    chan_order = [];
    animal_ID  = 50;
    
    chan_order = chan_order  +1;
elseif strcmp(animal_name,'ALK071')   % photometry,
    chan_order = [];
    animal_ID  = 51;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK074')   % photometry,
    chan_order = [];
    animal_ID  = 53;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK075')   % photometry,
    chan_order = [];
    animal_ID  = 55;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK078')   % photometry,
    chan_order = [];
    animal_ID  = 56;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'MMM001')   % photometry,
    chan_order = [];
    animal_ID  = 57;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'MMM002')   % photometry,
    chan_order = [];
    animal_ID  = 59;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'MMM003')   % photometry,
    chan_order = [];
    animal_ID  = 62;
    
    chan_order = chan_order  +1;
    
    
elseif strcmp(animal_name,'ALK084')   % photometry,
    chan_order = [];
    animal_ID  = 64;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'ALK083')   % photometry,
    chan_order = [];
    animal_ID  = 63;
    
    chan_order = chan_order  +1;
elseif strcmp(animal_name,'MMM004')   % photometry,
    chan_order = [];
    animal_ID  = 65;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name,'MMM005')   % photometry,
    chan_order = [];
    animal_ID  = 66;
    
    chan_order = chan_order  +1;
    
elseif strcmp(animal_name, 'MMM006') % photometry
    chan_order = [];
    animal_ID = 68;
    
    chan_order = chan_order+1;
elseif strcmp(animal_name, 'ALK085') % photometry
    chan_order=[];
    animal_ID = 69;
    
    chan_order = chan_order +1;
end