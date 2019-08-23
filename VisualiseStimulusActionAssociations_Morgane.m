close all 
clear all

%DMS
region = 'DMS';
Animals = [53, 62, 63, 71,72]
load('BehPhotoM_Exp23_DMS')
[IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors();
totalChannels = getTotalChanN(region);


PSTH_L = NaN(2,13100);
PATH_L = NaN(2,13100);
PSTH_R = NaN(2,13100);
PATH_R = NaN(2,13100);

c = 1;
animal_count = 0;
for iAnimal = Animals
    
    animal_count = animal_count + 1;
    ChanN = 0;
    
    if ~isempty(BehPhotoM(iAnimal).GrandSummaryL)
        ChanN = ChanN + 1;
    end
    
    if ~isempty(BehPhotoM(iAnimal).GrandSummaryR)
        ChanN = ChanN + 1;
    end
    
    for iChan = 1:ChanN
        
        if isempty(BehPhotoM(iAnimal).GrandSummaryR)
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
            hem = 'l';
        elseif isempty(BehPhotoM(iAnimal).GrandSummaryL)
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
            hem = 'r';
        elseif ~isempty(BehPhotoM(iAnimal).GrandSummaryL) && ~isempty(BehPhotoM(iAnimal).GrandSummaryR)
            if iChan == 1
                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
                hem = 'l';
            elseif iChan==2
                BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
                hem = 'r';
            end
        end
 

               
        
         
        if strcmp(hem, 'l')
            tempStimData = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeCorrect;
            tempStimData = tempStimData([1 7],:);
            PSTH_L = PSTH_L + tempStimData;
            
            tempActionData = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeCorrect;
            tempActionData = tempActionData([1 7],:);
            PATH_L = PATH_L + tempActionData;
            
        elseif strcmp(hem, 'r')
            tempStimData = BehPhotoM(iAnimal).GrandSummary.StimRasterLargeCorrect;
            tempStimData = tempStimData([1 7],:);
            PSTH_R = PSTH_R + tempStimData;
            
            tempActionData = BehPhotoM(iAnimal).GrandSummary.ActionRasterLargeCorrect;
            tempActionData = tempActionData([1 7],:);
            PATH_R = PATH_R + tempActionData;
            
        end

        c = c+1;
    end
    
end

% to do: 
% need to separate by ContraContra, ContraIpsi, IpsiContra, IpsiIpsi

PSTH_Contra = PSTH_L(2,:) + PSTH_R(1,:);
PSTH_Ipsi   = PSTH_L(1,:) + PSTH_R(2,:);
PATH_Contra = PATH_L(2,:) + PATH_R(1,:);
PATH_Ipsi   = PATH_L(1,:) + PATH_R(2,:);

PSTH_Contra = PSTH_Contra ./ totalChannels;
PSTH_Ipsi = PSTH_Ipsi ./ totalChannels;
PATH_Contra = PATH_Contra ./ totalChannels;
PATH_Ipsi = PATH_Ipsi ./ totalChannels;
% --- plotting 





function [IpsiContraColor, ErrorCorrectColor, SmallLargeColor] = getColors()
IpsiContraColor = [ 0 0 0
                    153/255 51/255 1
                    ];
                
ErrorCorrectColor = [   'r'
                        'g'
                        ];
                    
SmallLargeColor = [ 0.15 0.75 0
                    0.15 0.4 0
                    ];

end

function totalChannels = getTotalChanN(region)

if strcmpi(region, 'DMS')
    totalChannels = 7;
elseif strcmpi(region, 'NAC')
    totalChannels = 5;
elseif strcmpi(region, 'VTA')
    totalChannels = 4;
end
end

