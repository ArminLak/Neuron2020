% This code is for visualisation of outcome responses in exp 38
% Morgane November 2019

close all 
clearvars -except BehPhotoM

%DMS
Animals = [53, 62, 63, 71,72];
% load('BehPhotoM_Exp38_DMS')


Stimz2plot = [1 7]; % 2 and 6 are 0.25 contrast; 1 and 7 are 0.5; 

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 200; 

TimingVisualise = [-0.2 0.8
                   -0.7, 0.7
                   -0.2, 0.8]; % stim, action, reward in s
               
StartTime = 3700; % saved in the database.


%% processing 
[IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors();

GrandPopRwdLargeError     = zeros(2,13100);
GrandPopRwdLargeCorrect   = zeros(2,13100);
GrandPopRwdSmallCorrect   = zeros(2,13100);



c = 1;
animal_count = 0;
chan_count = 0;

for iAnimal = Animals
    
    animal_count = animal_count + 1;
    
    if isempty(BehPhotoM(iAnimal).GrandSummaryR)        % uni hem animals L
        ChanN = 1;
    elseif isempty(BehPhotoM(iAnimal).GrandSummaryL)    % uni hem animals R
        ChanN = 2;  
    else                                                % bi hem animals
        ChanN = [1 2];
    end
 
    for iChan = ChanN
    chan_count = chan_count + 1;
        BehPhotoM(iAnimal).GrandSummary     = [];
        SingleAnimalRwdTraceLargeCorrect   = [];
        SingleAnimalRwdTraceSmallCorrect   = [];
        SingleAnimalRwdTraceLargeError     = [];

        if iChan == 1
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryL;
        elseif iChan == 2
            BehPhotoM(iAnimal).GrandSummary = BehPhotoM(iAnimal).GrandSummaryR;
        end

        
        
        SingleAnimalRwdTraceLargeCorrect       = BehPhotoM(iAnimal).GrandSummary.Rew2DACorr;
        SingleAnimalRwdTraceSmallCorrect       = BehPhotoM(iAnimal).GrandSummary.RewAwayDACorr;
        SingleAnimalRwdTraceLargeError         = BehPhotoM(iAnimal).GrandSummary.Rew2DAErr;
        
                SingleAnimalNormStimDenom = max( [max(SingleAnimalRwdTraceLargeCorrect), max(SingleAnimalRwdTraceSmallCorrect), max(SingleAnimalRwdTraceLargeError)] );
        SingleAnimalRwdTraceLargeCorrect         = SingleAnimalRwdTraceLargeCorrect ./ SingleAnimalNormStimDenom;                
        SingleAnimalRwdTraceSmallCorrect       = SingleAnimalRwdTraceSmallCorrect ./ SingleAnimalNormStimDenom;
        SingleAnimalRwdTraceLargeError       = SingleAnimalRwdTraceLargeError ./ SingleAnimalNormStimDenom;

        
        if iChan == 1
            GrandPopRwdLargeCorrect(:,:,chan_count)    = SingleAnimalRwdTraceLargeCorrect;
            GrandPopRwdSmallCorrect(:,:,chan_count)    = SingleAnimalRwdTraceSmallCorrect;
            GrandPopRwdLargeError(:,:,chan_count)      = SingleAnimalRwdTraceLargeError;

        elseif iChan == 2
            GrandPopRwdLargeCorrect(:,:,chan_count)    = flipud(SingleAnimalRwdTraceLargeCorrect);
            GrandPopRwdSmallCorrect(:,:,chan_count)    = flipud(SingleAnimalRwdTraceSmallCorrect);
            GrandPopRwdLargeError(:,:,chan_count)      = flipud(SingleAnimalRwdTraceLargeError);       
               
        end
    end   
end

AvgGrandPopRwdLargeCorrect    = mean(GrandPopRwdLargeCorrect, 3);
AvgGrandPopRwdSmallCorrect    = mean(GrandPopRwdSmallCorrect, 3);
AvgGrandPopRwdLargeError      = mean(GrandPopRwdLargeError, 3);



%% plots 




%% functions =)
function [IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors()
IpsiContraColor = [ 0 0 0
                    0.2 0.2 0.2
                    0.4 0.4 0.4
                    0.6 0.6 0.7
                    0.6 0.4 0.8
                    0.6 0.2 0.9
                    0.6 0 1
                    ];
                
                IpsiContraColor2 = [ 0 0 0
                    0.6 0 1
                    ];
                
ErrorCorrectColor = [   'r'
                        'g'
                        ];
                    
SmallLargeColor = [ 0.15 0.75 0
                    0.15 0.4 0
                    ];

end

