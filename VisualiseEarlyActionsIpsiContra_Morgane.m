% Morgane October 2019

close all 
clear all

%DMS
% Animals = [53, 62, 63, 71,72];
Animals = [53];
load('BehPhotoM_Exp23_DMS')
load('MiceExpInfoPhotoM.mat')

%NAC
% Animals = [56 57 59 66];
% load('BehPhotoM_Exp23_NAc')
% load('MiceExpInfoPhotoM.mat')



Stimz2plot = 1:7; % or e.g. [1 7] (max level ipsi and contra)

min_RT = 0.9;
max_RT = 2;

sampleRate = 1200;
StartTime = 3700;
smooth_factor = 200;

TimingVisualise = [-0.2 0.8
    -0.7, 0.7
    -0.2, 0.8]; % stim, action, reward in s

StartTime = 3700; % saved in the database.


% --- get and organise data -----------------------------------------------
[IpsiContraColor, IpsiContraColor2, ErrorCorrectColor, SmallLargeColor] = getColors();

earlyBehDataLeft = [];
earlyBehDataRight = [];
earlyActionRasterLeft = [];
earlyActionRasterRight = [];
lateBehDataLeft = [];
lateBehDataRight = [];
lateActionRasterLeft = [];
lateActionRasterRight = [];

c = 1;
animal_count = 0;
chan_count = 0;
for iAnimal = Animals
    animal_name = MiceExpInfo.mice(iAnimal).name;
    [SessionList] = getSessionList_photoM(animal_name, '23');
    
    animal_count = animal_count + 1;
    
    if isempty(BehPhotoM(iAnimal).GrandSummaryR)        % uni hem animals L
        ChanN = 1;
    elseif isempty(BehPhotoM(iAnimal).GrandSummaryL)    % uni hem animals R
        ChanN = 2;
    else                                                % bi hem animals
        ChanN = [1 2];
    end
    
    for iSession = 1:length(BehPhotoM(iAnimal).Session)
        
        tempBehData = [];
        tempBehData = MiceExpInfo.mice(iAnimal).session(SessionList(iSession)).TrialTimingData;
        earlyActionTrials = (tempBehData(:,15)>tempBehData(:,10)) & (tempBehData(:,10)-tempBehData(:,13)>min_RT) & (tempBehData(:,10)-tempBehData(:,13)<max_RT);
        lateActionTrials = (tempBehData(:,15)<tempBehData(:,10)) & (tempBehData(:,10)-tempBehData(:,13)<max_RT);
        largeCorrectTrials = (tempBehData(:,9)==1 & ( (tempBehData(:,8)==1 & tempBehData(:,3)==-1 ) | ( tempBehData(:,8)==2 & tempBehData(:,3)==1 )));
        
        for iChan = ChanN
            
            if iChan == 1
                
                tempActionData = [];
                tempActionData = BehPhotoM(iAnimal).Session(iSession).NeuronActionL(and(earlyActionTrials,largeCorrectTrials),:);
                tempActionData = tempActionData ./ max(max(tempActionData));
                earlyActionRasterLeft = [earlyActionRasterLeft; tempActionData];
                earlyBehDataLeft = [earlyBehDataLeft; tempBehData(and(earlyActionTrials,largeCorrectTrials), :)];
                
                tempActionData = [];
                tempActionData = BehPhotoM(iAnimal).Session(iSession).NeuronActionL(and(lateActionTrials,largeCorrectTrials),:);
                tempActionData = tempActionData ./ max(max(tempActionData));
                lateActionRasterLeft = [lateActionRasterLeft; tempActionData];
                lateBehDataLeft = [lateBehDataLeft; tempBehData(and(lateActionTrials,largeCorrectTrials), :)];                
                
                
            elseif iChan == 2
                
                tempActionData = [];
                tempActionData = BehPhotoM(iAnimal).Session(iSession).NeuronActionR(and(earlyActionTrials,largeCorrectTrials),:);
                tempActionData = tempActionData ./ max(max(tempActionData));
                earlyActionRasterRight = [earlyActionRasterRight; tempActionData];
                earlyBehDataRight = [earlyBehDataRight; tempBehData(and(earlyActionTrials,largeCorrectTrials), :)];
                
                tempActionData = [];
                tempActionData = BehPhotoM(iAnimal).Session(iSession).NeuronActionR(and(lateActionTrials,largeCorrectTrials),:);
                tempActionData = tempActionData ./ max(max(tempActionData));
                lateActionRasterRight = [lateActionRasterRight; tempActionData];
                lateBehDataRight = [lateBehDataRight; tempBehData(and(lateActionTrials,largeCorrectTrials), :)];
                
            end
            
        end
            
    end
    
end
%%
earlyBehDataIpsi = [];
earlyBehDataContra = [];
lateBehDataIpsi = [];
lateBehDataContra = [];

earlyActionRasterIpsi = [];
earlyActionRasterContra = [];
earlyActionRasterZero = [];
lateActionRasterIpsi = [];
lateActionRasterContra = [];
lateActionRasterZero = [];

if contains(num2str(ChanN),'1')
    earlyBehDataIpsi = [earlyBehDataIpsi; earlyBehDataLeft(earlyBehDataLeft(:,2)<0,:)];
    earlyBehDataContra = [earlyBehDataContra; earlyBehDataLeft(earlyBehDataLeft(:,2)>0,:)];
    lateBehDataIpsi = [lateBehDataIpsi; lateBehDataLeft(lateBehDataLeft(:,2)<0,:)];
    lateBehDataContra = [lateBehDataContra; lateBehDataLeft(lateBehDataLeft(:,2)>0,:)];
    
    earlyActionRasterIpsi = [earlyActionRasterIpsi; earlyActionRasterLeft(earlyBehDataLeft(:,2)<0,:)];
    earlyActionRasterContra = [earlyActionRasterContra; earlyActionRasterLeft(earlyBehDataLeft(:,2)>0,:)];
    earlyActionRasterZero = [earlyActionRasterZero; earlyActionRasterLeft(earlyBehDataLeft(:,2)==0,:)];
    lateActionRasterIpsi = [lateActionRasterIpsi; lateActionRasterLeft(lateBehDataLeft(:,2)<0,:)];
    lateActionRasterContra = [lateActionRasterContra; lateActionRasterLeft(lateBehDataLeft(:,2)>0,:)];
    lateActionRasterZero = [lateActionRasterZero; lateActionRasterLeft(lateBehDataLeft(:,2)==0,:)];
end

if contains(num2str(ChanN),'2')
    earlyBehDataIpsi = [earlyBehDataIpsi; earlyBehDataRight(earlyBehDataRight(:,2)>0, :)];
    earlyBehDataContra = [earlyBehDataContra; earlyBehDataRight(earlyBehDataRight(:,2)<0, :)];
    lateBehDataIpsi = [lateBehDataIpsi; lateBehDataRight(lateBehDataRight(:,2)>0, :)];
    lateBehDataContra = [lateBehDataContra; lateBehDataRight(lateBehDataRight(:,2)<0, :)];
    
    earlyActionRasterIpsi = [earlyActionRasterIpsi; earlyActionRasterRight(earlyBehDataRight(:,2) >0,:)];
    earlyActionRasterContra = [earlyActionRasterContra; earlyActionRasterRight(earlyBehDataRight(:,2) <0,:)];
    earlyActionRasterZero = [earlyActionRasterZero; earlyActionRasterRight(earlyBehDataRight(:,2)==0,:)];
    lateActionRasterIpsi = [lateActionRasterIpsi; lateActionRasterRight(lateBehDataRight(:,2) >0,:)];
    lateActionRasterContra = [lateActionRasterContra; lateActionRasterRight(lateBehDataRight(:,2) <0,:)];
    lateActionRasterZero = [lateActionRasterZero; lateActionRasterRight(lateBehDataRight(:,2)==0,:)];
end


figure; 

subplot(2,2,1)
c = 1;
for iStim = [0.5 0.25 0.12]
plot( smooth(mean(earlyActionRasterIpsi(abs(earlyBehDataIpsi(:,2))==iStim,:)),smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
hold on;
c = c+1;
end

plot( smooth( mean(earlyActionRasterZero), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
c = c+1;

for iStim = [0.12 0.25 0.5];
plot( smooth ( mean ( earlyActionRasterContra(abs(earlyBehDataContra(:,2))==iStim,:)), smooth_factor), 'color', IpsiContraColor(5,:), 'LineWidth', 2)
c = c+1;
end

set(gca, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime + (TimingVisualise(2,1)*sampleRate), StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTickLabel', {'-0.7','0','0.7'},'TickDir','out','Box','off')


    
    
    
subplot(2,2,2)
c = 1;
for iStim = [0.5 0.25 0.12]
plot( smooth(mean(lateActionRasterIpsi(abs(lateBehDataIpsi(:,2))==iStim,:)),smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
hold on;
c = c+1;
end

plot( smooth( mean(lateActionRasterZero), smooth_factor), 'color', IpsiContraColor(c,:), 'LineWidth', 2)
c = c+1;

for iStim = [0.12 0.25 0.5];
plot( smooth ( mean ( lateActionRasterContra(abs(lateBehDataContra(:,2))==iStim,:)), smooth_factor), 'color', IpsiContraColor(5,:), 'LineWidth', 2);
c = c+1;
end

set(gca, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime + (TimingVisualise(2,1)*sampleRate), StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTickLabel', {'-0.7','0','0.7'},'TickDir','out','Box','off')
    
    subplot(2,2,3)

plot( smooth(mean(earlyActionRasterIpsi),smooth_factor), 'color', IpsiContraColor2(1,:), 'LineWidth', 2)
hold on; 
% plot( smooth( mean(earlyActionRasterZero), smooth_factor), 'color', IpsiContraColor(4,:), 'LineWidth', 2)
plot( smooth ( mean ( earlyActionRasterContra), smooth_factor), 'color', IpsiContraColor2(2,:), 'LineWidth', 2)


set(gca, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime + (TimingVisualise(2,1)*sampleRate), StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTickLabel', {'-0.7','0','0.7'},'TickDir','out','Box','off')
    
    
        subplot(2,2,4)

plot( smooth(mean(lateActionRasterIpsi),smooth_factor), 'color', IpsiContraColor2(1,:), 'LineWidth', 2)
hold on; 
% plot( smooth( mean(earlyActionRasterZero), smooth_factor), 'color', IpsiContraColor(4,:), 'LineWidth', 2)
plot( smooth ( mean ( lateActionRasterContra), smooth_factor), 'color', IpsiContraColor2(2,:), 'LineWidth', 2)
plot( smooth ( mean ( earlyActionRasterContra), smooth_factor), 'color', IpsiContraColor(5,:), 'LineWidth', 2)


set(gca, 'ylim', [-0.4 1], 'xlim', [StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime + (TimingVisualise(2,1)*sampleRate), StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTickLabel', {'-0.7','0','0.7'},'TickDir','out','Box','off')


    figure; 
    
    stimz = [0.12 0.25 0.5];
    for i = 1:length(stimz)
byContraStimActionResponses(i,:) = mean(earlyActionRasterContra(abs(earlyBehDataContra(:,2))==stimz(i),:), 1);
byIpsiStimActionResponses(i,:) = mean(earlyActionRasterIpsi(abs(earlyBehDataIpsi(:,2))==stimz(i),:), 1);
    end

plot( smooth(mean(byContraStimActionResponses),smooth_factor), 'color', IpsiContraColor2(2,:), 'LineWidth', 2)
hold on;
plot( smooth(mean(byIpsiStimActionResponses),smooth_factor), 'color', IpsiContraColor2(1,:), 'LineWidth', 2)
    
    
set(gca, 'ylim', [0 0.5], 'xlim', [StartTime + (TimingVisualise(2,1)*sampleRate) StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'YTick', [0 0.5 1], 'XTick', [StartTime + (TimingVisualise(2,1)*sampleRate), StartTime,  StartTime + (TimingVisualise(2,2)*sampleRate)], ...
        'XTickLabel', {'-0.7','0','0.7'},'TickDir','out','Box','off')
    

%% functions
 
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
    
    function totalChannels = getTotalChanN(region)
    
    if strcmpi(region, 'DMS')
        totalChannels = 7;
    elseif strcmpi(region, 'NAC')
        totalChannels = 5;
    elseif strcmpi(region, 'VTA')
        totalChannels = 4;
    end
    end

