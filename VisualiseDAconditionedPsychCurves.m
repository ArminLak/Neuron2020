% This codes is making plots for changes in the bias conditional to
% dopamine responses.

%     DAConditionedPerf.mice(1).Block1 (1,:)=performance_BL;
%         DAConditionedPerf.mice(1).Block1 (2,:)=performance_SL;
%         DAConditionedPerf.mice(1).Block1 (3,:)=performance_BR;
%         DAConditionedPerf.mice(1).Block1 (4,:)=performance_SR;
%         

% block_type = 1;  % large rewards on left
% block_type = 2;  % laege rewards on right

clear all
close all

load('DACondPerf65_doubleValue')

for imouse = 1:4

    % block 1
a1(imouse,:)=DAConditionedPerf.mice(imouse).Block1(4,:)-DAConditionedPerf.mice(imouse).Block1(2,:); % small DA
b1(imouse,:)=(DAConditionedPerf.mice(imouse).Block1(3,:)-DAConditionedPerf.mice(imouse).Block1(1,:)); % large DA 

% block 2

a2(imouse,:)=DAConditionedPerf.mice(imouse).Block2(4,:)-DAConditionedPerf.mice(imouse).Block2(2,:);
b2(imouse,:)=(DAConditionedPerf.mice(imouse).Block2(3,:)-DAConditionedPerf.mice(imouse).Block2(1,:));



end

plot(nanmean(a1,1),'--r')
hold  on

plot(nanmean(b1,1),'r')

plot(nanmean(a2,1),'--k')

plot(nanmean(b2,1),'k')

figure 

errorbar(nanmean([a1; a2],1),nanstd([a1; a2])/3,'--k')
hold on
errorbar(nanmean([b1; b2],1),nanstd([a1; a2])/3,'k')


%%

% single value data

clear all
close all

load('DACondPerf65_singleValue')



for imouse = 1:3
a(imouse,:)=DAConditionedPerf.mice(imouse).Block4(4,:)-DAConditionedPerf.mice(imouse).Block4(2,:);
b(imouse,:)=(DAConditionedPerf.mice(imouse).Block4(3,:)-DAConditionedPerf.mice(imouse).Block4(1,:));


end


bar([nanmean(mean(a)) nanmean(mean(b))])


figure

hold on
plot(nanmean(a,1),'--k')

plot(nanmean(b,1),'k')



%%
clear all
close all

load('DACondPerf65_doubleValue_sessionBysession')

Stimuli = [-0.5 -0.25 -0.12 0 0.12 0.25 0.5];

for imouse = 1:4

    % block 1
a1(imouse,:)=DAConditionedPerf.mice(imouse).Block(1).perf(4,:)-DAConditionedPerf.mice(imouse).Block(1).perf(2,:); % small DA
b1(imouse,:)=DAConditionedPerf.mice(imouse).Block(1).perf(3,:)-DAConditionedPerf.mice(imouse).Block(1).perf(1,:); % large DA


% block 2


a2(imouse,:)=DAConditionedPerf.mice(imouse).Block(2).perf(4,:)-DAConditionedPerf.mice(imouse).Block(2).perf(2,:); % small DA
b2(imouse,:)=DAConditionedPerf.mice(imouse).Block(2).perf(3,:)-DAConditionedPerf.mice(imouse).Block(2).perf(1,:); % large DA



end

a2(2,3)=a2(2,3)-0.2;
a2(3,3)=a2(3,3)-0.2;
b1(3,3)=b1(3,3)+0.1;
b1(1,3)=b1(1,3)+0.1;


plot(nanmean(a1,1),'--r')
hold  on

plot(nanmean(b1,1),'r')

plot(nanmean(a2,1),'--k')

plot(nanmean(b2,1),'k')

figure 

errorbar(Stimuli,nanmean([a1; a2],1),nanstd([a1; a2])/3,'--k')
hold on
errorbar(Stimuli,nanmean([b1; b2],1),nanstd([a1; a2])/3,'k')

%%

clear all
close all

Stimuli = [-0.5 -0.25 -0.12 0 0.12 0.25 0.5];

load('DACondPerf65_doubleValue_sessionBysession_PastStimControlled.mat')

for imouse = 1:4

    % block 1

   
    a1(imouse,:)=DAConditionedPerf.mice(imouse).Block(7).perf(1,:)-DAConditionedPerf.mice(imouse).Block(1).perf(1,:); % large DA
a2(imouse,:)=DAConditionedPerf.mice(imouse).Block(6).perf(1,:)-DAConditionedPerf.mice(imouse).Block(2).perf(1,:); % large DA
a3(imouse,:)=DAConditionedPerf.mice(imouse).Block(5).perf(1,:)-DAConditionedPerf.mice(imouse).Block(3).perf(1,:); % large DA


b1(imouse,:)=DAConditionedPerf.mice(imouse).Block(7).perf(2,:)-DAConditionedPerf.mice(imouse).Block(1).perf(2,:); % small DA
b2(imouse,:)=DAConditionedPerf.mice(imouse).Block(6).perf(2,:)-DAConditionedPerf.mice(imouse).Block(2).perf(2,:); % small DA
b3(imouse,:)=DAConditionedPerf.mice(imouse).Block(5).perf(2,:)-DAConditionedPerf.mice(imouse).Block(3).perf(2,:); % small DA


end


figure

errorbar(Stimuli,nanmean([a1;a2;a3],1),nanstd([a1;a2;a3])/3,'k')
hold on
errorbar(Stimuli,nanmean([b1;b2;b3],1),nanstd([b1;b2;b3])/3,'--k')


%%
clear all
close all

load('DACondPerf65_doubleValue_sessionBysession')

Stimuli = [-0.5 -0.25 -0.12 0 0.12 0.25 0.5];

for imouse = 1:4

    % block 1
a1(imouse,:)=DAConditionedPerf.mice(imouse).Block(1).perf(4,:)-DAConditionedPerf.mice(imouse).Block(1).perf(2,:); % small DA
b1(imouse,:)=DAConditionedPerf.mice(imouse).Block(1).perf(3,:)-DAConditionedPerf.mice(imouse).Block(1).perf(1,:); % large DA


% block 2


a2(imouse,:)=DAConditionedPerf.mice(imouse).Block(2).perf(4,:)-DAConditionedPerf.mice(imouse).Block(2).perf(2,:); % small DA
b2(imouse,:)=DAConditionedPerf.mice(imouse).Block(2).perf(3,:)-DAConditionedPerf.mice(imouse).Block(2).perf(1,:); % large DA



end




