function [TrialTimingData_aligned,lag] = Salvatore_AlignBehPhotoM(TrialTimingData,SynchData,sample_rate)

%Align times from Behaviour and photoM datafile (regression of reward analog time and reward time in the beh block file)
% Armin Dec 2017



reward_Times = SynchData'; 

% Think about this timing issue

if length(reward_Times) > length(TrialTimingData(TrialTimingData(:,9)==1,14))
reward_Times(end) = [];                                                             % we ignore the last trial
end


figure; hold on
plot(diff(reward_Times),'r')
plot(0.1 + diff((TrialTimingData(TrialTimingData(:,9)==1,14))),'b')
legend('analog','TrialTimingData')
ylabel('inter event interval')

%a=TrialTimingData(TrialTimingData(:,9)==1,14);
%s=regstats(reward_Times(1:100),a(1:100,:));         % find the time lag between beh file and neuronal file (regress water delivery time)

s=regstats(reward_Times,TrialTimingData(TrialTimingData(:,9)==1,14));         % find the time lag between beh file and neuronal file (regress water delivery time)

lag = s.beta(1)                                                                % neuronal data preceding beh data with this lag


TrialTimingData_aligned(:,[10, 12:14, 16]) = TrialTimingData(:,[10, 12:14, 16]) + lag;

%unacounted_lag = nanmean(reward_Times(1:100)' - a(1:100))

unacounted_lag = nanmean(reward_Times' - TrialTimingData_aligned(TrialTimingData(:,9)==1,14))



% to here

TrialTimingData_aligned(:,[10, 12:14,16]) = TrialTimingData_aligned(:,[10, 12:14,16]); % - unacounted_lag;

TrialTimingData_aligned(:,[1:9 , 11]) = TrialTimingData(:,[1:9 , 11]);
TrialTimingData_aligned(:,15) = TrialTimingData(:,15);




% fig for aligment check
figure
% plot(AnalogTiming(1:5000000),'r')
% hold on

z = zeros(1, floor(reward_Times(end)) * sample_rate);
 
z(floor(reward_Times* sample_rate))= 100;

plot(z(1:5000000) * -20)

zz=TrialTimingData_aligned(TrialTimingData(:,9)==1,14);
z = zeros(1, floor(zz(end) * sample_rate));

z(floor(zz* sample_rate))= -100;

plot(z(1:5000000)* -40 , 'k')
 



