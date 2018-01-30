% sandbox for looking at preliminary photometry Data

clear all

samplingRate = 12000;

PhotoMdata=readtable('\\zserver.cortexlab.net\Data2\Subjects\ALK068\2017-12-22\photoM\ALK068_2017-12-22_1_F.csv');

%data=csvread('\\zserver.cortexlab.net\Data2\Subjects\ALK068\2017-12-18\photoM\ALK068_2017-12-18_5.csv',1);

TimeStamps=PhotoMdata.Time_s_;
DeltaFoverF = PhotoMdata.AnalogIn_2_dF_F0;
WaterUp = PhotoMdata.TTL_1;


%figure; subplot(3,1,1);plot(data(:,1),data(:,4)); hold on

subplot(3,1,1);plot(TimeStamps, DeltaFoverF)
subplot(3,1,2);plot(TimeStamps, WaterUp)


a = diff(WaterUp);
index=find(a>0.5);





c = 1;
for i=index'
response(c,:) =DeltaFoverF(i-12000:i+24000);
c=c+1;
end

subplot(3,1,3);plot([-12000:24000] / samplingRate,mean(response))