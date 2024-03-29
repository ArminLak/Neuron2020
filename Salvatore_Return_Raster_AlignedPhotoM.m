function [Raster_Matrix]=Salvatore_Return_Raster_AlignedPhotoM(TimeStamps,event_times,DeltaFoverF,start,stop,downsampleScale)

% generate raster for photoM
% Armin 2017-12-21, London


%PSTH_Matrix=zeros(length(event_times),abs(stop - start)* 1000);

trace = 1;
%[y sort_index]=sort(sort_event_times-event_times);
DeltaFoverF = DeltaFoverF(:);
TimeStamps = TimeStamps(:);

DeltaFoverF = smooth(DeltaFoverF,downsampleScale); 
DeltaFoverF = downsample(DeltaFoverF,10);  % downsampling

TimeStamps = downsample(TimeStamps,downsampleScale);

limit = (stop - start ) * (12000/downsampleScale)- 100;   % 1200 comes because samplingRate =12000 and downsampling is 10
for ievent=event_times' %(sort_index)'

    index  = find(TimeStamps >  (ievent + start) & TimeStamps <  (ievent + stop));

    if ~isempty(index)

    if length(index) < limit
        DeltaFoverF=[DeltaFoverF;zeros(limit-length(index),1)];
        Raster_Matrix(trace,1:length(index)) = DeltaFoverF(index(1:length(index)),:) ;
    
    else
        
        Raster_Matrix(trace,:) = DeltaFoverF(index(1:limit),:) ;
        
    end
    
end
 
    trace = trace +1;
    
end


%Raster_Matrix = logical(PSTH_Matrix);

%PSTH_Matrix = PSTH_Matrix * 1000;


end





