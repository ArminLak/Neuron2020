function [eventTimes, eventValues,windows]=Salvatore_SelectKernelsPhotoM(ModelArrangment,beep_time,stim_time, action_onsetTime, outcome_time)

%ModelArrangment==8 Stim, Action, Outcome

%ModelArrangment==9  Stim, Action

%ModelArrangment==10 Action

%ModelArrangment==11  Stim, Outcome

%ModelArrangment==12  Action, outcome


%ModelArrangment==13  Stim 

%ModelArrangment==14  outcome 



if ModelArrangment ==8
    
    eventTimes = {beep_time;stim_time; action_onsetTime; outcome_time};
    eventValues = {[];[];[];[]};
    
    windows = {[-400 800];[-200 1600];[-600 200];[200 2000]}; % not in s but sample rate of 1200
    windows = {[-400 800];[-400 2600];[-1000 200];[-400 2000]}; %
    
   % windows = {[-400 800];[-400 2600];[-1000 200];[-400 4000]}; % for plotting only
    
    
elseif ModelArrangment ==9
    
      eventTimes = {stim_time; action_onsetTime};
    eventValues = {[];[]};
    
        windows = {[-400 2600];[-400 2000]}; %


elseif ModelArrangment ==10
    
    eventTimes = {action_onsetTime};
    eventValues = {[]};
    
    windows = {[-1000 200]}; %
    
    
    
% stim,outcome

elseif ModelArrangment ==11
    
    eventTimes = {beep_time;stim_time; outcome_time};
    eventValues = {[];[];[]};
    
    windows = {[-400 800];[-400 2600];[-400 2000]}; %
   windows = {[-400 1200];[200 2600];[0 2000]}; %
 

    
end