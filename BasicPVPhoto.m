%basic photometry session analysis

%load a bpod session, ie struct SessionData
basepath = 'C:\Data\DataPostdoc\PV-Photometry';
animal = 'tp30';
session = 'TP30_NosePoke_Feb15_2019_Session1.mat';
OUTNAME = 'Photometry3.pdf';

%params
channel = 1; %1=green, 2=red
RealignLeaving=false;
SortByWT = true;%for leaaving plot
SortByReward=true;%for reward plot
minWT = 2; %look only at waiting time higher than that (for all plots)
minRewardDelay = 1; %for average plot, consider only trials with minimum reward delay

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
file = fullfile(basepath,animal,'NosePoke','Session Data',session);
load(file)

nTrials = SessionData.nTrials;

TaskParameters = SessionData.Settings;

if channel == 1
    chname = 'Green';
    PhotoData = SessionData.NidaqData;%green chnannel
    modFreq = TaskParameters.GUI.LED1_Freq;
    modAmp = TaskParameters.GUI.LED1_Amp;
    
elseif channel ==2
    chname = 'Red';
    PhotoData = SessionData.Nidaq2Data;%red chnannel
    modFreq = TaskParameters.GUI.LED2_Freq;
    modAmp = TaskParameters.GUI.LED2_Amp;
end

%assemble de-modulate df/f data for each trial, aligned to an event
DemodPhotoDataReward=cell(1,nTrials);
DemodPhotoData0=cell(1,nTrials);
RewardStart = nan(1,nTrials);
ResponseStart = nan(1,nTrials);
ResolutionTime= nan(1,nTrials);
Baseline= nan(1,nTrials);
% BaselineCho= nan(1,nTrials);
for n = 1:nTrials
    
    data = PhotoData{n};
    statetimes = SessionData.RawEvents.Trial{n}.States;
    % rewarded
    if SessionData.Custom.Rewarded(n)
        if   SessionData.Custom.ChoiceLeft(n)==1
            name = 'water_L';
        elseif SessionData.Custom.ChoiceLeft(n) == 0
            name = 'water_R';
        else
            DemodPhotoDataReward{n}=nan(150,3);
        end
        
        timetozero = statetimes.(name)(1);
        RewardStart(n) = timetozero;
        ResponseStart(n) = statetimes.wait_Sin(2);
        Mod = Nidaq_modulation(modAmp,modFreq,TaskParameters);
        
        [currentNidaq1, rawNidaq1]=NidaqDemod(data(:,1),Mod,modFreq,modAmp,timetozero,TaskParameters);
        
        DemodPhotoDataReward{n} = currentNidaq1;
        
    end
    
    %leaving

        if RealignLeaving
            timetozero = 0;
        else
            timetozero=statetimes.ITI(1)-TaskParameters.GUI.GracePeriod;
        end
        ResolutionTime(n)=timetozero;
        
        Mod = Nidaq_modulation(modAmp,modFreq,TaskParameters);
        
        [currentNidaq1, rawNidaq1,base]=NidaqDemod(data(:,1),Mod,modFreq,modAmp,timetozero,TaskParameters);
        
        DemodPhotoData0{n} = currentNidaq1;
%         Baseline(n) = base;
        Baseline(n)= quantile(DemodPhotoData0{n}(1:20,2),0.25);
 
    
end

%renormalize?
for n =1:nTrials
    if ~isempty(DemodPhotoData0{n})
        DemodPhotoData0{n}(:,3)= (DemodPhotoData0{n}(:,2) - Baseline(n)) / Baseline(n) * 100;
    end
    if ~isempty(DemodPhotoDataReward{n})
     DemodPhotoDataReward{n}(:,3)= (DemodPhotoDataReward{n}(:,2) - Baseline(n)) / Baseline(n) * 100;
    end
end

%% determine 'leaving decision'
if RealignLeaving
 leaveT = 1; %minimum leaving time require to count as 'leave decision'
 
wt=nan(1,nTrials);
wt_start=nan(1,nTrials);
time_in_next_trial=nan(1,nTrials);
for n = 1:nTrials-1
    
     if ~SessionData.Custom.Rewarded(n) && ~isnan(SessionData.Custom.ChoiceLeft(n))
         if SessionData.Custom.ChoiceLeft(n)
             ChoicePortOut = 'Port1Out';
             ChoicePortIn = 'Port1In';
         else
             ChoicePortOut = 'Port3Out';
             ChoicePortIn = 'Port3In';
         end
         
         events = SessionData.RawEvents.Trial{n}.Events;
         
         in =  events.(ChoicePortIn)(events.(ChoicePortIn) > events.Port2In(1));
         if isfield(events,ChoicePortOut)
         out =  events.(ChoicePortOut)(events.(ChoicePortOut) > events.Port2In(1));
         else
             out=[];
         end
         

         if length(out)>1
         diff = in(2:length(out)) - out(1:end-1);
         %determine the first leave that was 'long' (potentially longer then
         %the grace period)
         ii = find(diff > leaveT,1,'first');
         else
             ii=[]; %only one leaving time --> have to check next trial
         end
        
         if length(in) == length(out) && isempty(ii) %regular case
         %tricky cases are trials where the last event is a poke in (given ITI>leaveT, which can be assumed,
         %and neglecting the case where a trial is ending *during* the
         %generaly short grace period)
             ii = length(out);
         end
         
         if ~isempty(ii)
             finalout = out(ii);
         else
             
             %check next trial
             %careful, there might be missing events in between trial
             %updating - we assume we don't miss too much
             nextevents = SessionData.RawEvents.Trial{n+1}.Events;
             if isfield(nextevents,ChoicePortIn)
                 nextin = nextevents.(ChoicePortIn);
                 if isfield(nextevents,'Port2In')
                    nextin = nextevents.(ChoicePortIn)(nextevents.(ChoicePortIn) < nextevents.Port2In(1));
                 end
                 nextin = nextin - SessionData.TrialStartTimestamp(n) +  SessionData.TrialStartTimestamp(n+1);
             else
                 nextin=[];
             end
             if isfield(nextevents,ChoicePortOut)
                 nextout =  nextevents.(ChoicePortOut);
                 if isfield(nextevents,'Port2In')
                nextout =  nextevents.(ChoicePortOut)(nextevents.(ChoicePortOut) < nextevents.Port2In(1));
                 end
                 nextout = nextout - SessionData.TrialStartTimestamp(n) +  SessionData.TrialStartTimestamp(n+1);
             else
                 nextout=[];
             end
             
             if length(nextout)>length(nextin)
                 nextin = [NaN, nextin];
             end
             
             if length(nextout)>1
             diff = nextin(2:end) - nextout(1:end-1);
              ii = find(diff > leaveT,1,'first');
              if isempty(ii)
                  ii=length(nextout);
              end
             else
                 ii=length(nextout);
             end
             
             if ii == 0
                 if ~isempty(out)
                 finalout = out(end); %default, leave it at old trial
                 else
                     finalout=NaN;
                 end
             else
                finalout = nextout(ii);
                time_in_next_trial(n) = finalout + SessionData.TrialStartTimestamp(n) -  SessionData.TrialStartTimestamp(n+1);
             end
             
         end
         
         %final determing waiting time
         wt(n) = finalout-in(1);
         wt_start(n) = in(1);
         
     end
end

%% align and patch together photometry data to 'leaving decision'
DemodPhotoData2=cell(1,nTrials);
for n = 1:nTrials-1
    if ~isnan(wt(n))
         tleave = wt_start(n) + wt(n); %default case
        if ~isnan( time_in_next_trial(n))
            %here we have to stitch together
            tnext = time_in_next_trial(n);
           
             DemodPhotoData2(n) = DemodPhotoData0(n+1);
             %exclude first data point due to acquisition start
              DemodPhotoData2{n}(1,2)=NaN;
              DemodPhotoData2{n}(1,3)=NaN;
             DemodPhotoData2{n}(:,1) = DemodPhotoData2{n}(:,1)-tnext;
             
             %re-baseline
             DemodPhotoData2{n}(:,3) = (DemodPhotoData2{n}(:,2)-Baseline(n))/Baseline(n);
             
             
        else
            DemodPhotoData2(n) = DemodPhotoData0(n);
             DemodPhotoData2{n}(:,1) = DemodPhotoData2{n}(:,1)-tleave;
        end
    end
end

else %if not realign
    DemodPhotoData2=DemodPhotoData0;
    ii= ~SessionData.Custom.Rewarded & ~isnan(SessionData.Custom.ChoiceLeft);
    DemodPhotoData2(~ii)=repmat({[]},1,sum(~ii));
    wt=nan(1,nTrials);
    for n =1:nTrials
        if ii(n)
            if SessionData.Custom.ChoiceLeft(n)
                FeedbackPortTimes = SessionData.RawEvents.Trial{n}.States.wait_L;
                wt(n) = FeedbackPortTimes(end,end)-FeedbackPortTimes(1,1);
            else
                FeedbackPortTimes = SessionData.RawEvents.Trial{n}.States.wait_R;
                wt(n) = FeedbackPortTimes(end,end)-FeedbackPortTimes(1,1);
            end
        end
    end
end

LeaveTrialIndex = 1:nTrials;
RewardTrialIndex=1:nTrials;
RewardDelay = SessionData.Custom.RewardDelay(1:nTrials);


%% remove empty trials
delI = cellfun(@isempty,DemodPhotoDataReward);
DemodPhotoDataReward(delI)=[];
RewardTrialIndex(delI)=[];

delI = cellfun(@isempty,DemodPhotoData2);
DemodPhotoData2(delI)=[];
LeaveTrialIndex(delI)=[];


%% prepare for plotting
%reward
minT = -6;
maxT = 2;
tt = (maxT-minT)*10;
time_reward = linspace(minT,maxT,tt);
PlotDataReward = nan(length(DemodPhotoDataReward),tt);
for n = 1:length(DemodPhotoDataReward)
    data = DemodPhotoDataReward{n};
    ii = find(data(:,1)>0,1,'first')-1;
    
    if ~isempty(ii) && ii>1
        
        mini = max([1,ii+minT*10]);
        maxi = min([size(data,1),ii+maxT*10]);
        
        lower = -minT*10 - (ii - mini)  + 1;
       upper =  (maxi-ii) - minT*10  ;
        
        PlotDataReward(n,lower:upper) = data(mini+1:maxi,3)';
    end
end

%leaving
minT = -6;
maxT = 2;
tt = (maxT-minT)*10;
time_leave = linspace(minT,maxT,tt);
PlotDataLeave = nan(length(DemodPhotoData2),tt);
for n = 1:length(DemodPhotoData2)
    data = DemodPhotoData2{n};
    ii = find(data(:,1)>0,1,'first')-1;
    
    if ~isempty(ii) && ii>1
        
        mini = max([1,ii+minT*10]);
        maxi = min([size(data,1),ii+maxT*10]);
        
        lower = -minT*10 - (ii - mini)  + 1;
        upper =  (maxi-ii) - minT*10  ;
        
        PlotDataLeave(n,lower:upper) = data(mini+1:maxi,3)';
    end
end

%delete trials with many NaNs or other critera
delI = sum(isnan(PlotDataReward),2)>30;
PlotDataReward(delI,:)=[];
RewardTrialIndex(delI)=[];

delI = sum(isnan(PlotDataLeave),2)>30 | wt(LeaveTrialIndex)'<minWT;
PlotDataLeave(delI,:)=[];
LeaveTrialIndex(delI)=[];

%normalize per trial 
% baselineidx = 1:40;
% PlotDataReward = PlotDataReward - repmat( nanmean(PlotDataReward(:,baselineidx),2),1,size(PlotDataReward,2));
% PlotDataLeave = PlotDataLeave - repmat( nanmean(PlotDataLeave(:,baselineidx),2),1,size(PlotDataLeave,2));

time = linspace(minT,maxT,tt);

%resort data?
if SortByReward
    [rewardsort,sortireward]=sort(RewardDelay(RewardTrialIndex));
    PlotDataReward=PlotDataReward(sortireward,:);
end
if SortByWT
    [wtsort,sortileave]=sort(wt(LeaveTrialIndex));
    PlotDataLeave=PlotDataLeave(sortileave,:);
end

%% plot

figure('Color',[1,1,1],'Position',[      668   336   705   515])

%reward
ax1=subplot(3,2,[1,3]); hold on

t1 = min(time_reward);
t0 = 0;
t2 = max(time_reward);
i_t1 = 1;
i_t0 = find(time_reward>0,1,'first');
i_t2=length(time_reward);

imagesc(PlotDataReward)
set(gca,'XTick',[i_t1,i_t0,i_t2],'XTickLabel',[t1,t0,t2],'YLim',[0,size(PlotDataReward,1)],'XLim',[1,size(PlotDataReward,2)])
ax1.YAxis.Direction='reverse';
cc = colorbar();
ylabel(cc,'dF/F (%)')
xlabel('Time from reward (s)')
ylabel('Trials')
uicontrol('Style','text','String',strrep(session(1:end-4),'_','-'),'FontName','Arial','Position',[10,580,200,10],'BackgroundColor',[1,1,1])
% title(chname)

%ticks reward
axes(ax1)
PlotDataRewardCorr=PlotDataReward;
responsealigned = ResponseStart(RewardTrialIndex) - RewardStart(RewardTrialIndex);
responsealigned= responsealigned(sortireward);
responsealigned_xidx = round((responsealigned+6) * 10)+1;
responsealigned_xidx(responsealigned_xidx<=0)=NaN;
for n =1:length(responsealigned_xidx)
    line([responsealigned_xidx(n),responsealigned_xidx(n)],[n-1,n],'Color','k','LineWidth',0.25)
    PlotDataRewardCorr(n,1:min([responsealigned_xidx(n),size(PlotDataRewardCorr,2)]))=NaN;
end


ax2=subplot(3,2,5);
y=nanmean(PlotDataRewardCorr);
plot(time_reward,y,'-r')
set(gca,'XTick',[t1,t0,t2]);
set(gca,'YLim',[min(y)-0.01,max(y)+0.01])
ylabel('dF/F (%)')
xlabel('Time from reward (s)')
pos1 = get(ax1,'Position');
pos=get(gca,'Position');
set(gca,'Position',[pos(1),pos(2),pos1(3),pos(4)]);
line([0,0],get(gca,'YLim'),'Color',[0,0,0],'LineStyle','--')

%leaving
ax_leave=subplot(3,2,[2,4]);

t1 = min(time_leave);
t0 = 0;
t2 = max(time_leave);
i_t1 = 1;
i_t0 = find(time_leave>0,1,'first');
i_t2=length(time_leave);

imagesc(PlotDataLeave)
set(gca,'XTick',[i_t1,i_t0,i_t2],'XTickLabel',[t1,t0,t2])
cc=colorbar();
ylabel(cc,'dF/F (%)')
xlabel('Time from leaving (s)')
% ylabel('dF/F (%)')

%ticks in leave plot
axes(ax_leave)
PlotDataLeaveCorr=PlotDataLeave;
responsealigned = -wt(LeaveTrialIndex) ;
responsealigned= responsealigned(sortileave);
responsealigned_xidx = round((responsealigned+6) * 10)+1;
responsealigned_xidx(responsealigned_xidx<=0)=NaN;
for n =1:length(responsealigned_xidx)
    line([responsealigned_xidx(n),responsealigned_xidx(n)],[n-1,n],'Color','k','LineWidth',0.25)
    PlotDataLeaveCorr(n,1:min([responsealigned_xidx(n),size(PlotDataLeaveCorr,2)]))=NaN;
end

subplot(3,2,6)
y=nanmean(PlotDataLeaveCorr);
plot(time_leave,y,'-r')
set(gca,'XTick',[t1,t0,t2]);
set(gca,'YLim',[min(y)-0.01,max(y)+0.01])
ylabel('dF/F (%)')
xlabel('Time from leaving (s)')
pos1 = get(ax1,'Position');
pos=get(gca,'Position');
set(gca,'Position',[pos(1),pos(2),pos1(3),pos(4)]);
line([0,0],get(gca,'YLim'),'Color',[0,0,0],'LineStyle','--')

RedoTicks(gcf)
setfontline(8,2,'Arial')




if ~isempty(OUTNAME)
    writefigs(gcf,fullfile(basepath,animal,OUTNAME))
end

%% Plot waiting time distribution
if channel ==1
figure('Color',[1,1,1])
histogram(wt,50,'FaceColor',[0,0,0])

ylabel('n');xlabel('Waiting time (s)')

RedoTicks(gcf)
if ~isempty(OUTNAME)
    writefigs(gcf,fullfile(basepath,animal,OUTNAME))
end
end

