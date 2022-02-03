%basic photometry session analysis

%load a bpod session, ie struct SessionData
basepath = 'C:\Data\DataPostdoc\PV-Photometry';
animal = 'JP01';
session = 'JP01_Dual2AFC_Apr15_2019_Session1.mat';
OUTNAME = '_1.pdf';

%params
channel = 1; %1=green, 2=red
% RealignLeaving=false;
SortByWT = true;%for leaaving plot
SortByReward=true;%for reward plot
ChoiceInvestmentOnly = true; %for choice-aligned plot, consider only investment (no-reward) trials
minWT = 2; %look only at waiting time higher than that (for all plots)
minRewardDelay = 1; %for average plot, consider only trials with minimum reward delay

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
file = fullfile(basepath,animal,'Dual2AFC','Session Data',session);
load(file)

nTrials = SessionData.nTrials-1;

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
DemodPhotoDataChoice=cell(1,nTrials);
DemodPhotoDataReward=cell(1,nTrials);
DemodPhotoDataLeave=cell(1,nTrials);
RewardStart = nan(1,nTrials);
ResponseStart = nan(1,nTrials);
BaselineChoice= nan(1,nTrials);
for n = 1:nTrials
    
    data = PhotoData{n};
    statetimes = SessionData.RawEvents.Trial{n}.States;
    
    %choice
    if ~isnan(SessionData.Custom.ChoiceLeft(n))
        timetozero = statetimes.wait_Sin(2);
        ResponseStart(n) = statetimes.wait_Sin(2);
        Mod = Nidaq_modulation(modAmp,modFreq,TaskParameters);
        
        [currentNidaq1, rawNidaq1]=NidaqDemod(data(:,1),Mod,modFreq,modAmp,timetozero,TaskParameters);
        
        DemodPhotoDataChoice{n} = currentNidaq1;  
        
        BaselineChoice(n)= quantile(DemodPhotoDataChoice{n}(1:20,2),0.25);
        
    end
    
    % reward
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
        
        Mod = Nidaq_modulation(modAmp,modFreq,TaskParameters);
        
        [currentNidaq1, rawNidaq1]=NidaqDemod(data(:,1),Mod,modFreq,modAmp,timetozero,TaskParameters);
        
        DemodPhotoDataReward{n} = currentNidaq1;
        
    end
    
    %leaving
    if SessionData.Custom.CatchTrial(n)==1 || SessionData.Custom.ChoiceCorrect(n) == 0
        
        timetozero=SessionData.Custom.ResolutionTime(n);
      
        Mod = Nidaq_modulation(modAmp,modFreq,TaskParameters);
        
        [currentNidaq1, rawNidaq1,base]=NidaqDemod(data(:,1),Mod,modFreq,modAmp,timetozero,TaskParameters);
        
        DemodPhotoDataLeave{n} = currentNidaq1;
        Baseline(n) = base;
    end
    
end

%renormalize?
for n =1:nTrials
    if ~isempty(DemodPhotoDataChoice{n})
        DemodPhotoDataChoice{n}(:,3)= (DemodPhotoDataChoice{n}(:,2) - BaselineChoice(n)) / BaselineChoice(n) * 100;
    end
    if ~isempty(DemodPhotoDataReward{n})
     DemodPhotoDataReward{n}(:,3)= (DemodPhotoDataReward{n}(:,2) - BaselineChoice(n)) / BaselineChoice(n) * 100;
    end
    if ~isempty(DemodPhotoDataLeave{n})
     DemodPhotoDataLeave{n}(:,3)= (DemodPhotoDataLeave{n}(:,2) - BaselineChoice(n)) / BaselineChoice(n) * 100;
    end
end

ChoiceTrialIndex = 1:nTrials;
LeaveTrialIndex = 1:nTrials;
RewardTrialIndex=1:nTrials;

RewardDelay = SessionData.Custom.FeedbackTime(1:nTrials);


%% remove empty trials
delI = cellfun(@isempty,DemodPhotoDataReward);
DemodPhotoDataReward(delI)=[];
RewardTrialIndex(delI)=[];

delI = cellfun(@isempty,DemodPhotoDataLeave);
DemodPhotoDataLeave(delI)=[];
LeaveTrialIndex(delI)=[];

delI = find(cellfun(@isempty,DemodPhotoDataChoice));
if ChoiceInvestmentOnly
    delI2 = find(~(SessionData.Custom.CatchTrial(1:nTrials) == 1 | SessionData.Custom.ChoiceCorrect(1:nTrials)==0));
    delI=union(delI,delI2);
end
DemodPhotoDataChoice(delI)=[];
ChoiceTrialIndex(delI)=[];


%% prepare for plotting
%choice
minT=-2;
maxT=6;
tt = (maxT-minT)*10;
time_choice = linspace(minT,maxT,tt);
PlotDataChoice = nan(length(DemodPhotoDataChoice),tt);
for n = 1:length(DemodPhotoDataChoice)
    data = DemodPhotoDataChoice{n};
    ii = find(data(:,1)>0,1,'first')-1;
    
    if ~isempty(ii) && ii>1
        
        mini = max([1,ii+minT*10]);
        maxi = min([size(data,1),ii+maxT*10]);
        
        lower = -minT*10 - (ii - mini)  + 1;
        upper =  (maxi-ii) - minT*10  ;
        
        PlotDataChoice(n,lower:upper) = data(mini+1:maxi,3)';
    end
end

%reward
minT = -2;
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
PlotDataLeave = nan(length(DemodPhotoDataLeave),tt);
for n = 1:length(DemodPhotoDataLeave)
    data = DemodPhotoDataLeave{n};
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
delI = sum(isnan(PlotDataChoice),2)>30;
PlotDataChoice(delI,:)=[];
ChoiceTrialIndex(delI)=[];

delI = sum(isnan(PlotDataReward),2)>30;
PlotDataReward(delI,:)=[];
RewardTrialIndex(delI)=[];

delI = sum(isnan(PlotDataLeave),2)>30 | RewardDelay(LeaveTrialIndex)'<minWT;
PlotDataLeave(delI,:)=[];
LeaveTrialIndex(delI)=[];

%normalize per trial 

% baselineidx = 1:20;
% PlotDataChoice = PlotDataChoice - repmat( nanmean(PlotDataChoice(:,baselineidx),2),1,size(PlotDataChoice,2));
% PlotDataReward = PlotDataReward - repmat( nanmean(PlotDataReward(:,baselineidx),2),1,size(PlotDataReward,2));
% PlotDataLeave = PlotDataLeave - repmat( nanmean(PlotDataLeave(:,baselineidx),2),1,size(PlotDataLeave,2));




%resort data?
if SortByReward
    [rewardsort,sortireward]=sort(RewardDelay(RewardTrialIndex));
    PlotDataReward=PlotDataReward(sortireward,:);
end
if SortByWT
    [wtsort,sortileave]=sort(RewardDelay(LeaveTrialIndex));
    PlotDataLeave=PlotDataLeave(sortileave,:);
    
    if ChoiceInvestmentOnly
    [wtsort,sortichoice]=sort(RewardDelay(ChoiceTrialIndex));
    PlotDataChoice=PlotDataChoice(sortichoice,:);    
    end
end

%% plot

figure('Color',[1,1,1],'Position',[      416         432        1022         485])

%choice
ax_choice=subplot(3,3,[1,4]); hold on

t1 = min(time_choice);
t0 = 0;
t2 = max(time_choice);
i_t1 = 1;
i_t0 = find(time_choice>0,1,'first');
i_t2=length(time_choice);

imagesc(PlotDataChoice)
set(gca,'XTick',[i_t1,i_t0,i_t2],'XTickLabel',[t1,t0,t2],'YLim',[0,size(PlotDataChoice,1)],'XLim',[1,size(PlotDataChoice,2)])
ax_choice.YAxis.Direction='reverse';
cc = colorbar();
ylabel(cc,'dF/F (%)')
xlabel('Time from choice (s)')
ylabel('Trials')
uicontrol('Style','text','String',strrep(session(1:end-4),'_','-'),'FontName','Arial','Position',[10,580,200,10],'BackgroundColor',[1,1,1])
% title(chname)

%ticks in choice plot
axes(ax_choice)
PlotDataChoiceCorr=PlotDataChoice;
responsealigned =  RewardDelay(ChoiceTrialIndex);
responsealigned= responsealigned(sortichoice);
responsealigned_xidx = round((responsealigned+2) * 10)+1;
responsealigned_xidx(responsealigned_xidx<=0)=NaN;
for n =1:length(responsealigned_xidx)
    line([responsealigned_xidx(n),responsealigned_xidx(n)],[n-1,n],'Color','k','LineWidth',0.25)
    if ~isnan(responsealigned_xidx(n))
    PlotDataChoiceCorr(n,min([responsealigned_xidx(n),size(PlotDataChoiceCorr,2)]):end)=NaN;
    end
end

%mean plot
ax2=subplot(3,3,7);
y=nanmean(PlotDataChoiceCorr);
plot(time_choice,y,'-r')
set(gca,'XTick',[t1,t0,t2]);
set(gca,'YLim',[min(y)-0.01,max(y)+0.01])
ylabel('dF/F (%)')
xlabel('Time from choice (s)')
pos1 = get(ax_choice,'Position');
pos=get(gca,'Position');
set(gca,'Position',[pos(1),pos(2),pos1(3),pos(4)]);
line([0,0],get(gca,'YLim'),'Color',[0,0,0],'LineStyle','--')

%reward
ax_reward=subplot(3,3,[2,5]); hold on

t1 = min(time_reward);
t0 = 0;
t2 = max(time_reward);
i_t1 = 1;
i_t0 = find(time_reward>0,1,'first');
i_t2=length(time_reward);

imagesc(PlotDataReward)
set(gca,'XTick',[i_t1,i_t0,i_t2],'XTickLabel',[t1,t0,t2],'YLim',[0,size(PlotDataReward,1)],'XLim',[1,size(PlotDataReward,2)])
ax_reward.YAxis.Direction='reverse';
cc = colorbar();
ylabel(cc,'dF/F (%)')
xlabel('Time from reward (s)')
ylabel('Trials')
uicontrol('Style','text','String',strrep(session(1:end-4),'_','-'),'FontName','Arial','Position',[10,580,200,10],'BackgroundColor',[1,1,1])
% title(chname)

%ticks in reward plot
axes(ax_reward)
PlotDataRewardCorr=PlotDataReward;
responsealigned = ResponseStart(RewardTrialIndex) - RewardStart(RewardTrialIndex);
responsealigned= responsealigned(sortireward);
responsealigned_xidx = round((responsealigned+2) * 10)+1;
responsealigned_xidx(responsealigned_xidx<=0)=NaN;
for n =1:length(responsealigned_xidx)
    line([responsealigned_xidx(n),responsealigned_xidx(n)],[n-1,n],'Color','k','LineWidth',0.25)
    if ~isnan(responsealigned_xidx(n))
    PlotDataRewardCorr(n,1:min([responsealigned_xidx(n),size(PlotDataRewardCorr,2)]))=NaN;
    end
end


ax2=subplot(3,3,8);
y=nanmean(PlotDataRewardCorr);
plot(time_reward,y,'-r')
set(gca,'XTick',[t1,t0,t2]);
set(gca,'YLim',[min(y)-0.01,max(y)+0.01])
ylabel('dF/F (%)')
xlabel('Time from reward (s)')
pos1 = get(ax_reward,'Position');
pos=get(gca,'Position');
set(gca,'Position',[pos(1),pos(2),pos1(3),pos(4)]);
line([0,0],get(gca,'YLim'),'Color',[0,0,0],'LineStyle','--')

%leaving
ax_leave=subplot(3,3,[3,6]);

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
responsealigned = ResponseStart(LeaveTrialIndex) - SessionData.Custom.ResolutionTime(LeaveTrialIndex) ;
responsealigned= responsealigned(sortileave);
responsealigned_xidx = round((responsealigned+6) * 10)+1;
responsealigned_xidx(responsealigned_xidx<=0)=NaN;
for n =1:length(responsealigned_xidx)
    line([responsealigned_xidx(n),responsealigned_xidx(n)],[n-1,n],'Color','k','LineWidth',0.25)
    if ~isnan(responsealigned_xidx(n))
    PlotDataLeaveCorr(n,1:min([responsealigned_xidx(n),size(PlotDataLeaveCorr,2)]))=NaN;
    end
end

subplot(3,3,9)
y=nanmean(PlotDataLeave);
plot(time_leave,y,'-r')
set(gca,'XTick',[t1,t0,t2]);
set(gca,'YLim',[min(y)-0.01,max(y)+0.01])
ylabel('dF/F (%)')
xlabel('Time from leaving (s)')
pos1 = get(ax_leave,'Position');
pos=get(gca,'Position');
set(gca,'Position',[pos(1),pos(2),pos1(3),pos(4)]);
line([0,0],get(gca,'YLim'),'Color',[0,0,0],'LineStyle','--')

RedoTicks(gcf)
setfontline(8,2,'Arial')



if ~isempty(OUTNAME)
    writefigs(gcf,fullfile(basepath,animal,[session(1:end-4),OUTNAME]))
end

figure
hold on
trials1=1:40;
trials2=41:70;
trials3=71:size(PlotDataChoiceCorr,1);
plot(nanmean(PlotDataChoiceCorr(trials1,:)));
plot(nanmean(PlotDataChoiceCorr(trials2,:)));
plot(nanmean(PlotDataChoiceCorr(trials3,:)));

%check reward response
RewardResponse = mean(PlotDataReward(:,21:31),2);
TimeInTrial = RewardStart(RewardTrialIndex)-ResponseStart(RewardTrialIndex) ;
figure
scatter(TimeInTrial,RewardResponse)

%check pre-leave peak
PlotDataLeaveCorr_smooth=[];
for n =1:size(PlotDataLeaveCorr,1)
PlotDataLeaveCorr_smooth(n,:) = smoothdata(PlotDataLeaveCorr(n,:),'gaussian',5);
end

PlotDataLeaveCorr_smooth(PlotDataLeaveCorr_smooth==0)=NaN;

d_Leave = diff(PlotDataLeaveCorr_smooth,1,2);
d_neg = d_Leave(:,1:end-2) < 0 & d_Leave(:,2:end-1)<0 & d_Leave(:,3:end)<0;
d_neg = [d_neg,false(size(d_Leave,1),2)];
d_Leave(~d_neg)=NaN;

[~,peak_idx]=min(d_Leave(:,1:60),[],2);

% [~,peak_idx]=max(PlotDataLeaveCorr_smooth(:,1:60),[],2);
WT = RewardDelay(LeaveTrialIndex);
WT=WT(sortileave);
figure
scatter(WT,-time_leave(peak_idx))
ylabel('max response')
xlabel('WT')
pp=polyfit(WT,-time_leave(peak_idx),1);
hold on
xx=linspace(min(WT),max(WT),100);
plot(xx,pp(1).*xx+pp(2))

figure('Color',[1,1,1])
BaselineChoice(1)=[];
iii = find(~isnan(BaselineChoice));
plot(iii,BaselineChoice(~isnan(BaselineChoice)),'-k');
xlabel('Trials')
ylabel('Baseline F')

% %% Plot waiting time distribution
% if channel ==1
% figure('Color',[1,1,1])
% histogram(wt,50,'FaceColor',[0,0,0])
% 
% ylabel('n');xlabel('Waiting time (s)')
% 
% RedoTicks(gcf)
% if ~isempty(OUTNAME)
%     writefigs(gcf,fullfile(basepath,animal,OUTNAME))
% end
% end

