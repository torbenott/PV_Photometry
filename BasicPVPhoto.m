%basic photometry session analysis

%load a bpod session, ie struct SessionData
basepath = 'C:\Data\DataPostdoc\PV-Photometry';
animal = 'T05';
session = 'T05_NosePoke_Feb07_2019_Session1.mat';
OUTNAME = 'Photometry.pdf';

%params
channel = 2; %1=green, 2=red

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
DemodPhotoData=cell(1,nTrials);
DemodPhotoData2=cell(1,nTrials);
for n = 1:nTrials
    
    data = PhotoData{n};
    stateteimes = SessionData.RawEvents.Trial{n}.States;
    % rewarded
    if SessionData.Custom.Rewarded(n)
        if   SessionData.Custom.ChoiceLeft(n)==1
            name = 'water_L';
        elseif SessionData.Custom.ChoiceLeft(n) == 0
            name = 'water_R';
        else
            DemodPhotoData{n}=nan(150,3);
        end
        
        timetozero = stateteimes.(name)(1);
        Mod = Nidaq_modulation(modAmp,modFreq,TaskParameters);
        
        [currentNidaq1, rawNidaq1]=NidaqDemod(data(:,1),Mod,modFreq,modAmp,timetozero,TaskParameters);
        
        DemodPhotoData{n} = currentNidaq1;
    end
    
    %leaving
    if ~SessionData.Custom.Rewarded(n) && ~isnan(SessionData.Custom.ChoiceLeft(n))
        
        
        timetozero = stateteimes.ITI(1);
        Mod = Nidaq_modulation(modAmp,modFreq,TaskParameters);
        
        [currentNidaq1, rawNidaq1]=NidaqDemod(data(:,1),Mod,modFreq,modAmp,timetozero,TaskParameters);
        
        DemodPhotoData2{n} = currentNidaq1;
    end
    
end

%remove empty trials
DemodPhotoData(cellfun(@isempty,DemodPhotoData))=[];
DemodPhotoData2(cellfun(@isempty,DemodPhotoData2))=[];

%prepare for plotting
%reward
minT = -2;
maxT = 2;
tt = (maxT-minT)*10;
PlotData = nan(length(DemodPhotoData),tt);
for n = 1:length(DemodPhotoData)
    data = DemodPhotoData{n};
    ii = find(data(:,1)>0,1,'first');
    
    if ~isempty(ii) && ii>1
        
        mini = max([1,ii+minT*10]);
        maxi = min([size(data,1),ii+maxT*10]);
        
        lower = -minT*10 - (ii - mini)  + 1;
        upper = maxi-ii +   maxT*10 ;
        
        PlotData(n,lower:upper) = data(mini+1:maxi,3)';
    end
end
%leaving
minT = -2;
maxT = 2;
tt = (maxT-minT)*10;
PlotData2 = nan(length(DemodPhotoData2),tt);
for n = 1:length(DemodPhotoData2)
    data = DemodPhotoData2{n};
    ii = find(data(:,1)>0,1,'first');
    
    if ~isempty(ii) && ii>1
        
        mini = max([1,ii+minT*10]);
        maxi = min([size(data,1),ii+maxT*10]);
        
        lower = -minT*10 - (ii - mini)  + 1;
        upper = maxi-ii +   maxT*10 ;
        
        PlotData2(n,lower:upper) = data(mini+1:maxi,3)';
    end
end

%delete trials with many NaNs
PlotData(sum(isnan(PlotData),2)>30,:)=[];
PlotData2(sum(isnan(PlotData2),2)>30,:)=[];

%normalize per trial
PlotData = PlotData - repmat( nanmean(PlotData,2),1,size(PlotData,2));
PlotData2 = PlotData2 - repmat( nanmean(PlotData2,2),1,size(PlotData2,2));

time = linspace(minT,maxT,tt);

%% plot

figure('Color',[1,1,1],'Position',[     680   381   660   597])
%reward
subplot(3,2,[1,3])
imagesc(PlotData)
set(gca,'XTick',[1,21,40],'XTickLabel',[-2,0,2])
cc = colorbar();
ylabel(cc,'dF/F')
xlabel('Time from reward (s)')
ylabel('Trials')
uicontrol('Style','text','String',strrep(session(1:end-4),'_','-'),'FontName','Arial','Position',[10,580,200,10],'BackgroundColor',[1,1,1])
title(chname)

subplot(3,2,5)

plot(time,nanmean(PlotData),'-r')
set(gca,'XTick',[-2,0,2]);
ylabel('dF/F')
xlabel('Time from reward (s)')

%leaving
subplot(3,2,[2,4])
imagesc(PlotData2)
set(gca,'XTick',[1,21,40],'XTickLabel',[-2,0,2])
cc=colorbar();
ylabel(cc,'dF/F')
xlabel('Time from leaving (s)')
% ylabel('dF/F')


subplot(3,2,6)

plot(time,nanmean(PlotData2),'-r')
set(gca,'XTick',[-2,0,2]);
ylabel('dF/F')
xlabel('Time from leaving (s)')

RedoTicks(gcf)
setfontline(8,2,'Arial')

if ~isempty(OUTNAME)
    writefigs(gcf,fullfile(basepath,animal,OUTNAME))
end