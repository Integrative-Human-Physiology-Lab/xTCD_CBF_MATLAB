function align_labchart_velocity (velocity, diameter, labchart, destination_folder_name)


%clear;

%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%destination_folder_name = 'E:\Min and Steph\subjects_data\final_data_2018';
% load velocity
velocity_data = load(velocity);
timeDVD = velocity_data.time_array ;
timeDVD = timeDVD';
timeDVD_original = timeDVD;
timeDVD = timeDVD(:,1);
datDvd = velocity_data.Velocity_forRRI;
sd=size(datDvd);
%remove NaN values before interpolation;
bad = find(isnan(timeDVD));
if isempty(bad) == 0
    timeDVD(bad)=[];
    datDvd(bad)=[];
end
bad = find(isnan(datDvd));
if isempty(bad) == 0
    timeDVD(bad)=[];
    datDvd(bad)=[];
end
% interpolate velocity data to median Sample rate of current signal to fill
% in any missing values
timeMatraw = timeDVD; %time signal produced from matlab velocity analysis
VelMatraw = datDvd; % Velocity derived from the matlab velocity analysis
timeMatdiff = median(diff(timeMatraw));
timeMat = min(timeMatraw):timeMatdiff:max(timeMatraw);
VelMat = interp1(timeMatraw,VelMatraw,timeMat,'linear');
Vel_SR = 1 / diff(timeDVD(1:2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load diameter
diameter_data = xlsread(diameter);
diam = diameter_data(:,6)'; %mm
time_diam = diameter_data(:,5)'; %sec
% assign variables that have been produced by brachial analyzer software
timeBA=time_diam; %time based on the brachial analzyer
%temporary correction since there seems to be a 150 ms delay for maui
%analysis to occur compared to matlab velcoity
timeBA = timeBA -0.15;
diamBA=diam/10; %diameter (in cm)  derived from brachial analyzer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load labchart
load(labchart);
SR = round(1./mean(diff(t)));
if SR == 100
    BVcm_int = BVcm(ikp,:); %added by BS, never indicated BVcm of interest only BV
elseif SR ==1000
end

t_test=t';
if length(ECG) ~= length(t_test)
    ECG = ECG(ikp);
end

ECG_new=ECG;
ECG_new=ECG_new./max(ECG_new);
ECG_new=((max(BP)-min(BP)).*ECG_new)+min(BP);

dumy_BP=BP';
ECG_new=ECG_new';
BP_delay_shift=0;
BP_shifter=0;
ECG_orig1 = ECG_new;

%% Plot all signals together
txtfontsize= 11;
t_name1 = ' ' ;
t_name2 = ' ' ;
figure(1)
clf
% Set default figure size to near maximum of screen size
scrsz = get(0,'ScreenSize');
position_default = [0.01*scrsz(3) 0.07*scrsz(4) 0.98*scrsz(3) 0.85*scrsz(4)];
set(gcf,'position',position_default)

% Plot parameters
fs = 'fontsize';            % save fontsize field property
fs_title = 14;              % fontsize for titles
fs_label = 12;              % fontsize for axis labels
lw = 'LineWidth';           % save linewidth field property
lw_thin = 0.3;              % linewidth for thin lines
lw_thick = 2;               % linewidth for thicker lines
ms = 'markersize';          % save markersize field property
ms_large = 10;              % markersize for large points
k =num_artery+3;

subplot(k,1,1)
hold on;
plot(t,bp_f, 'b', lw, lw_thin);
plot(tRRi,MBP,'k', lw, lw_thick);
plot(tRRi,SBP,'r', lw, lw_thick);
plot(tRRi,DBP,'r', lw, lw_thick);
hold off;
axis tight
yl = ylim;
yl(1) = yl(1)-yl(2)*.05;
yl(2) = yl(2)+yl(2)*.05;
ylim(yl);
ylabel({'Arterial';'Blood Pressure';'(mmHg)'}, fs, fs_label)
title(regexprep(plot_title, '_', ' '), fs, fs_title)

for figcont = 1:num_artery
    subplot(k,1,figcont+1)
    hold on;
    plot(t,bv_f(:,figcont), 'g', lw, lw_thin);
    plot(tRRi,MBV(:,figcont),'k', lw, lw_thick);
    plot(tRRi,SBV(:,figcont),'r', lw, lw_thick);
    plot(tRRi,DBV(:,figcont),'r', lw, lw_thick);
    line([0 t(end)], [100 100], 'Color','m');
    hold off;
    axis tight;
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    xlabel('Time (s)', fs, fs_label);
    ylabel({cell2mat(labels (figcont));'(%)'}, fs, fs_label);
end

subplot(k,1,num_artery+2)
plot(tRRi,1000./RRi*60,'r', lw, lw_thick);
axis tight
yl = ylim;
yl(1) = yl(1)-yl(2)*.05;
yl(2) = yl(2)+yl(2)*.05;
ylim(yl);
ylabel('HR (bpm)', fs, fs_label);

subplot(k,1,num_artery+3)
plot(tco2,CO2, 'g', lw, lw_thin), hold on;
plot(tetCO2,etCO2,'r', lw, lw_thick), hold off;
axis tight
yl = ylim;
yl(1) = yl(1)-yl(2)*.05;
yl(2) = yl(2)+yl(2)*.05;
ylim(yl);
ylabel('etCO2 (mmHg)', fs, fs_label);
if ~isempty(Notes_list)
    for no =1:length(Notes_list )
        Notes(no) =  Notes_list{no};
    end
else
    Notes = 'No notes';
end
Notes_cont = uicontrol(gcf,'Style','text','String',Notes,'units','normalized','Position',[.915 .11 .08 .85]);
% Resize figure for full page printing
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.25 0.25 8 10.5]);

%% Plot blood pressure and ECG data
% ECG should lead Blood Pressure (BP)
button_2=0;
dx=3;
ha = [0 t(end)];
basefilter1 = 0;

TIMEP_length = length(timeP);
TIMETEST_length = length(t_test);

while button_2 == 0
    handles.H = figure(2);
    clf
    a=gca;
    ha_c=0;
    dx_width=0;
    hplot=plot(t_test,dumy_BP,'b'); hold on;
    if length(ECG_new) == TIMEP_length
        plot(timeP,ECG_new,'g');
        t_test = timeP;
    elseif length(ECG_new) == TIMETEST_length
        plot(t_test,ECG_new,'g');
    end
    legend('BP','ECG');
    xlabel('Time');
    title('Check for ECG to BP Alignment'); %scrollplot
    
    set(gcf,'doublebuffer','on');
    % This avoids flickering when updating the axis
    xlim([ha]);
    % Generate constants for use in uicontrol initialization
    pos=get(a,'position');
    Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
    % This will create a slider which is just underneath the axis
    % but still leaves room for the axis labels above the slider
    xmax=max(t_test);
    S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
    pb1au = uicontrol(gcf,'Style','pushbutton','String','Shift+','units','normalized','Position',[.925 .90 .06 .03],'Callback','dumy_BP=circshift(dumy_BP,round(SR*.01));BP_delay_shift=BP_delay_shift+round(SR*.01);BP_shifter = BP_shifter+1; uiresume;');
    pb2au = uicontrol(gcf,'Style','pushbutton','String','Shift-','units','normalized','Position',[.925 .86 .06 .03],'Callback', 'dumy_BP=circshift(dumy_BP,-round(SR*.01));BP_delay_shift =BP_delay_shift-round(SR*.01);BP_shifter = BP_shifter-1;  uiresume;');
    pb3au = uicontrol(gcf,'Style','pushbutton','String','Shift++','units','normalized','Position',[.925 .82 .06 .03],'Callback','dumy_BP=circshift(dumy_BP,round(SR*.1));BP_delay_shift=BP_delay_shift+round(SR*.1);BP_shifter = BP_shifter+10; uiresume;');
    pb4au = uicontrol(gcf,'Style','pushbutton','String','Shift--','units','normalized','Position',[.925 .78 .06 .03],'Callback','dumy_BP=circshift(dumy_BP,-round(SR*.1));BP_delay_shift=BP_delay_shift-(round(SR*.1)); BP_shifter = BP_shifter-10; uiresume;');
    pb5au = uicontrol(gcf,'Style','pushbutton','String','ECG Invert','units','normalized','Position',[.925 .71 .06 .03],'Callback','ECG_new = ECG_new*-1; uiresume;');
    pb6au = uicontrol(gcf,'Style','pushbutton','String','ECG offset+','units','normalized','Position',[.925 .68 .06 .03],'Callback','ECG_new = ECG_new+1; uiresume;');
    pb7au = uicontrol(gcf,'Style','pushbutton','String','ECG offset-','units','normalized','Position',[.925 .65 .06 .03],'Callback','ECG_new = ECG_new-1; uiresume;');
    pb8au = uicontrol(gcf,'Style','pushbutton','String','Base Filter','units','normalized','Position',[.925 .60 .06 .03],'Callback','basefilter1=1; uiresume;');
    pb9au = uicontrol(gcf,'Style','pushbutton','String','Undo Filter','units','normalized','Position',[.925 .57 .06 .03],'Callback','ECG_new=ECG_orig1; uiresume;');
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Done','units','normalized','Position',[.925 .11 .06 .03],'Callback','button_2 = 1; uiresume;');
    pbreset = uicontrol(gcf,'Style','pushbutton','String','Reset','units','normalized','Position',[.925 .2 .06 .03],'Callback','ha_c=1; uiresume;');
    pbwin = uicontrol(gcf,'Style','pushbutton','String','Window +','units','normalized','Position',[.925 .3 .06 .03],'Callback','dx=dx+2; dx_width=2; uiresume;');
    pbwin = uicontrol(gcf,'Style','pushbutton','String','Window -','units','normalized','Position',[.925 .34 .06 .03],'Callback','dx=dx-2; dx_width=-2; uiresume;');
    h=uicontrol('style','slider','units','normalized','position',Newpos,'callback',S,'min',0,'max',xmax-dx);
    uiwait(gcf)
    if  ha_c==0
        ha=get(gca,'XLim');
        ha(2)=ha(2)+dx_width;
        xlim([ha]);
    elseif  ha_c==1
        ha=[0 t(end)];
    end
    if basefilter1 == 1
        basefilter1 = 0;
        [bb,aa] = butter(8,5/(SR/2),'high');
        ECG_new = filtfilt(bb,aa,ECG_new);
    end
    set(hplot,'Visible','off');
end
close(handles.H)
BP=dumy_BP';

%% align data based on marker if possible
disp('Proceeding to Peak Detection')
timemarker = timeMat; % Matlab Time marker
SignalMarker = VelMat; % Matlab Velocity marker Signal

%Align diameter with Velocity signal if using velocity signal from matlab
Vel = VelMat;
diam = diamBA;

n = size (diam);
marker = NaN(n);
markerChart = NaN(n);

% Align diam with Vel signal
tdiam = timeBA;

%Trim Data to just include aligned diam and Vel data
st1 = max([t_1000HZ(1) timemarker(1)]);
st2 = min([t_1000HZ(end) timemarker(end)]);
good = find(tdiam>=st1 & tdiam<=st2);
diam = diamBA;
diam = diam(good);
tdiam = tdiam(good);
clear good
good = find(timemarker>=st1 & timemarker<=st2);
Velocity_PeakDetect = SignalMarker(good);
timeMarker = timemarker(good);
clear good

m = 1; %value preset for use in basefilter (Highpass) function.
n = 1; %value preset for use in LPF function.
timeMat = timeMarker;
timeMatSound = timeMat; % Created new variable for the purpose of being able to select Velocity Sound Signal during loop.
tRRi_corr = tRRi;
SR = round(1./mean(diff(timeMat)));

points = ceil(0.2*SR);
thresh = nanstd(Velocity_PeakDetect);
maxthresh = nanmax(Velocity_PeakDetect)*1.1;
minthresh = min(0,min(Velocity_PeakDetect));
Value = 'value';
Remove_Group = 0;
button = 0;
basefilter = 0;
LPfilter = 0;

iRRiV = Peaksmsa(Velocity_PeakDetect,points,thresh);
XL = [0 max(timeMat)];
XL2 = [0 max(timeMat)];
redetect = 0;
HRV = 1./diff(timeMat(iRRiV))*60;
Velocityorig = Velocity_PeakDetect;
k=50;
SR = round(1./mean(diff(timeMat)));
cutoff_freq = SR/4;
thresh = nanstd(Velocity_PeakDetect);
points = ceil(0.2*SR);

figure(3)
while button==0
    
    if redetect == 1
        iRRiV = Peaksmsa(Velocity_PeakDetect,points,thresh);
        bad = find(diff(timeMat(iRRiV))<0.15);
        iRRiV(bad+1) = [];
        HRV = 1./diff(timeMat(iRRiV))*60;
        XL = [0 max(timeMat)];
        XL2 = [0 max(timeMat)];
        redetect = 0;
    end
    
    clf
    subplot(3,1,1:2)
    plot(timeMat,Velocity_PeakDetect)
    ha = get(gcf,'CurrentAxes');
    hold on
    hline = refline(0,thresh);
    set(hline,'Color','g');
    ylabel('Velocity')
    plot(timeMat(iRRiV),Velocity_PeakDetect(iRRiV),'ro','linewidth',2 )
    hold off
    title ('Peak Detection of Velocity');
    
    xlim(XL)
    YL1 = ylim;%gets the YLim or YLimMode property of an axes
    
    subplot(3,1,3)
    plot(timeMat(iRRiV(2:end)),HRV)
    ylabel('HRV')
    xlim(XL2)
    ha2 = get(gcf,'CurrentAxes');
    ss = get(gcf,'Position');
    
    pb1ut1 = uicontrol(gcf,'style','text','string',ceil(YL1(2)),'units','normalized','HorizontalAlignment','left','Position',[.96 .86 .025 .025]);
    pb1utime = uicontrol(gcf,'style','text','string',floor(YL1(1)),'units','normalized','HorizontalAlignment','left','Position',[.96 .5 .025 .025]);
    pb1ut3 = uicontrol(gcf,'style','text','string',round(thresh),'units','normalized','HorizontalAlignment','left','Position',[.96 .68 .025 .025]);
    pb1u = uicontrol(gcf,'Style','slider','Max',ceil(max(Velocity_PeakDetect)),'Min',floor(min(Velocity_PeakDetect)-max(Velocity_PeakDetect)),'Value',round(thresh),'SliderStep',[0.001 0.1],'units','normalized','Position',[.925 .5 .025 .4],'Callback','thresh = get(pb1u,Value); uiresume;');
    
    pb2u = uicontrol(gcf,'Style','pushbutton','String','Delay+','units','normalized','Position',[.925 .45 .05 .03],'Callback','points = points + round(points*.1); redetect = 1; uiresume;');
    pb2l = uicontrol(gcf,'Style','pushbutton','String','Delay-','units','normalized','Position',[.925 .41 .05 .03],'Callback','points = points - round(points*.1); redetect = 1; uiresume;');
    pb3u = uicontrol(gcf,'Style','pushbutton','String','Delay++','units','normalized','Position',[.925 .37 .05 .03],'Callback','points = points + round(points*.5); redetect = 1; uiresume;');
    pb3l = uicontrol(gcf,'Style','pushbutton','String','Delay--','units','normalized','Position',[.925 .33 .05 .03],'Callback','points = points - round(points*.5); redetect = 1; uiresume;');
    
    pb5u = uicontrol(gcf,'Style','pushbutton','String','Thresh+','units','normalized','Position',[.925 .29 .05 .03],'Callback','thresh = thresh + thresh*.05; redetect = 1; uiresume;');
    pb5l = uicontrol(gcf,'Style','pushbutton','String','Thresh-','units','normalized','Position',[.925 .25 .05 .03],'Callback','thresh = thresh - thresh*.05; redetect = 1; uiresume;');
    
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF','units','normalized','Position',[.925 .2 .05 .03],'Callback','n = n + 1;LPfilter = 1; redetect = 1; uiresume;');
    pb4l = uicontrol(gcf,'Style','pushbutton','String','BaseF','units','normalized','Position',[.925 .16 .05 .03],'Callback','m = m + 1; basefilter = 1; redetect = 1; uiresume;');
    pb6l = uicontrol(gcf,'Style','pushbutton','String','ReDetect','units','normalized','Position',[.925 .12 .05 .03],'Callback','redetect = 1; uiresume;');
    
    pb9 = uicontrol(gcf,'Style','pushbutton','String','Diameter','units','normalized','Position',[.8 .025 .05 .03],'Callback','Velocity_PeakDetect = diamBA; timeMat = timeBAvel; thresh = max(diamBA)*0.35; redetect = 1; uiresume;');
    pb7 = uicontrol(gcf,'Style','pushbutton','String','Velocity_Sound','units','normalized','Position',[.4 .025 .05 .03],'Callback','Velocity_PeakDetect = VelMat; timeMat = timeMatSound; thresh = max(VelMat)*0.35;  redetect = 1;uiresume;');
    pb8 = uicontrol(gcf,'Style','pushbutton','String','Velocity_Brachial Analyzer','units','normalized','Position',[.6 .025 .05 .03],'Callback','Velocity_PeakDetect = Velocity_BA; timeMat = timeBAvel; thresh = max(Velocity_BA)*0.35; redetect = 1; uiresume;');
    
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    
    uiwait(gcf)
    
    Velocityf = zeros(length(Velocity_PeakDetect),1);
    
    if basefilter == 1
        basefilter = 0;
        if m > 5
            m = 5;
        end
        [b,a] = butter(m,5/SR,'high');
        Velocityf = filtfilt(b,a,Velocity_PeakDetect);
        Velocity_PeakDetect = Velocityf;
    end
    
    if LPfilter == 1
        set(gcf,'Pointer','watch')
        if n > 8
            n = 8;
        end
        [b,a] = butter(n,cutoff_freq/(SR/2),'low');
        Velocityf = filtfilt(b,a,Velocity_PeakDetect);
        Velocity_PeakDetect = Velocityf;
        LPfilter = 0;
        cutoff_freq=cutoff_freq*.9;
        set(gcf,'Pointer','arrow')
    end
    XL = xlim(ha);
    
end
close(gcf)

% Add good beats, remove bad beats
closefig = get(0,'CurrentFigure');
close(closefig);
bad = find(diff(timeMat(iRRiV))<0.25);
iRRiV(bad+1) = [];
button = 0;
bad = [];
iRRiVold = iRRiV;
tRRiV = timeMat(iRRiV);
RRiV = diff(tRRiV)*1000;
HRV = 60000./RRiV;
timePV = timeMat(iRRiV)';
XL = [0 max(timeMat)];
XL2 = XL;
removed_peaks = iRRiV;
iRRiV_original = iRRiV;
notes=0;
Remove_Group = 0;
button = 0;
SR = round(1./mean(diff(timeMat)));

figure(4)
while button == 0
    
    clf
    subplot(3,1,1:2)
    plot(timeMat,Velocity_PeakDetect)
    hold on
    plot(timeMat(iRRiV),Velocity_PeakDetect(iRRiV),'ro','linewidth',2 )
    ylabel('Velocity')
    hold off
    title ('Add or remove peaks');
    
    xlim(XL)
    ha = get(gcf,'CurrentAxes');
    subplot(3,1,3)
    plot(timeMat(iRRiV(2:end)),HRV)
    ylabel('DVD RRI')
    xlim(XL2)
    ha2 = get(gcf,'CurrentAxes');
    ss = get(gcf,'Position');
    
    
    pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove','Position',[ss(3)*.925 ss(4)*.7 ss(3)*.05 ss(4)*.05],'Callback', ...
        '[x, y] = ginput(1); x = x(end); [mpt bad] = min(abs(timeMat(iRRiV)-x)); iRRiVold = iRRiV; iRRiV(bad) = []; HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
    pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','Position',[ss(3)*.925 ss(4)*.5 ss(3)*.05 ss(4)*.05],'Callback', ...
        '[x, y] = ginput(1); x = round(x(end).*SR); rows2 = find(x > iRRiV_original); iRRiVold = iRRiV; iRRiV(end+1) = iRRiV_original(rows2(end)); iRRiV = sort(iRRiV); iRRiV = unique(iRRiV); HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
    pb3u = uicontrol(gcf,'Style','pushbutton','String','Undo','Position',[ss(3)*.925 ss(4)*.3 ss(3)*.05 ss(4)*.05],'Callback', ...
        'iRRiV = iRRiVold; HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','AddClick','Position',[ss(3)*.925 ss(4)*.9 ss(3)*.05 ss(4)*.05],'Callback', ...
        '[x, y] = ginput(1); rows2 = find(x < timeMat); iRRiVold = iRRiV; iRRiV(end+1) = rows2(1); iRRiV = sort(iRRiV); iRRiV = unique(iRRiV); HRV = 1./diff(timeMat(iRRiV)).*60;  uiresume;');
    pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove Group','Position',[ss(3)*.025 ss(4)*.7 ss(3)*.1 ss(4)*.05],'Callback','Remove_Group = 1; uiresume;');
    pb6u = uicontrol(gcf,'Style','pushbutton','String','Undo All','Position',[ss(3)*.025 ss(4)*.3 ss(3)*.1 ss(4)*.05],'Callback', ...
        'iRRiV = iRRiV_original; HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
    pb7u = uicontrol(gcf,'Style','pushbutton','String','Notes','FontSize',12,'Position',[ss(3)*.025 ss(4)*.55 ss(3)*.1 ss(4)*.05],'BackgroundColor','g','Callback', ...
        'notes = 1; uiresume;');
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    
    uiwait(gcf)
    
    if notes ==1
        answer_notes = inputdlg('please enter your notes');
        Notes_list{length(Notes_list)+1}= answer_notes;
        notes =0;
    end
    
    if Remove_Group > 0
        iRRiVold = iRRiV;
        removed_peaks =iRRiV;
        k = waitforbuttonpress;
        point1 = get(gca,'CurrentPoint');    % button down detected
        finalRect = rbbox;                   % return figure units
        pointime = get(gca,'CurrentPoint');    % button up detected
        point1 = point1(1,1:2);              % extract x and y
        pointime = pointime(1,1:2);
        p1 = min(point1,pointime);             % calculate locations
        offset = abs(point1-pointime);         % and dimensions
        x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)]; % is the time
        y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        hold on
        axis manual
        fig = plot(x,y,'r','linewidth',2) ;
        
        Index_Start_Array=find(timeMat(iRRiV)<= min(x));
        Index_end_Array  = find(timeMat(iRRiV)<= max(x));
        
        if  (isempty(Index_Start_Array)) ||  (isempty(Index_end_Array) )
            display('');
        else
            Index_Start = Index_Start_Array(end);
            Index_end = Index_end_Array(end);
            
            removed_peaks(Index_Start:Index_end)=[];
            iRRiV = removed_peaks;
            
            HRV = 1./diff(timeMat(iRRiV)).*60;
            Remove_Group = 0;
            Index_Start_Array=[];
            Index_end_Array=[];
            removed_peaks=[];
        end
    end
    
    XL = xlim(ha);
    XL2 = xlim(ha2);
    
end
close(gcf)
ikpv = find(iRRiV<length(timeMat));
tRRiV = timeMat(iRRiV);
RRiV = diff(tRRiV)*1000;

HRV = 1./diff(tRRiV).*60;
SR = round(1./mean(diff(t)));
XL = [0 max(t)];
tHR = tRRi_old;
tHRV = timeMarker;
HRVShift = HRV;

%% plot heart rates to align
figure(5)
button = 0;
YL = [min([HR HRVShift]) max([HR HRVShift])];

while button == 0
    clf
    plot(tHR,HR, 'bo-')
    hold on
    plot(tHRV(iRRiV(1:end-1)),HRVShift, 'ro-')
    hold off
    title ('Align Using Heart Rates');
    xlim(XL);
    ylim(YL);
    xlabel('Time(Seconds)')
    legend('RRI LabChart', 'RRI Velocity');
    
    ha = get(gcf,'CurrentAxes');
    YL1 = ylim;
    
    ss = get(gcf,'Position');
    pb1u = uicontrol(gcf,'Style','slider','Max',round(max([HR HRVShift])),'Min',round(min([HR HRVShift])),'Value',round(YL(2)),'SliderStep',[0.01 0.1],'units','normalized','Position',[.05 .6 .025 .3],'Callback','YL(2) = get(pb1u,Value); ylim(YL); uiresume;');
    pb1ut3 = uicontrol(gcf,'style','text','string',round(YL(2)),'units','normalized','HorizontalAlignment','left','Position',[.05 .68 .025 .025]);
    pb1l = uicontrol(gcf,'Style','slider','Max',round(max([HR HRVShift])),'Min',round(min([HR HRVShift])),'Value',round(YL(1)),'SliderStep',[0.01 0.1],'units','normalized','Position',[.05 .2 .025 .3],'Callback','YL(1) = get(pb1l,Value); ylim(YL); uiresume;');
    pb1lt3 = uicontrol(gcf,'style','text','string',round(YL(1)),'units','normalized','HorizontalAlignment','left','Position',[.05 .28 .025 .025]);
    pb1au = uicontrol(gcf,'Style','pushbutton','String','Shift+','units','normalized','Position',[.925 .9 .05 .03],'Callback','tHRV = tHRV+0.1; uiresume;');
    pb2au = uicontrol(gcf,'Style','pushbutton','String','Shift-','units','normalized','Position',[.925 .85 .05 .03],'Callback', 'tHRV = tHRV-0.1;  uiresume;');
    pb3au = uicontrol(gcf,'Style','pushbutton','String','Shift+','units','normalized','Position',[.925 .775 .05 .03],'Callback','tHRV = tHRV+1; uiresume;');
    pb4au = uicontrol(gcf,'Style','pushbutton','String','Shift-','units','normalized','Position',[.925 .725 .05 .03],'Callback', 'tHRV = tHRV-1;  uiresume;');
    pb5au = uicontrol(gcf,'Style','pushbutton','String','Shift++','units','normalized','Position',[.925 .65 .05 .03],'Callback','tHRV = tHRV+10; uiresume;');
    pb6au = uicontrol(gcf,'Style','pushbutton','String','Shift--','units','normalized','Position',[.925 .60 .05 .03],'Callback','tHRV = tHRV-10; uiresume;');
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    
    uiwait(gcf)
    XL = xlim(ha);
end
close(gcf)

%% Further Align Velocity(or Diameter) and Signal based on previous shift
% decided to comment out since don't think we need a second alignment

button = 0;
delay_shift=0;

BP_Shiftnorm = zscore(bp_f);
Velocity_Shiftnorm = zscore(Velocity_PeakDetect);
figure(6)

while button == 0
    clf
    plot(t, BP_Shiftnorm, 'b-')
    hold on
    plot(tHRV,Velocity_Shiftnorm, 'r-')
    hold off
    xlim(XL);
    title('Velocity and BP Plot')
    xlabel('Time(Seconds)')
    legend('BP', 'Velocity');
    
    ha = get(gcf,'CurrentAxes');
    ss = get(gcf,'Position');
    
    pb5au = uicontrol(gcf,'Style','pushbutton','String','Shift+','units','normalized','Position',[.925 .85 .05 .03],'Callback','tHRV = tHRV+0.01;uiresume;');
    pb6au = uicontrol(gcf,'Style','pushbutton','String','Shift-','units','normalized','Position',[.925 .80 .05 .03],'Callback', 'tHRV = tHRV-0.01;  uiresume;');
    pb7au = uicontrol(gcf,'Style','pushbutton','String','Shift++','units','normalized','Position',[.925 .60 .05 .03],'Callback','tHRV = tHRV+0.05; uiresume;');
    pb8au = uicontrol(gcf,'Style','pushbutton','String','Shift--','units','normalized','Position',[.925 .55 .05 .03],'Callback','tHRV = tHRV-0.05; uiresume;');
    pb9au = uicontrol(gcf,'Style','pushbutton','String','Shift+++','units','normalized','Position',[.925 .35 .05 .03],'Callback','tHRV = tHRV+1.0; uiresume;');
    pb10au = uicontrol(gcf,'Style','pushbutton','String','Shift---','units','normalized','Position',[.925 .30 .05 .03],'Callback','tHRV = tHRV-1.0; uiresume;');
    pb11 = uicontrol(gcf,'Style','pushbutton','String','ECG inverse','units','normalized','Position',[.925  .15 .05 .03],'Callback','Signal = -Signal; uiresume;');
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    uiwait(gcf)
    XL = xlim(ha);
end

close(gcf)
timeAligned = tHRV;

%check for debugging purposes
time_dur_vel=tHRV(end)-tHRV(1);
time_dur_diam=timeBA(end)-timeBA(1);
ShiftTime = nanmean(tHRV(1)-t(1));
timeBA_orig = timeBA;
timeBA = timeBA + ShiftTime;

button = 0;
Diam_Shiftnorm = zscore(diamBA);

figure(7)
% Added 9/5/2019
Diam_Shiftnorm = interp1(timeBA,Diam_Shiftnorm,tHRV);
while button == 0
    
    
    clf
    plot(t, BP_Shiftnorm, 'b-')
    hold on
    % Added 9/5/2019 //////////////////////////////////////////////////////////
    plot(tHRV, Diam_Shiftnorm,'r')
    % /////////////////////////////////////////////////////////////////////////
    hold off
    axis tight
    xlim(XL);
    title('Diameter and BP Plot')
    xlabel('Time(Seconds)')
    legend('BP', 'Diameter');
    
    ha = get(gcf,'CurrentAxes');
    ss = get(gcf,'Position');
    
    pb5au = uicontrol(gcf,'Style','pushbutton','String','Shift+','units','normalized','Position',[.925 .85 .05 .03],'Callback','timeBA = timeBA+0.01;uiresume;');
    pb6au = uicontrol(gcf,'Style','pushbutton','String','Shift-','units','normalized','Position',[.925 .80 .05 .03],'Callback', 'timeBA = timeBA-0.01;  uiresume;');
    pb7au = uicontrol(gcf,'Style','pushbutton','String','Shift++','units','normalized','Position',[.925 .60 .05 .03],'Callback','timeBA = timeBA+0.05; uiresume;');
    pb8au = uicontrol(gcf,'Style','pushbutton','String','Shift--','units','normalized','Position',[.925 .55 .05 .03],'Callback','timeBA = timeBA-0.05; uiresume;');
    pb9au = uicontrol(gcf,'Style','pushbutton','String','Shift+++','units','normalized','Position',[.925 .35 .05 .03],'Callback','timeBA = timeBA+1.0; uiresume;');
    pb10au = uicontrol(gcf,'Style','pushbutton','String','Shift---','units','normalized','Position',[.925 .30 .05 .03],'Callback','timeBA = timeBA-1.0; uiresume;');
    pb11 = uicontrol(gcf,'Style','pushbutton','String','ECG inverse','units','normalized','Position',[.925  .15 .05 .03],'Callback','Signal = -Signal; uiresume;');
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    uiwait(gcf)
    XL = xlim(ha);
    
end
close(gcf)

%% zscore all signals - do not zscore signals (Dr Serrador)

%Trim data down to just overlapped data
%Use letter o to indicate overlapped data
startbin = max([t(1) tHRV(1)]);
stopbin = min([t(end) tHRV(end)]);
int = find(tHR>=startbin & tHR<=stopbin);
tHRo = tHR(int);
HRo = HR(int);

%take clean data waveform data and create overlap variables
intwf1 = find(t>=startbin & t<=stopbin);
to = t(intwf1);
bp_fo = bp_f(intwf1);
%     for i=1:num_artery;
if exist('bv_f') == 1
    bv_fo = bv_f(intwf1,:);
    bv_fo_cm = BVcm_int(intwf1,:); %added by BS 2/16/18 in order to include tcd cm/s
else
    bv_f = BVcm_int;
    bv_fo = BVcm_int;
    bv_fo_cm = BVcm_int;
end
clear i

intwf2 = find(time>=startbin & time<=stopbin);
tECGo = time(intwf2);
ECGo = ECG;

intRRi_temp = find(iRRi>=intwf2(1) & iRRi<=intwf2(end));
iRRio = iRRi(intRRi_temp);
iRRio = iRRio - intwf2(1)+1;

if exist('Carotid_clean')
    intwf3 = find(ti>=startbin & ti<=stopbin);
    tio = ti(intwf1);
    Carotid_cleano = Carotid_clean(intwf1);
end

intwf4 = find(tco2>=startbin & tco2<=stopbin);
tco2o = tco2(intwf4);
CO2o = CO2(intwf4);

intbbb1 = find(tRRi>=startbin & tRRi<=stopbin);
tRRio = tRRi(intbbb1);
MBPo = MBP(intbbb1);
SBPo = SBP(intbbb1);
DBPo = DBP(intbbb1);
MBVo = MBV(intbbb1,:);
SBVo = SBV(intbbb1,:);
DBVo = DBV(intbbb1,:);
MBVcmo = MBVcm(intbbb1,:);
SBVcmo = SBVcm(intbbb1,:);
DBVcmo = DBVcm(intbbb1,:);

if exist('M_Carotid_Clean')
    M_Carotid_Cleano = M_Carotid_Clean(intbbb1);
    S_Carotid_Cleano = S_Carotid_Clean(intbbb1);
    D_Carotid_Cleano = D_Carotid_Clean(intbbb1);
end

intbbb2 = find(tetCO2>=startbin & tetCO2<=stopbin);
tetCO2o = tetCO2(intbbb2);
etCO2o = etCO2(intbbb2);

Velocity_Shift = Velocity_PeakDetect;

int2 = find(tHRV>=startbin & tHRV<=stopbin);
tHRVo = tHRV(int2);
MatVelo = Velocity_Shift(int2);

int5 = find(timeBA>=startbin & timeBA<=stopbin);
timeBAo = timeBA(int5);
DiamBAo = diamBA(int5);

int3 = find(tRRi>=startbin & tRRi<=stopbin);
tRRio = tRRi(int3);
iRRio = iRRi(int3);

tRRiV = timeMat(iRRiV);
int4 = find(tRRiV>=startbin & tRRiV<=stopbin);
iRRiVo = iRRiV(int4);
tRRiVo = tRRiV(int4);

endaligntime = max([max(tHRo) max(tHRVo) max(tRRio) max(tRRiVo) max(to) max(tco2o) max(tetCO2)]);

%% Plotting all the channels

figure(8)
clf
XL = [0 endaligntime];

% Set default figure size to near maximum of screen size
scrsz = get(0,'ScreenSize');
position_default = [0.01*scrsz(3) 0.07*scrsz(4) 0.98*scrsz(3) 0.85*scrsz(4)];
set(gcf,'position',position_default)

%slashes = findstr(infile,'\');
slashes = infile;
plot_title = [infile ' (Original Waveforms) analyzed  ' datestr(now)];

% Plot parameters
fs = 'fontsize';            % save fontsize field property
fs_title = 14;              % fontsize for titles
fs_label = 12;              % fontsize for axis labels
lw = 'LineWidth';           % save linewidth field property
lw_thin = 0.3;              % linewidth for thin lines
lw_thick = 2;               % linewidth for thicker lines
ms = 'markersize';          % save markersize field property
ms_large = 10;              % markersize for large points

if exist('Carotid_cleano')
    k =num_artery+6;
else
    k =num_artery+5;
end

button=0;

while button==0
    
    subplot(k,1,1)
    plot(to,bp_fo, 'b', lw, lw_thin), hold on;
    plot(tRRio,MBPo,'k', lw, lw_thick),
    plot(tRRio,SBPo,'r', lw, lw_thick),
    plot(tRRio,DBPo,'r', lw, lw_thick), hold off;
    axis tight
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    ylabel({'Arterial';'Blood Pressure';'(mmHg)'}, fs, fs_label)
    title(regexprep(plot_title, '_', ' '), fs, fs_title)
    
    for figcont = 1:num_artery
        subplot(k,1,figcont+1)
        plot(to,bv_fo(:,figcont), 'g', lw, lw_thin), hold on;
        plot(tRRio,MBVo(:,figcont),'k', lw, lw_thick),
        plot(tRRio,SBVo(:,figcont),'r', lw, lw_thick),
        plot(tRRio,DBVo(:,figcont),'r', lw, lw_thick),
        line([to(1) to(end)], [100 100], 'Color','m'),hold off;
        axis tight
        yl = ylim;
        yl(1) = yl(1)-yl(2)*.05;
        yl(2) = yl(2)+yl(2)*.05;
        ylim(yl);
        xlabel('Time (s)', fs, fs_label);
        ylabel({cell2mat(labels (figcont));'(%)'}, fs, fs_label);
    end
    
    subplot(k,1,num_artery+2)
    plot(tHRVo,MatVelo,'b', lw, lw_thick);
    axis tight
    yl(2) = max(MatVelo((round(length(MatVelo)/4)):((round(length(MatVelo)/4))*3)))*1.2;
    yl(1) = min(MatVelo((round(length(MatVelo)/4)):((round(length(MatVelo)/4))*3)))*0.8;
    ylim(yl);
    ylabel({'Carotid';'(cm/s)'}, fs, fs_label)
    
    subplot(k,1,num_artery+3)
    plot(timeBAo,DiamBAo,'r', lw, lw_thick);
    axis tight
    meanDiamBAo= mean(DiamBAo((round(length(DiamBAo)/4)):((round(length(DiamBAo)/4))*3)));
    stdDiamBAo = std(DiamBAo((round(length(DiamBAo)/4)):((round(length(DiamBAo)/4))*3)));
    yl(1) = meanDiamBAo-stdDiamBAo*4;
    yl(2) = meanDiamBAo+stdDiamBAo*4;
    ylim(yl); xlim([tHRVo(1) tHRVo(end)]);
    ylabel({'Diameter';'(cm)'}, fs, fs_label)
    
    subplot(k,1,num_artery+4)
    plot(tHRo,HRo,'r', lw, lw_thick);
    axis tight
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    ylabel({'HR';'(bpm)'}, fs, fs_label)
    
    subplot(k,1,num_artery+5)
    plot(tco2o,CO2o, 'g', lw, lw_thin), hold on;
    plot(tetCO2o,etCO2o,'r', lw, lw_thick), hold off;
    axis tight
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    ylabel({'etCO2';'(mmHg)'}, fs, fs_label)
    
    % this part to plot Carotid if signal recorded on chart
    if exist('Carotid_cleano')
        
        subplot(k,1,num_artery+6)
        plot(to,Carotid_cleao, 'g', lw, lw_thin), hold on;
        plot(tRRio,M_Carotid_Cleano,'k', lw, lw_thick),
        plot(tRRio,S_Carotid_Cleano,'r', lw, lw_thick),
        plot(tRRio,D_Carotid_Cleano,'r', lw, lw_thick),hold off;
        axis tight
        %    for txt=1:length(markers)
        %     tetx_handle = text (t(markers_Index(txt)) ,bv_f(markers_Index(txt),4),markers{txt});
        %     set(tetx_handle,'Colo','r','FontSize',txtfontsize,'Rotation',90)
        %    end
        yl = ylim;
        yl(1) = yl(1)-yl(2)*.05;
        yl(2) = yl(2)+yl(2)*.05;
        ylim(yl);
        xlabel('Time (s)', fs, fs_label);
        ylabel({cell2mat(labels (end));'(cm/s)'}, fs, fs_label);
        
    end
    
    
    % Resize figure for full page printing
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0.25 0.25 8 10.5]);
    
    pb = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','selblock = 1; button = 1; uiresume;');
    
    uiwait(gcf)
    
end
close(gcf)


%% Select region(time period) of interest

figure(9)
clf

start = min([min(tHRVo) min(timeBAo)]);
stop = max([max(tHRVo) max(timeBAo)]);
button = 0;
resety = 0;

while button==0
    clf
    
    subplot(2,1,1)
    plot(tHRVo,MatVelo)
    title('Velocity');
    axpos = get(gca,'position');
    xl = xlim;
    for i=1:1
        bbox(1) = (start(i)-xl(1))/(xl(2)-xl(1))*axpos(3)+axpos(1); %left position
        bbox(3) = ((stop(i)-xl(1))/(xl(2)-xl(1))-(start(i)-xl(1))/(xl(2)-xl(1)))*axpos(3); %width of box
        bbox(2) = axpos(2); %position from bottom
        bbox(4) = axpos(4); %height of box
        annotation('rectangle','Position',bbox,'facecolor',[0.5 0.5 0.5],'FaceAlpha',.3);
    end
    
    subplot(2,1,2)
    plot(timeBAo,DiamBAo)
    title('Diameter');
    axpos = get(gca,'position');
    xl = xlim;
    for i=1:1
        bbox(1) = (start(i)-xl(1))/(xl(2)-xl(1))*axpos(3)+axpos(1); %left position
        bbox(3) = ((stop(i)-xl(1))/(xl(2)-xl(1))-(start(i)-xl(1))/(xl(2)-xl(1)))*axpos(3); %width of box
        bbox(2) = axpos(2); %position from bottom
        bbox(4) = axpos(4); %height of box
        annotation('rectangle','Position',bbox,'facecolor',[0.5 0.5 0.5],'FaceAlpha',.3);
    end
    
    pb1u = uicontrol(gcf,'Style','pushbutton','String','Base+','units','normalized','Position',[.925 .9 .05 .05],'Callback','start = start+0.5; stop = stop+0.5; uiresume;');
    pb1l = uicontrol(gcf,'Style','pushbutton','String','Base-','units','normalized','Position',[.925 .8 .05 .05],'Callback','start = start-0.5; stop = stop-0.5;uiresume;');
    pb2u = uicontrol(gcf,'Style','pushbutton','String','Base++','units','normalized','Position',[.925 .7 .05 .05],'Callback','start = start+5; stop = stop+5; uiresume;');
    pb2l = uicontrol(gcf,'Style','pushbutton','String','Base--','units','normalized','Position',[.925 .6 .05 .05],'Callback','start = start-5; stop = stop-5; uiresume;');
    pb3u = uicontrol(gcf,'Style','pushbutton','String','BaseL+','units','normalized','Position',[.925 .5 .05 .05],'Callback','stop = stop+0.5; uiresume;');
    pb3l = uicontrol(gcf,'Style','pushbutton','String','BaseL-','units','normalized','Position',[.925 .4 .05 .05],'Callback','stop = stop-0.5; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','BaseL++','units','normalized','Position',[.925 .3 .05 .05],'Callback','stop = stop+5; uiresume;');
    pb4l = uicontrol(gcf,'Style','pushbutton','String','BaseL--','units','normalized','Position',[.925 .2 .05 .05],'Callback','stop = stop-5; uiresume;');
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    
    uiwait(gcf)
    
end
close(gcf)

block_data_start = start;
block_data_stop = stop;
startime = start(1);
endtime = stop(1);

%% save old interpolated data _old variable and replace with new overlapped
%interopolated data

BPi_old = BPi;
BVi_old = BVi;
CO2i_old = CO2i;
clear BPi BVi CO2i

%resample everything to 100 Hz for cleaning and further analysis
tint = startime:0.01:endtime;
if length(tECGo)~=length(ECGo) %added if statement since if data was collected at higher sampling caused error
    ECGo=interp1(t,ECGo,tECGo,'linear');
end

ECGi = interp1(tECGo,ECGo,tint,'linear');
BPi = interp1(to,bp_fo,tint,'linear');
for i=1:num_artery
    BVi(:,i) = interp1(to.',bv_fo(:,i),tint.','linear');
end

MatVeli = interp1(tHRVo,MatVelo,tint,'linear');
DiamBAi = interp1(timeBAo,DiamBAo,tint,'linear');
CO2i = interp1(tco2o,CO2o,tint,'linear');

%find times of ECG beats that fall within interpolated area
good = find(time(iRRio)>=startime & time(iRRio)<=endtime);
tECGPeak = time(iRRio(good));

%find all the beat by beat values that are within overlap range
intbbb1 = find(tRRio>=startime & tRRio<=endtime);
tRRio = tRRio(intbbb1);
iRRio = iRRio(intbbb1);
MBPo = MBPo(intbbb1);
SBPo = SBPo(intbbb1);
DBPo = DBPo(intbbb1);
MBVo = MBVo(intbbb1,:);
SBVo = SBVo(intbbb1,:);
DBVo = DBVo(intbbb1,:);

if exist('M_Carotid_Clean')
    M_Carotid_Cleano = M_Carotid_Cleano(intbbb1);
    S_Carotid_Cleano = S_Carotid_Cleano(intbbb1);
    D_Carotid_Cleano = D_Carotid_Cleano(intbbb1);
end

%% start cleaning of velocity and diameter

% Plot parameters
warning('off');
data_ch =data;
fs = 'fontsize';            % save fontsize field property
fs_title = 14;              % fontsize for titles
fs_label = 12;              % fontsize for axis labels
lw = 'LineWidth';           % save linewidth field property
lw_thin = 0.3;              % linewidth for thin lines
lw_thick = 2;               % linewidth for thicker lines
ms = 'markersize';          % save markersize field property
ms_large = 10;              % markersize for large points

CreateStruct.WindowStyle = 'modal';
CreateStruct.Interpreter = 'tex';
AddOpts.Resize = 'on';
AddOpts.WindowStyle = 'modal';
AddOpts.Interpreter = 'tex';
Prompt = {'Output File Name','BP Channel','ICH side Channel','Contralteral Channel','ECG Channel'};
Title = 'Velocity/Diameter Cleaning';
LineNo = 1;
Question = 'Output File Already Exists, Overwrite?';
Qtitle = 'Transfer Function Analysis';

% check to see if you need to low pass filter the velocity signal
button = 0;
LPfilter = 0;
BVfilt = MatVeli;
BVfiltorig = BVfilt;
MatVeli_orig = MatVeli;
figure(10)
n=1;
while button==0
    
    clf
    plot(BVfilt);
    title('Lowpass filter velocity or not');
    
    pb1u = uicontrol(gcf,'Style','pushbutton','String','LPF-3','units','normalized','Position',[.925 .9 .05 .03],'Callback','n = 3; LPfilter = 1; redetect = 1; uiresume;');
    pb2u = uicontrol(gcf,'Style','pushbutton','String','LPF-7','units','normalized','Position',[.925 .8 .05 .03],'Callback','n = 7; LPfilter = 1; redetect = 1; uiresume;');
    pb3u = uicontrol(gcf,'Style','pushbutton','String','LPF-11','units','normalized','Position',[.925 .7 .05 .03],'Callback','n = 11; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-21','units','normalized','Position',[.925 .6 .05 .03],'Callback','n = 21; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-33','units','normalized','Position',[.925 .5 .05 .03],'Callback','n = 33; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-66','units','normalized','Position',[.925 .4 .05 .03],'Callback','n = 66; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-99','units','normalized','Position',[.925 .3 .05 .03],'Callback','n = 99; LPfilter = 1; redetect = 1; uiresume;');
    pb1du = uicontrol(gcf,'Style','pushbutton','String','Undo','units','normalized','Position',[.925 .2 .05 .03],'Callback','BVfilt = BVfiltorig; uiresume;');
    
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    
    uiwait(gcf)
    
    if LPfilter == 1
        BVf = medfilt1(BVfilt,n);
        BVfilt = BVf;
        LPfilter = 0;
    end
    
    
end
MatVeli = BVfilt;


% check to see if you need to low pass filter the diameter signal
button = 0;
LPfilter = 0;
BVfilt = DiamBAi;
BVfiltorig = BVfilt;
MatVeli_orig = DiamBAi;
n=1;
while button==0
    
    clf
    plot(BVfilt);
    title('Lowpass filter diameter or not');
    pb1u = uicontrol(gcf,'Style','pushbutton','String','LPF-3','units','normalized','Position',[.925 .9 .05 .03],'Callback','n = 3; LPfilter = 1; redetect = 1; uiresume;');
    pb2u = uicontrol(gcf,'Style','pushbutton','String','LPF-7','units','normalized','Position',[.925 .8 .05 .03],'Callback','n = 7; LPfilter = 1; redetect = 1; uiresume;');
    pb3u = uicontrol(gcf,'Style','pushbutton','String','LPF-11','units','normalized','Position',[.925 .7 .05 .03],'Callback','n = 11; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-21','units','normalized','Position',[.925 .6 .05 .03],'Callback','n = 21; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-33','units','normalized','Position',[.925 .5 .05 .03],'Callback','n = 33; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-66','units','normalized','Position',[.925 .4 .05 .03],'Callback','n = 66; LPfilter = 1; redetect = 1; uiresume;');
    pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF-99','units','normalized','Position',[.925 .3 .05 .03],'Callback','n = 99; LPfilter = 1; redetect = 1; uiresume;');
    pb1du = uicontrol(gcf,'Style','pushbutton','String','Undo','units','normalized','Position',[.925 .2 .05 .03],'Callback','BVfilt = BVfiltorig; uiresume;');
    
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    
    uiwait(gcf)
    
    if LPfilter == 1
        BVf = medfilt1(BVfilt,n);
        BVfilt = BVf;
        LPfilter = 0;
    end
    
    
end
DiamBAi = BVfilt;
close(gcf)

%% calculate beat by beat data

%determine number of beats to be analayzed
w = length(tECGPeak);

%create arrays for compiling average waveforms of velocity and diameter
%making it 200 points long (2 seconds) since will be longer than any single
%beat and thus won't have issues with assigning values later
MatVeliPlot = NaN(w-3,200);
DiamBAiPlot = NaN(w-3,200);

for i = 1:w(1)-1
    
    %cacluate velocity based on clean data times of ECG peaks
    intbeat = find(tint>=tECGPeak(i) & tint<tECGPeak(i+1));
    MVelo(i) = nanmean(MatVeli(intbeat));
    [SVelo(i) SVeloPt(i)] = nanmax(MatVeli(intbeat(1:30)));  % calculate systolic value and point
    [DVelo(i) DVeloPt(i)] = nanmin(MatVeli(intbeat(end-30:end))); %find diastolic point between systolic point calculated plus beat length
    SVeloPt(i) = SVeloPt(i) + intbeat(1)-1;
    DVeloPt(i) = DVeloPt(i) + intbeat(1) + length(intbeat) - 30 -1;
    
    MDiamo(i) = nanmean(DiamBAi(intbeat));
    [SDiamo(i) SDiamPto(i)] = nanmax(DiamBAi(intbeat(:)));
    SDiamPto(i) = SDiamPto(i) + intbeat(1)-1;
    [DDiamo(i) DDiamPto(i)] = nanmin(DiamBAi(intbeat(end-30:end)));
    DDiamPto(i) = DDiamPto(i) + intbeat(1) + length(intbeat) - 30 -1;
    
    if intbeat(1)<31
        startptplot = intbeat(1);
    else
        startptplot = intbeat(1)-30;
    end
    
    if (intbeat(end)+50)>length(tint)
        endptplot = length(tint);
    else
        endptplot = intbeat(end)+50;
    end
    
    %plot individual beats
    figure(11)
    subplot(2,2,1)
    plot(tint(startptplot:endptplot),MatVeli(startptplot:endptplot))
    hold on
    plot(tint(SVeloPt(i)),MatVeli(SVeloPt(i)),'go')
    plot(tint(DVeloPt(i)),MatVeli(DVeloPt(i)),'ro')
    hold off
    xlim([tint(startptplot) tint(endptplot)])
    title('Velocity');
    
    subplot(2,2,3)
    plot(tint(startptplot:endptplot),DiamBAi(startptplot:endptplot))
    hold on
    plot(tint(SDiamPto(i)),DiamBAi(SDiamPto(i)),'go')
    plot(tint(DDiamPto(i)),DiamBAi(DDiamPto(i)),'ro')
    hold off
    xlim([tint(startptplot) tint(endptplot)])
    title('Diameter');
    
    if i>1 && i<(w-1)
        
        %starting with second beat save each individual value
        MatVeliPlot(i,1:length(MatVeli((intbeat(1)-30):(intbeat(end)+50)))) = MatVeli((intbeat(1)-30):(intbeat(end)+50));
        DiamBAiPlot(i,1:length(DiamBAi((intbeat(1)-30):(intbeat(end)+50)))) = DiamBAi((intbeat(1)-30):(intbeat(end)+50));
        ECGiPlot(i,1:length(ECGi((intbeat(1)-30):(intbeat(end)+50)))) = ECGi((intbeat(1)-30):(intbeat(end)+50));
        BPiPlot(i,1:length(BPi((intbeat(1)-30):(intbeat(end)+50)))) = BPi((intbeat(1)-30):(intbeat(end)+50));
        if num_artery>0
            if min(size(BVi))==3
                BV1Plot(i,1:length(BVi((intbeat(1)-30):(intbeat(end)+50),1))) = BVi((intbeat(1)-30):(intbeat(end)+50),1);
                BV2Plot(i,1:length(BVi((intbeat(1)-30):(intbeat(end)+50),2))) = BVi((intbeat(1)-30):(intbeat(end)+50),2);
                BV3Plot(i,1:length(BVi((intbeat(1)-30):(intbeat(end)+50),3))) = BVi((intbeat(1)-30):(intbeat(end)+50),3);
            end
            
            if min(size(BVi))==2
                BV1Plot(i,1:length(BVi((intbeat(1)-30):(intbeat(end)+50),1))) = BVi((intbeat(1)-30):(intbeat(end)+50),1);
                BV2Plot(i,1:length(BVi((intbeat(1)-30):(intbeat(end)+50),2))) = BVi((intbeat(1)-30):(intbeat(end)+50),2);
            end
            
            if min(size(BVi))==1
                BV1Plot(i,1:length(BVi((intbeat(1)-30):(intbeat(end)+50),1))) = BVi((intbeat(1)-30):(intbeat(end)+50),1);
            end
        end
        subplot(2,2,2)
        hold on
        plot(MatVeliPlot(i,:))
        plot((SVeloPt(i)-(intbeat(1)-30)),MatVeliPlot(i-1,(SVeloPt(i)-(intbeat(1)-30))),'go')
        plot((DVeloPt(i)-(intbeat(1)-30)),MatVeliPlot(i-1,(DVeloPt(i)-(intbeat(1)-30))),'ro')
        title('Velocity beat by beat');
        
        subplot(2,2,4)
        hold on
        plot(DiamBAiPlot(i,:))
        plot((SDiamPto(i)-(intbeat(1)-30)),DiamBAiPlot(i-1,(SDiamPto(i)-(intbeat(1)-30))),'go')
        plot((DDiamPto(i)-(intbeat(1)-30)),DiamBAiPlot(i-1,(DDiamPto(i)-(intbeat(1)-30))),'ro')
        title('Diameter beat by beat');
    end
    
    
    if exist('BV1')==1
        [LIA,LOCB] = ismember(SBV(i,1),BV1(:,1));
        NDX_Max(i,m) = LOCB;
        
        [LIA,LOCB] = ismember(DBV(i,1),BV1(:,1));
        NDX_Min(i,m) = LOCB;
    end
    
end

MatVeliPlotMedian = nanmedian(MatVeliPlot);
DiamBAiPlotMedian = nanmedian(DiamBAiPlot);
SingleBeatFlow = MatVeliPlotMedian.*(pi()*(DiamBAiPlotMedian/2).^2).*60;

figure(12)
subplot(3,1,1)
plot(MatVeliPlotMedian(1:150))
title('Velocity Median');
axis tight
subplot(3,1,2)
plot(DiamBAiPlotMedian(1:150))
title('Diameter Median');
axis tight
subplot(3,1,3)
plot(SingleBeatFlow(1:150))
axis tight
title('Blood Flow');

%% make beat by beat data from cleaned data same as velocity and diameter
ikeep = 1:length(MVelo);
tRRio = tRRio(ikeep);
MBPo = MBPo(ikeep);
SBPo = SBPo(ikeep);
DBPo = DBPo(ikeep);
MBVo = MBVo(ikeep,:);
SBVo = SBVo(ikeep,:);
DBVo = DBVo(ikeep,:);
MBVcmon = MBVcmo(ikeep,:);
SBVcmon = SBVcmo(ikeep,:);
DBVcmon = DBVcmo(ikeep,:);
if exist('M_Carotid_Clean')
    M_Carotid_Cleano = M_Carotid_Cleano(ikeep);
    S_Carotid_Cleano = S_Carotid_Cleano(ikeep);
    D_Carotid_Cleano = D_Carotid_Cleano(ikeep);
end

%% Plot Velocity,and Diameter Sys, Dias,and mean

figure(13)

button = 0;
BVel1(:,1) = MatVeli;
BVel1(:,2) = DiamBAi;

tVelBeat = tint(SVeloPt);
SBVel(:,1) = SVelo;
MBVel(:,1) = MVelo;
DBVel(:,1) = DVelo;

tDiamBeat = tint(SDiamPto);
SBVel(:,2) = SDiamo;
MBVel(:,2) = MDiamo;
DBVel(:,2) = DDiamo;

% question what is tRRIO?
% tVDBeat(:,1) = tRRio;
% tVDBeat(:,2) = tRRio;
% or should it be:
tVDBeat(:,1) = tVelBeat;
tVDBeat(:,2) = tDiamBeat;

title_name = ['Velocity', 'Diameter'];
while button == 0
    
    for nm=1:2
        subplot(2,1,nm)
        hold on
        plot(tint,BVel1(:,nm), 'b', lw, lw_thin)
        plot(tVDBeat(:,nm),SBVel(:,nm),'ro', lw, lw_thick)
        plot(tVDBeat(:,nm),MBVel(:,nm),'ko', lw, lw_thick)
        plot(tVDBeat(:,nm),DBVel(:,nm),'go', lw, lw_thick)
        hold off
        axis tight
        yl = ylim;
        yl(1) = yl(1)-yl(2)*.05;
        yl(2) = yl(2)+yl(2)*.05;
        ylim(yl);
        title(title_name(nm));
        xlabel('Time (s)', fs, fs_label)
    end
    
    b4u = uicontrol(gcf,'Style','pushbutton','String',' Continue ','BackgroundColor','g','FontSize',15,'units','normalized','Position',[.025 .05 .1 .05],'Callback', 'button = 1 ; uiresume;');
    uiwait(gcf)
end
close(gcf)

%% Drop Filter Procedure
% Preallocate vector for removed beats

tVel_bad= nan(length(tRRio),2);
MBVel_tempf_bad = nan(length(tRRio),2);
SBVel_tempf_bad = nan(length(tRRio),2);
DBVel_tempf_bad = nan(length(tRRio),2);

for m = 1:2
    
    figure(14)
    scrsz = get(0,'ScreenSize');
    set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
    clf
    
    plot(tint,BVel1(:,m), 'b', lw, lw_thin), hold on;
    plot(tVDBeat(:,nm),MBVel(:,m),'ko', lw, lw_thick),
    plot(tVDBeat(:,nm),SBVel(:,m),'ro', lw, lw_thick),
    plot(tVDBeat(:,nm),DBVel(:,m),'go', lw, lw_thick), hold off;
    axis tight
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    xlabel('Time (s)', fs, fs_label)
    
    button = questdlg('Do you need to interpolate for bad beats?','Continue Operation','Yes','No','Help','No');
    if strcmp(button,'Yes')
        disp('Performing drop filter')
    elseif strcmp(button,'No')
        disp('Continue')
    elseif strcmp(button,'Help')
        disp('Sorry, no help available')
        break
    end
    
    if strcmp(button,'Yes')
        
        MBV_temp = MBVel(:,m);
        SBV_temp = SBVel(:,m);
        DBV_temp = DBVel(:,m);
        tVDBeat_temp = tVDBeat(:,m);
        
        bv_f_temp = BVel1(:,m);
        
        %Remove servocorrects from bp signal
        PP = (SBV_temp-DBV_temp)./MBV_temp;
        stdPP = nanstd(PP);
        meanPP = nanmean(PP);
        LPPthreshold=meanPP-2*stdPP;
        UPPthreshold=meanPP+2*stdPP;
        good = find(PP>=LPPthreshold & PP<=UPPthreshold);
        
        button = 0;
        Remove_Group=0;
        
        XL = [min(time(1:length(bv_f_temp))) max(time(1:length(bv_f_temp)))];
        XL2 = XL;
        
        while button == 0
            
            instruct1 = ['Remove physiocals by looking for beats with reduced pulse pressure (PP). Adjust upper threshold (red line) with UThresh buttons (+ increases a little, ++ increases a lot, - decreases a little, -- decreases a lot). Adjust lower threshold (green) using the LThresh buttons. Once you are happy, press Continue'];
            
            clf
            good = find(PP>=LPPthreshold & PP<=UPPthreshold);
            
            clf
            subplot(2,1,1)
            plot(tint,bv_f_temp, 'b', lw, lw_thin), hold on;
            plot(tVDBeat(good),MBV_temp(good),'ko', lw, lw_thick),
            plot(tVDBeat(good),SBV_temp(good),'ro', lw, lw_thick),
            plot(tVDBeat(good),DBV_temp(good),'ro', lw, lw_thick), hold off;
            yl(1) = yl(1)-yl(2)*.05;
            yl(2) = yl(2)+yl(2)*.05;
            ylim(yl);
            ylabel('Velocity', fs, fs_label)
            xlabel('Time (s)', fs, fs_label)
            xlim(XL)
            
            subplot(2,1,2)
            plot(tVDBeat(1:length(PP),m),PP, 'bo')
            hold on
            hlineU = refline(0,LPPthreshold);
            set(hlineU,'Color','r')
            hlineL = refline(0,UPPthreshold);
            set(hlineL,'Color','g')
            hold off
            ylim;
            xlim(XL2)
            ylabel('Pulse Pressure  (mmHg)', fs, fs_label)
            
            pb1u = uicontrol(gcf,'Style','pushbutton','String','UThresh+','units','normalized','Position',[.925 .9 .05 .05],'Callback', ...
                'UPPthreshold=UPPthreshold+meanPP*0.1; uiresume;');
            pb2u = uicontrol(gcf,'Style','pushbutton','String','UThresh-','units','normalized','Position',[.925 .8 .05 .05],'Callback', ...
                'UPPthreshold=UPPthreshold-meanPP*0.1; uiresume;');
            pb3u = uicontrol(gcf,'Style','pushbutton','String','UThresh++','units','normalized','Position',[.925 .7 .05 .05],'Callback', ...
                'UPPthreshold=UPPthreshold+meanPP*0.5; uiresume;');
            pb4u = uicontrol(gcf,'Style','pushbutton','String','UThresh--','units','normalized','Position',[.925 .6 .05 .05],'Callback', ...
                'UPPthreshold=UPPthreshold-meanPP*0.5; uiresume;');
            pb5u = uicontrol(gcf,'Style','pushbutton','String','LThresh+','units','normalized','Position',[.925 .5 .05 .05],'Callback', ...
                'LPPthreshold=LPPthreshold+meanPP*0.1; uiresume;');
            pb6u = uicontrol(gcf,'Style','pushbutton','String','LThresh-','units','normalized','Position',[.925 .4 .05 .05],'Callback', ...
                'LPPthreshold=LPPthreshold-meanPP*0.1; uiresume;');
            pb7u = uicontrol(gcf,'Style','pushbutton','String','LThresh++','units','normalized','Position',[.925 .3 .05 .05],'Callback', ...
                'LPPthreshold=LPPthreshold+meanPP*0.5; uiresume;');
            pb8u = uicontrol(gcf,'Style','pushbutton','String','LThresh--','units','normalized','Position',[.925 .2 .05 .05],'Callback', ...
                'LPPthreshold=LPPthreshold-meanPP*0.5; uiresume;');
            pb7l = uicontrol(gcf,'Style','pushbutton','String','Instructions','units','normalized','Position',[.925 .125 .05 .05],'Callback','questdlg(instruct1,''Instructions'',''Continue'',''Continue''); uiresume;');
            pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
            
            uiwait(gcf)
        end
        
        MBV_temp_old = MBV_temp;
        SBV_temp_old = SBV_temp;
        DBV_temp_old = DBV_temp;
        
        tbp = tVDBeat_temp(good);
        
        MBV_tempf = MBV_temp(good);
        SBV_tempf = SBV_temp(good);
        DBV_tempf = DBV_temp(good);
        
        
        button = 0;
        bad = [];
        
        PPf = (SBV_tempf-DBV_tempf)./MBV_tempf;
        
        XL = [min(time(1:length(bv_f_temp))) max(time(1:length(bv_f_temp)))];
        XL2 = XL;
        
        while button == 0
            
            instruct1 = ['Remove any bad beats. Use remove button and place cross hairs over beat that is not correct. Once you are happy, press Continue'];
            
            clf
            
            subplot(3,1,1:2)
            plot(tint,bv_f_temp, 'b', lw, lw_thin), hold on;
            plot(tbp,MBV_tempf,'ko', lw, lw_thick),
            plot(tbp,SBV_tempf,'ro', lw, lw_thick),
            plot(tbp,DBV_tempf,'ro', lw, lw_thick), hold off;
            xlim(XL)
            ha = get(gcf,'CurrentAxes');
            subplot(3,1,3)
            plot(tbp,PPf,'ro-')
            ylabel('HR')
            xlim(XL2)
            ha2 = get(gcf,'CurrentAxes');
            ss = get(gcf,'Position');
            
            pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove','units','normalized','Position',[.925 .8 .05 .05],'Callback', ...
                '[x,y,butnum] = ginput(1); [mpt bad] = min(abs(tbp-x)); tbpold = tbp; MBV_tempfold = MBV_tempf; SBV_tempfold = SBV_tempf; DBV_tempfold = DBV_tempf; PPfold = PPf; tbp(bad) = []; MBV_tempf(bad) = []; SBV_tempf(bad) = []; DBV_tempf(bad) = []; PPf(bad)=[]; uiresume;');
            pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','units','normalized','Position',[.925 .6 .05 .05],'Callback', ...
                '[x, y] = ginput(1); x = x(end); [minpt addpt] = min(abs(tRRio-x)); tbpold = tbp; MBV_tempfold = MBV_tempf; SBV_tempfold = SBV_tempf; DBV_tempfold = DBV_tempf; PPfold = PPf; tbp(end+1) = tRRio(addpt); MBV_tempf(end+1) = MBV_temp(addpt); SBV_tempf(end+1) = SBV_temp(addpt); DBV_tempf(end+1) = DBV_temp(addpt); [tbp ix] = sort(tbp); MBV_tempf = MBV_tempf(ix); SBV_tempf = SBV_tempf(ix); DBV_tempf = DBV_tempf(ix);  PPf = SBV_tempf-DBV_tempf; uiresume;');
            
            pb3u = uicontrol(gcf,'Style','pushbutton','String','Undo','units','normalized','Position',[.925 .4 .05 .05],'Callback', ...
                'tbp = tbpold; MBV_tempf = MBV_tempfold; SBV_tempf = SBV_tempfold; DBV_tempf = DBV_tempfold; PPf = PPfold; uiresume;');
            pb7l = uicontrol(gcf,'Style','pushbutton','String','Instructions','units','normalized','Position',[.925 .2 .05 .05],'Callback','questdlg(instruct1,''Instructions'',''Continue'',''Continue''); uiresume;');
            pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .1 .05 .05],'Callback','button = 1; uiresume;');
            pbRG1u = uicontrol(gcf,'Style','pushbutton','String','Remove Group','Position',[ss(3)*.025 ss(4)*.7 ss(3)*.1 ss(4)*.05],'Callback','Remove_Group = 1; uiresume;');
            
            uiwait(gcf)
            
            if Remove_Group > 0
                
                k = waitforbuttonpress;
                point1 = get(gca,'CurrentPoint');    % button down detected
                finalRect = rbbox;                   % return figure units
                point2 = get(gca,'CurrentPoint');    % button up detected
                point1 = point1(1,1:2);              % extract x and y
                point2 = point2(1,1:2);
                p1 = min(point1,point2);             % calculate locations
                offset = abs(point1-point2);         % and dimensions
                x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)]; % is the time
                y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
                hold on
                axis manual
                fig = plot(x,y,'r','linewidth',2) ;
                
                Index_Start_Array=find(tbp<= min(x));
                Index_end_Array  = find(tbp<= max(x));
                
                if  (isempty(Index_Start_Array)) ||  (isempty(Index_end_Array) )
                    
                    display('');
                    
                else
                    Index_Start = Index_Start_Array(end);
                    Index_end = Index_end_Array(end);
                    
                    tbpold = tbp; MBV_tempfold = MBV_tempf; SBV_tempfold = SBV_tempf; DBV_tempfold = DBV_tempf; PPfold = PPf;
                    tbp(Index_Start:Index_end) = [];
                    MBV_tempf(Index_Start:Index_end) = []; SBV_tempf(Index_Start:Index_end) = []; DBV_tempf(Index_Start:Index_end) = []; PPf(Index_Start:Index_end)=[];
                    
                    Remove_Group = 0;
                end
            end
            
            XL = xlim(ha);
            XL2 = xlim(ha2);
            
        end
        
        MBVelold(:,m) = MBVel(:,m);
        SBVelold(:,m) = SBVel(:,m);
        DBVelold(:,m) = DBVel(:,m);
        
        %Determine the indeces and values for the good beats (beats kept)
        MBVel(:,m) = interp1(tbp,MBV_tempf,tRRio,'linear');
        SBVel(:,m) = interp1(tbp,SBV_tempf,tRRio,'linear');
        DBVel(:,m) = interp1(tbp,DBV_tempf,tRRio,'linear');
        
        tRRis = tint(SVeloPt);
        tRRid =  tint(DVeloPt);
        
        %Determine the indeces and values for the bad beats (beats removed)
        inds_bp = find(~ismember( tRRio, tbp));
        tbp_bad(inds_bp,m) = tRRio(inds_bp);
        MBVel_tempf_bad(inds_bp,m) = MBVelold(inds_bp,m);
        SBVel_tempf_bad(inds_bp,m) = SBVelold(inds_bp,m);
        DBVel_tempf_bad(inds_bp,m) = DBVelold(inds_bp,m);
    end
    close(gcf)
end


%% Calculate Absolute Flow Based on mean velocity and diameter values
BpFlow = BPi;
CBV = BVel1(:,1);
DiamFlow = BVel1(:,2);

%Brain size calculations
% prompt = {'Subjects Weight in Pounds(lbs):','Subjects Height in
% inches(in):','Subjects Sex(M/F):'}; % May be used later if necessary
prompt = {'Subjects Age(Yrs):','Subjects Sex (M/F):'};

dlg_title = 'Input';
num_lines = 1;
def = {'','',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
% BrainSize = answer(1);
% Height = answer(2); %Not sure what we need to do with this variable
Age = str2double (answer(1));
Sex = cell2mat(answer(2));

% A= cell2mat(BrainSize);
% A= str2double (A);

% Brain size calculations are based on modified Birch and Riddle equations
% (Borzage et al, 2012 "Equations to describe brain size across the
% continuum of human lifespan")
% BrainWeight = 6*(0.12*((A * 453.59237) ).^(2/3));%Brain mass in g but
% this is wrong

%The piecewise function for male brain weight (MBW), in Kilograms, is as follows
%for age N > 1.58;MBW(age) = 1:134 × age^0.1068 × exp (? age/252.3) (Riddle
%et al, 2010 "Modeling brain tissue volumes over the lifespan: quantitative analysis of
%postmortem weights and in vivo MR images")
%P1 = 1.134; P2 = 0.1068; P3 = 252.3

%The piecewise function for Female brain weight (FBW), in Kilograms, is as follows
%for age N > 1.18;FBW(age) = 0.955 × age^0.1357 × exp (? age/212.1) (Riddle
%et al, 2010 "Modeling brain tissue volumes over the lifespan: quantitative analysis of
%postmortem weights and in vivo MR images")
%P1 = 1.18; P2 = 0.1357; P3 = 212.1

% dy1/dx = (P1 * y1(P2-y1))/(P2-y1 + P3*Y1) % Original Birch equation
% y = y1*(1+P4*x + P5*x^2); % x = age, y = brain size, P1-5 are parameters, y1 = brain size for the original unmodified Birch
% BrainWeight = y;

% Set the Riddle model parameter values based on Sex input responses
if strcmpi(Sex,'M')
    P1 = 399;
    P2 = 594;
    P3 = 1030;
    P4 = 0.147;
    P5 = 194;
elseif strcmpi(Sex,'F')
    P1 = 400;
    P2 = 584;
    P3 = 1010;
    P4 = 0.116;
    P5 = 229;
end

%Total cerebral weight (grams)
BrainWeight = P3*Age^(P4)*exp(-Age/P5);

%Mean Carotid Diameter Calculations
DiamFlow_cm = DiamFlow;  %DiamAvg is in cm

%Cerebral Vascular Resistance calculations
%http://www.vhlab.umn.edu/atlas/physiology-tutorial/blood-flow.shtml
CBF = pi()*(((DiamFlow_cm)/2).^2).* CBV * 60;

SBD_flow = pi().*((SBVel(:,2)./2).^2).* SBVel(:,1).* 60;
MBD_flow = pi().*((MBVel(:,2)./2).^2).* MBVel(:,1).* 60;
DBD_flow = pi().*((DBVel(:,2)./2).^2).* DBVel(:,1).* 60;

CBF_mean = nanmean(MBD_flow);
%MBD = MBV(:,3);
Diam_mean = nanmean(MBVel(:,2));

%Cerebreal Blood Flow per hundred gram of brain tissue is in mL/g.min
CBF_Hgbt = (CBF./ BrainWeight)* 100;
% mean_CBF = nanmean(CBF);
% mean_CBF_Hgbt = nanmean(CBF_Hgbt);

%Calculate the Systolic, mean, and diastolic per hundred gram of brain tissue is in mL/g.min
SBD_flow_Hgbt = (SBD_flow./ BrainWeight)* 100;
MBD_flow_Hgbt = (MBD_flow./ BrainWeight)* 100;
DBD_flow_Hgbt = (DBD_flow./ BrainWeight)* 100;

%MAP = MBV(:,1);

CVR = MBPo./ MBD_flow;
% CVR_mean = nanmean(CVR);
CVR_mean = nanmean(MBPo)/mean(MBD_flow); %changed --> took means first then divided

%% Distensibility Calculations (%)
% Assessment of flow-mediated dilation in humans: a methodological and physiological guideline
% DHJ Thijssen, MA Black, KE Pyke - American Journal of , 2011 - Am Physiological Soc
% Arterial distensibility** Relative change in diameter (or area) for a given pressure change;
% inverse of elastic modulus
% DD/(DP3D) (mmHg 1)

% J Ultrasound Med. Author manuscript; available in PMC 2009 May 5.
% 1. Strain as the amount of deformation relative to the unstressed state and expressed as
% percent change in the arterial diameter: strain = (SD  DD)/DD, where SD was the
% systolic and DD the diastolic CCA diameter;
% 2. Stiffness (?) as stress (SBP  DBP)-to-strain ratio: ? = ln(SBP/DBP)/strain, where
% SBP and DBP were brachial BPs measured in the systolic and diastolic cardiac cycles,
% respectively;
% 3. Distensibility as 1/? and adjusted to IMT: 1/?, = 1/[ln(SBP/DBP)/strain × IMT];
%    Also calculated as Distensibility = (2*(Ds-Dd) / Dd)/(SBP-DBP)
% 4. Pressure-strain Youngs elastic modulus (E): E = K (SBP  DBP)/strain, where K =
% 133.3 was the conversion factor for mm Hg to Nm?2

%first calculate strain
Ds =SBVel(:,2).*10; %systolic diameter converted to mm
Dd =DBVel(:,2).*10; %diastolic diameter converted to mm
Strain = (Ds - Dd)./Dd;
Stiffness = log(SBPo./DBPo)./Strain;
Distensibility_new = 1/Stiffness;
Distensibility = ((2*(Ds - Dd))./Dd) ./ (SBPo-DBPo);

mean_Distensibility = nanmean(Distensibility);
Stiffness_mean = nanmean(Stiffness);

% Mean values for BP
Sys_BP_Mean = nanmean(SBPo);
Mean_BP_Mean = nanmean(MBPo);
Dias_BP_Mean = nanmean(DBPo);
% Mean values for Velocity
Sys_Vel_Mean = nanmean(SBVel(:,1));
Mean_Vel_Mean = nanmean(MBVel(:,1));
Dias_Vel_Mean = nanmean(DBVel(:,1));
% Mean values for Diameter
Sys_Diam_Mean = nanmean(SBVel(:,2));
Mean_Diam_Mean = nanmean(MBVel(:,2));
Dias_Diam_Mean = nanmean(DBVel(:,2));
% Mean values for CBF
Sys_CBF_Mean = nanmean(SBD_flow);
Mean_CBF_Mean = nanmean(MBD_flow);
Dias_CBF_Mean = nanmean(DBD_flow);
% Mean values for CBF_Hgbt
Sys_CBF_Hgbt_Mean = nanmean(SBD_flow_Hgbt);
Mean_CBF_Hgbt_Mean = nanmean(MBD_flow_Hgbt);
Dias_CBF_Hgbt_Mean = nanmean(DBD_flow_Hgbt);

%% Plot of Diameter, Velocity, and Flow

button = 0;
figure(15)
clf('reset')
while button == 0
    hold on;
    % Plot of BP
    subplot(3,2,1)
    plot(tint,BpFlow,'b');hold on;
    plot(tRRio,SBPo,'ro');
    plot(tRRio,MBPo,'yo');
    plot(tRRio,DBPo,'go');
    hold off
    axis tight
    ylim([min(DBPo)*.9 max(SBPo)*1.1])
    title('Blood Pressure')
    ylabel('mmHg')
    %scrollplot
    
    % Plot of Diameter
    subplot(3,2,2)
    plot(tint,DiamFlow_cm,'b');hold on;
    plot(tRRio,SBVel(:,2),'ro');
    plot(tRRio,MBVel(:,2),'yo');
    plot(tRRio,DBVel(:,2),'go');
    hold off
    axis tight
    ylim([min(DBVel(:,2))*.9 max(SBVel(:,2))*1.1])
    title('Diameter')
    ylabel('cm')
    %scrollplot
    
    % Plot of Cerebral Velocity
    subplot(3,2,3)
    plot(tint,CBV, lw, lw_thin),
    hold on
    plot(tRRio,SBVel(:,1),'ro');
    plot(tRRio,MBVel(:,1),'yo');
    plot(tRRio,DBVel(:,1),'go');
    hold off
    axis tight
    ylim([min(DBVel(:,1))*.9 max(SBVel(:,1))*1.1])
    %     ylim([min(DBV)/10*.9 max(SBV)/10*1.1])
    title('Carotid Velocity')
    ylabel('cm/sec')
    %scrollplot
    
    % Plot of Cerebral Blood Flow
    subplot(3,2,4)
    plot(tint,CBF)
    hold on
    plot(tRRio,SBD_flow,'ro');
    plot(tRRio,MBD_flow,'yo');
    plot(tRRio,DBD_flow,'go');
    hold off
    axis tight
    ylim([min(DBD_flow)*.9 max(SBD_flow)*1.1])
    title('Carotid Blood Flow')
    ylabel('Flow mL/min')
    %scrollplot
    
    % Plot of Blood Flow per Hundred Grams of Brain Tissue
    subplot(3,2,5)
    hold on
    plot(tint,CBF_Hgbt)
    plot(tRRio,SBD_flow_Hgbt,'ro');
    plot(tRRio,MBD_flow_Hgbt,'yo');
    plot(tRRio,DBD_flow_Hgbt,'go');
    hold off
    axis tight
    ylim([min(DBD_flow_Hgbt)*.9 max(SBD_flow_Hgbt)*1.1])
    title('Carotid Flow per Hundred Grams of Brain Tissue ')
    ylabel('Flow mL/g.min')
    %scrollplot
    hold off
    
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    uiwait(gcf)
    
end
close(gcf)

button = 0;
figure(16)
clf('reset')
while button == 0
    
    subplot(421)
    plot(tRRio,MBPo,'-b*');
    title('Mean Blood Pressure')
    ylabel('mmHg')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    
    subplot(422)
    plot(tRRio,MBVo,'-*');
    title('TCD Blood Velocity')
    ylabel('%')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    subplot(423)
    plot(tRRio,MBVel(:,2),'-b*');
    title('Mean Diameter')
    ylabel('Cm')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    
    subplot(424)
    plot(tRRio,MBVel(:,1),'-k*');
    title('Mean Blood Velocity')
    ylabel('Cm/s')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    subplot(425)
    plot(tRRio,MBD_flow,'-r*')
    title('Mean Blood flow')
    ylabel('Flow mL/min ')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    
    subplot(426)
    plot(tRRio,Distensibility,'-g*')
    title('Distensibility')
    ylabel('')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    subplot(427)
    plot(tRRio,Stiffness,'-b*')
    hold on
    title('Stiffness')
    ylabel('')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    subplot(428)
    plot(tRRio,CVR,'-m*')
    title('Resistance')
    ylabel('')
    axis tight
    YL = ylim;
    ylim([YL(1)*.95 YL(2)*1.05])
    %scrollplot
    
    
    pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    uiwait(gcf)
end
close(gcf)

%% Save data

% save filename
[pathstr, subject_file_name, ext] = fileparts(regexprep(infile, '_', ' '));
subject_file_dir = [destination_folder_name, '\', subject_file_name];
final_subject_file_name = [subject_file_name, 'AlignedClean.mat'];

% checking if file originally exist
check_exist = exist(subject_file_dir, 'file');

if check_exist ~= 0
    
    mkdir(subject_file_dir);
    save([subject_file_dir '\' final_subject_file_name])
    
    
elseif check_exist ==0
    save([subject_file_dir '\' final_subject_file_name])
end


end