% BrainFlowAnalysis04112016.m Revised to process autonomic, diameter, and velocity data

% Functions called:  cell2str, Minsma, screensize, fastsmooth, nanmedian,
%                    scrollplot, hline, Peaksma, siglocation, HRV, refline,
%                    zscore

% Testing:           09/22/2016  Dr. Serrador
% Written:          Dolu Obatusin, Bishoy Samy & Dr. Serrador
%-----------------------------------------------------------------------------------------
% This script will align datasets from LabChart and Velocity profiles from
% A DVD recording of Carotid arteries
%% Import LabChart Files
% Clear variables and functions from memory.
close all
clc
clearvars
clear all;

set(0, 'DefaulttextInterpreter', 'none') %change text interpreter so no longer creates subscripts in titles

initials = inputdlg('please enter you initials');%Operator enters initials for future purposes

TR= cell2mat(initials);


%% Adding list option to choose type of files user will input

file_str = {'CleanData','MATLAB Velocity','MAUI Diameter','Brachial Analyzer Velocity','Brachial Analyzer Diameter'};

[file_value,v] = listdlg('PromptString','Choose files to import:',...
    'SelectionMode','multiple',...
    'ListString',file_str);

% END of edits 1/10/17 by BS
%%

lc_data = find(file_value==1); 
lc_exist = isempty(lc_data);                                                         %checking to see if labchart data was initially chosen
if lc_exist == 1
    msgbox('Labchart data import not chosen, skipping Labchart data import');
else
    Notes_list ={};
    %    Prompt = {'Output File Name','ACA Channel','MCA Channel','VA Channel','MCA Channel','ECG Channel','CO2 Channel'};
    
    filterspec = '*clean.mat';
    Title = 'Pick file with CleanData';
    original_path = pwd;
    [infile2,pathname] = uigetfile(filterspec,Title); %Standard open file dialog box.
    
    if infile2 == 0
        %        errordlg('File not found','File Error');
        return;
    end
    infile1 = [pathname infile2];
    fid = fopen('tcd_ich.prf','r');
    if fid == -1
        %   filterspec = 'C:\IHPL\Research Studies\VRC\*.wdq';
        Def = {'.txt','2','3','4','5','6','7'};
    else
        load tcd_ich.prf -MAT
        if isempty(Def)
            Def = {'output.txt','2','3','4','5','6','7'};
        end
        if ~isempty(Def)~5;
            Def = {'output.txt','2','3','4','5','6','7'};
        end
        %    if isempty(filterspec)
        %       filterspec = 'C:\IHPL\Research Studies\VRC\*.wdq';
        %    end
        fclose(fid);
    end
    LineNo=1;
    load(infile1);
    %assignin ('base','d',d)
    
    %%5 Loaded data now assign all data to one variable called "data"
    
    %     if isfield (d,'data_block1') % lab chart 5 first case
    %         data =d.data_block1;
    %         titles=d.titles_block1;
    %         ticktimes=d.ticktimes_block1;
    %         t = ticktimes;
    %         comtext = d.comtext_block1;
    %         comtick=d.comtick_block1;
    %         index=d.index_block1;
    %     elseif isfield (d,'data') && isfield (d,'ticktimes')  % Lab chart 5 2nd case
    %         %fn = fieldnames(d);
    %         data =d.data;
    %         titles=d.titles;
    %         ticktimes=d.ticktimes;
    %         t = ticktimes;
    %         comtext=d.comtext;
    %         comtick=d.comtick;
    %         index=d.index;
    %     else                % Lab chart 7
    %         data = d.data;
    %         datastart=d.datastart;
    %         dataend =d.dataend;
    %         titles =d.titles ;
    %         % creating index
    %         index = 1 :length(titles);
    %         index=index';
    %         %________________
    %         total_time_sec= d.blocktimes(1);
    %         originaldata=data;
    %         size_of_data = size(datastart);
    %         numchanel = size_of_data (1);
    %         numblocks = size_of_data(2);
    %         data_channel =[];
    %         elo =numchanel;
    %         numchanel = length (dataend);
    %         for i = 1: numchanel
    %             for k = 1 : numblocks
    %                 data_block = data((datastart(i,k): dataend (i,k)));
    %                 data_channel = cat (2, data_channel,data_block);
    %             end
    %             arrat_raw_data (i) = mat2cell(data_channel);
    %             data_channel =[];
    %         end
    %
    %         channel1 = cell2mat (arrat_raw_data(1));
    %         channel2 = cell2mat (arrat_raw_data(1));
    %         channel3 = cell2mat (arrat_raw_data(3));
    %         channel4 = cell2mat (arrat_raw_data(4));
    %         channel5 = cell2mat (arrat_raw_data(5));
    %         channel6 = cell2mat (arrat_raw_data(6));
    %         %channel7 = cell2mat (arrat_raw_data(7));
    %         %channel8 = cell2mat (arrat_raw_data(8));
    %         t = 1: (total_time_sec-1)/(length(channel1)-1): total_time_sec;
    %         data =  ones(numchanel, length(channel1));
    %         for i =1 : numchanel
    %             data(i,:) = cell2mat(arrat_raw_data(i));
    %         end
    %
    %         t(end)=[];
    %     end
    %
    %     timeChart = t';
    %     timeECG =t';
    %     timeBP = t';
    %
    %     % for titles
    %     % if titles
    %
    %     titles_str = cellstr(titles);
    %     for i = 1: length(titles_str)
    %         if strcmp (titles_str(i,:),'BP')
    %             bpChart = data(i,:); %Measured in mmHg
    %         end
    %     end
    %
    %     for i = 1: length(titles_str)
    %         if strcmp (titles_str(i,:),'ECG')
    %             ecgChart = data(i,:);%Measured in beats per minute
    %         end
    %     end
    %
    %     for i = 1: length(titles_str)
    %         if strcmp (titles_str(i,:),'CO2')
    %             co2Chart = data(i,:);%Measured in mmHg
    %         end
    %     end
    %
    %     for i = 1: length(titles_str)
    %         if strcmp (titles_str(i,:),'MCA1')
    %             mcaoneChart = data(i,:);%Measured in cm/s
    %             MCA1 = mcaoneChart;
    %
    %
    %         elseif strcmp (titles_str(i,:),'Comment')
    %             cmtChart = data(4,:);
    %         end
    %
    %     end
    %
    %     for i = 1: length(titles_str)
    %         if strcmp (titles_str(i,:),'MCA2')
    %             mcatwoChart = data(i,:);%Measured in cm/s
    %             MCA2 = mcatwoChart;
    %         elseif strcmp (titles_str(i,:),'Marker')
    %             markerChart = data(5,:);
    %         end
    %
    %     end
    %
    %     % co2Chart = data(2,:);
    %     % ecgChart = data(3,:);
    %     % cmtChart = data(4,:);
    %     % markerChart = data(5,:);
    %     % end
    %     Chart_SR = 1 / diff(ticktimes(1:2));
    %     disp(['Chart signal sampling rate is ', num2str(Chart_SR), ' Hz'])
    %     ScriptName = mfilename;
    
end

%% Import DVD Velocity data (matlab velocity analysis)
%Read .mat files from Velocity Data
%FILTERSPEC ={'*.wav;*.avi;*.wma'};

mlv_data = find(file_value==2); mlv_exist = isempty(mlv_data);                                           %checking to see if matlab velocity data was initially chosen
if mlv_exist == 1
    msgbox('MATLAB velocity data import not chosen, skipping MATLAB velocity import');
else
    FILTERSPEC =('*.mat');
    TITLE = 'Pick file with Velocity Data';
    original_path = pwd;
    [Filename, pathname] = uigetfile(FILTERSPEC,TITLE);
    if Filename==0
        %     break;
    end
    if isequal(Filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, Filename)])
    end
    DVD = load([pathname Filename]);
    titlesDvd='Velocity';
    timeDVD = DVD.time_array ;% The same as DVD.T
    timeDVD = timeDVD';
    timeDVD_original = timeDVD;
    timeDVD = timeDVD(:,1);
    datDvd = DVD.Velocity_forRRI;
    % iRRi = DVD.iRRi;
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
    disp(['Velocity signal sampling rate is ', num2str(Vel_SR), ' Hz'])
    
end
% %% Import Excel(Diameter, Velocity) Data
% % Updated 07/12/2016 to prevent importing incorrect Velocity data from
% % Brachial Analyzer Excel form
%
% filterspec = '*.xlsx';
% Title = 'Pick the Diameter Data';
% [infile,pathname] = uigetfile(filterspec,Title);
%
% if infile == 0
%     return;
% end
% infile1 = [pathname infile];
% [STATUS,SHEETS] = xlsfinfo(infile1);
%
% choice = questdlg('Which program is used for the Diameter Tracking?', ...
%     'Software Identifier', ...
%     'Brachial Analyzer','MAUI','MAUI');
%
% switch choice
%     case 'Brachial Analyzer'
% % Diameter data is located in the second sheet
%
% Sheets = char(SHEETS(2));
% diam_exc = xlsread(infile1,Sheets);
%
% Sheets = char(SHEETS(2));
% [~, ~, diam_exc] = xlsread(infile1,Sheets);
% diam_exc2 = xlsread(infile1,Sheets);
%
% diam_exc(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),diam_exc)) = {''};
%
% % Allocate imported array to column variable names
% col = size (diam_exc);
% try
%     for i = 1: col(2)
%         cellVector = diam_exc(:,i);
%         % combinedStr = strcat(s1,s2)
%         field = ['VarName' int2str(i)];  %# Create the name string
%         value = cellVector;%# Create the Value cell array
%         s.(field) = struct(field,{value});
%         for j = 1: col(2)
%             header = diam_exc(1,j);
%             header = [header{1}];
%             header = strrep(header, ' ', '_');
%             header = strrep(header, '/', '_');
%             header = strrep(header, '(', '_');
%             header = strrep(header, ')', '_');
%             header = strrep(header, '-', '_');
%             header = strrep(header, '.', '_');
%             header = strrep(header, '%', '_');
%             header = strrep(header, ':', '_');
%
%             %         units = diam_exc(2,j);
%             %         units = [units{1}];
%             %         units = strrep(units, ' ', '_');
%             %         units = strrep(units, '/', '_');
%             %         units = strrep(units, '(', '_');
%             %         units = strrep(units, ')', '_');
%             %         units = strrep(units, '-', '_');
%             %         units = strrep(units, '.', '_');
%             %         units = strrep(units, '%', '_');
%             %         units = strrep(units, ':', '_');
%
%             cellVector = diam_exc(4:end,j);
%             data.(header) = struct('value',{cellVector});
%             %         data.(header).units= struct('units',{units});
%         end
%
%     end
% catch err
%     if (strcmp(err.identifier,'MATLAB:AddField:InvalidFieldName'))
%
%         uiwait( errordlg('This file does not have Velocity data. Please try again.'));
%         return;
%     end
% end
%
% try
% %	Identify the column of Diameter and Velocity Data
% pat1 = 'B\w*M';
% pat2 = 'F\w*G';
% pat3 = 'M\w*C';
%
% datastrct = struct2nv(data);
% DiamIndx = regexpi(datastrct, pat1);
% VelIndx = regexpi(datastrct, pat2);
% timeIndx = regexpi(datastrct, pat3);
%
% % identify the location of the Signal
% DiamIndx = find(not(cellfun('isempty', DiamIndx)));
% VelIndx = find(not(cellfun('isempty', VelIndx)));
% timeIndx = find(not(cellfun('isempty', timeIndx)));
%
% %Data from Brachial Analyzer
% time_diam = diam_exc2(:,timeIndx)/1000; %msec--> divide by 1000 to convert from msec to sec
% diam = diam_exc2(:,DiamIndx); %mm
%
%
%
% % Find the right velocity column
% Velocity_BA = diam_exc2(:,VelIndx); %m/s
% Velocity_BA = Velocity_BA * 100; % Convert to cm/s
%
% % Remove all blanks found in either diameter or velocity
% %Find empty time cells, if any
% T_nan2 =isnan(time_diam);
% [  I2 ] = find(T_nan2>= 1);
% time_diam(I2)=[];
% diam(I2)=[];
% Velocity_BA(I2)=[];
% %Find empty diameter cells, if any
% D_nan =isnan(diam);
% [ D ] = find(D_nan>= 1);
% time_diam(D)= [];
% diam(D)= [];
% Velocity_BA(D)=[];
% %Find empty Velocity cells, if any
% V_nan =isnan(Velocity_BA);
% [ V ] = find(V_nan>= 1);
% Velocity_BA(V)=[];
% time_diam(V)= [];
% diam(V)= [];
%
% % assign variables that have been produced by brachial analyzer software
% % analysis
% timeBA=time_diam; %time based on the brachial analzyer
% diamBA=diam/10; %diameter (in cm)  derived from brachial analyzer
%
% % Check if length of velocity and length of diamter from BA match....
% %Error message will end the program.
% if length(diamBA)  ~= length(Velocity_BA)
%     uiwait( errordlg('Diameter and Velocity Size mismatch, Please call an Engineer.'));
%     return;
% end
%
% diam_SR = 1 / diff(time_diam(1:2));
% disp(['Brachial Analyzer signal sampling rate is ', num2str(diam_SR), ' Hz'])
% catch err
%
%     if (strcmp(err.identifier,'MATLAB:catenate:dimensionMismatch'))
%
%         try
%             % Assign Signals and position coded values
%             %Common = 1;
%             %Internal= 2;
%             %Seated =1;
%             %Supine=2;
%
%             import = importdata(infile);
%             T = readtable(infile);
%             Procedure = T{3,2};
%
%
%             x = cell2mat(Procedure);
%             y = length(x);
%             Position_length = y - 5;%Six because both supine and seated are lengths of six
%             Position_str = x(Position_length:end);
%             % Procedure = cell2mat(Procedure);
%
%             % Match regular expression, ignoring case.
%             Procedure;
%             pat1 = 'C\w*n';
%             pat2 = 'I\w*l';
%             pat3 = 'S\w*d';
%             pat4 = 'S\w*e';
%
%             if cell2mat(regexpi(Procedure, pat1)) > 0 ;
%                 Artery = 1;
%             elseif cell2mat(regexpi(Procedure, pat2))> 0;
%                 Artery = 2;
%             end
%
%             if cell2mat(regexpi(Procedure, pat3)) > 0 ;
%                 Position = 1;
%             elseif cell2mat(regexpi(Procedure, pat4))> 0;
%                 Position = 2;
%             end
%
%         catch exception2
%             rethrow(exception2);
%
%         end
%     end
%
% end
%% Import Excel(Diameter, Velocity) Data

bav_data = find(file_value==4); bav_exist = isempty(bav_data);                 %checking to see if BA velocity data was initially chosen
ba_dia_data = find(file_value==5); ba_dia_exist = isempty(ba_dia_data);        %checking to see if BA diameter data was initially chosen

if bav_exist == 1 && ba_dia_exist == 1
    msgbox('BA velocity and diameter data not chosen to import, skipping import!');
elseif bav_exist == 0 && ba_dia_exist == 0
    
    filterspec = {'*.xls;*.csv;*.xlsx;'};
    Title = 'Pick the BA Diameter/Velocity Data';
    [infile,pathname] = uigetfile(filterspec,Title);
    
    if infile == 0
        return;
    end
    infile1 = [pathname infile];
    [STATUS,SHEETS] = xlsfinfo(infile1);
    
    % Diameter data is located in the second sheet
    
    Sheets = char(SHEETS(1));
    [T1,T2,T] = xlsread(infile1,Sheets,'B:B');
    
    Sheets = char(SHEETS(2));
    diam_exc = xlsread(infile1,Sheets);
    
    %Data from Brachial Analyzer
    time_diam = diam_exc(:,8)/1000; %msec--> divide by 1000 to convert from msec to sec
    diam = diam_exc(:,2); %mm
    Velocity_BA = diam_exc(:,end); %m/s
    Velocity_BA = Velocity_BA * 100; % Convert to cm/s
    
    % Remove all blanks found in either diameter or velocity
    %Find empty time cells, if any
    T_nan2 =isnan(time_diam);
    [  I2 ] = find(T_nan2>= 1);
    time_diam(I2)=[];
    diam(I2)=[];
    Velocity_BA(I2)=[];
    %Find empty diameter cells, if any
    D_nan =isnan(diam);
    [ D ] = find(D_nan>= 1);
    time_diam(D)= [];
    diam(D)= [];
    Velocity_BA(D)=[];
    %Find empty Velocity cells, if any
    V_nan =isnan(Velocity_BA);
    [ V ] = find(V_nan>= 1);
    Velocity_BA(V)=[];
    time_diam(V)= [];
    diam(V)= [];
    
    % assign variables that have been produced by brachial analyzer software
    % analysis
    timeBA=time_diam; %time based on the brachial analzyer
    diamBA=diam/10; %diameter (in cm)  derived from brachial analyzer
    
    % Check if length of velocity and length of diamter from BA match....
    %Error message will end the program.
    if length(diamBA)  ~= length(Velocity_BA)
        uiwait( errordlg('Diameter and Velocity Size mismatch, Please call an Engineer.'));
        return;
    end
    
    diam_SR = 1 / diff(time_diam(1:2));
    disp(['Brachial Analyzer signal sampling rate is ', num2str(diam_SR), ' Hz'])
    
elseif bav_exist == 0 && ba_dia_exist == 1                               %BA vel exists but no diameter
    
    filterspec = {'*.xls;*.csv;*.xlsx;'};
    Title = 'Pick the BA Velocity Data';
    [infile,pathname] = uigetfile(filterspec,Title);
    
    if infile == 0
        return;
    end
    infile1 = [pathname infile];
    [STATUS,SHEETS] = xlsfinfo(infile1);
    
    % Diameter data is located in the second sheet
    
    Sheets = char(SHEETS(1));
    [T1,T2,T] = xlsread(infile1,Sheets,'B:B');
    
    Sheets = char(SHEETS(2));
    diam_exc = xlsread(infile1,Sheets);
    
    %Data from Brachial Analyzer
    time_diam = diam_exc(:,8)/1000; %msec--> divide by 1000 to convert from msec to sec
    %diam = diam_exc(:,2); %mm
    Velocity_BA = diam_exc(:,end); %m/s
    Velocity_BA = Velocity_BA * 100; % Convert to cm/s
    
    % Remove all blanks found in either diameter or velocity
    %Find empty time cells, if any
    T_nan2 =isnan(time_diam);
    [  I2 ] = find(T_nan2>= 1);
    time_diam(I2)=[];
    %diam(I2)=[];
    Velocity_BA(I2)=[];
    %Find empty diameter cells, if any
    %D_nan =isnan(diam);
    %[ D ] = find(D_nan>= 1);
    %time_diam(D)= [];
    %diam(D)= [];
    %Velocity_BA(D)=[];
    %Find empty Velocity cells, if any
    V_nan =isnan(Velocity_BA);
    [ V ] = find(V_nan>= 1);
    Velocity_BA(V)=[];
    time_diam(V)= [];
    %diam(V)= [];
    
    % assign variables that have been produced by brachial analyzer software
    % analysis
    timeBA=time_diam;           %time based on the brachial analzyer
    timeBAvel=time_diam;
    %diamBA=diam/10; %diameter (in cm)  derived from brachial analyzer
    
    % Check if length of velocity and length of diamter from BA match....
    %Error message will end the program.
    if length(timeBA)  ~= length(Velocity_BA)
        uiwait( errordlg('Velocity and Time Size mismatch, Please call an Engineer.'));
        return;
    end
    
    diam_SR = 1 / diff(time_diam(1:2));
    disp(['Brachial Analyzer signal sampling rate is ', num2str(diam_SR), ' Hz'])
    
    
elseif bav_exist == 1 && ba_dia_exist == 0                               %BA diam exists but no velocity
    
    filterspec = {'*.xls;*.csv;*.xlsx;'};
    Title = 'Pick the BA Diameter Data';
    [infile,pathname] = uigetfile(filterspec,Title);
    
    if infile == 0
        return;
    end
    infile1 = [pathname infile];
    [STATUS,SHEETS] = xlsfinfo(infile1);
    
    % Diameter data is located in the second sheet
    
    Sheets = char(SHEETS(1));
    [T1,T2,T] = xlsread(infile1,Sheets,'B:B');
    
    Sheets = char(SHEETS(2));
    diam_exc = xlsread(infile1,Sheets);
    
    %Data from Brachial Analyzer
    time_diam = diam_exc(:,8)/1000; %msec--> divide by 1000 to convert from msec to sec
    diam = diam_exc(:,2); %mm
    %Velocity_BA = diam_exc(:,end); %m/s
    %Velocity_BA = Velocity_BA * 100; % Convert to cm/s
    
    % Remove all blanks found in either diameter or velocity
    %Find empty time cells, if any
    T_nan2 =isnan(time_diam);
    [  I2 ] = find(T_nan2>= 1);
    time_diam(I2)=[];
    diam(I2)=[];
    %Velocity_BA(I2)=[];
    %Find empty diameter cells, if any
    D_nan =isnan(diam);
    [ D ] = find(D_nan>= 1);
    time_diam(D)= [];
    diam(D)= [];
    %Velocity_BA(D)=[];
    %Find empty Velocity cells, if any
    %V_nan =isnan(Velocity_BA);
    %[ V ] = find(V_nan>= 1);
    %Velocity_BA(V)=[];
    %time_diam(V)= [];
    %diam(V)= [];
    
    % assign variables that have been produced by brachial analyzer software
    % analysis
    timeBA=time_diam; %time based on the brachial analzyer
    diamBA=diam/10; %diameter (in cm)  derived from brachial analyzer
    
    % Check if length of velocity and length of diamter from BA match....
    %Error message will end the program.
    if length(timeBA)  ~= length(diamBA)
        uiwait( errordlg('Diameter and Time Size mismatch, Please call an Engineer.'));
        return;
    end
    
    diam_SR = 1 / diff(time_diam(1:2));
    disp(['Brachial Analyzer signal sampling rate is ', num2str(diam_SR), ' Hz'])
    
end

% MAUI diameter data import
mauid_data = find(file_value==3); mauid_exist = isempty(mauid_data);        %checking to see if MAUI diameter data was initially chosen

if mauid_exist == 1
    msgbox('MAUI diam data not chosen to import, skipping import!');
else
    filterspec = {'*.xls;*.csv;*.xlsx;'};
    Title = 'Pick the MAUI Data';
    [infile,pathname] = uigetfile(filterspec,Title);
    
    if infile == 0
        return;
    end
    infile1 = [pathname infile];
    [STATUS,SHEETS] = xlsfinfo(infile1);
    
    % Diameter data is located in the second sheet
    Sheets = char(SHEETS(1));
    diam_exc = xlsread(infile1,Sheets);
    
    diam = diam_exc(:,6)'; %mm
    time_diam = diam_exc(:,5)'; %sec
    
    % assign variables that have been produced by brachial analyzer software
    % analysis
    timeBA=time_diam; %time based on the brachial analzyer
    
    %temporary correction since there seems to be a 150 ms delay for maui
    %analysis to occur compared to matlab velcoity
    timeBA = timeBA -0.15;
    
    
    diamBA=diam/10; %diameter (in cm)  derived from brachial analyzer
end
%% User prompt to define Artery and position during protocol

% Assign Signals and position coded values
%Common = 1;
%Internal= 2;
%Seated =1;
%Supine=2;
str = {'Common Seated','Internal Seated', 'Common Supine', 'Internal Supine', 'Common Standing', 'Internal Standing'};
choice = listdlg('PromptString','Select Artery and Position:',...
    'SelectionMode','single',...
    'ListString',str);
Procedure = cell2mat(str(choice));
switch cell2mat(str(choice))
    case 'Common Seated'
        Artery = 1;
        Position = 1;
    case 'Internal Seated'
        Artery = 2;
        Position = 1;
    case 'Common Supine'
        Artery = 1;
        Position = 2;
    case 'Internal Supine'
        Artery = 2;
        Position = 2;
    case 'Common Standing'
        Artery = 1;
        Position = 3;
    case 'Internal Standing'
        Artery = 2;
        Position = 3;
end

%% Change Directory
name = infile1;
folder_name = [infile(1:end-5), '_', str2mat(cell2mat(initials))];

for i = 1:length(name)
    if name(end-i) ~= '\'
    elseif name(end-i) == '\'
        cont =i;
        break
    end
end
name(end-cont:end)=[];
cd(name);

folder = mkdir (name,folder_name);
new_dir=[name '\' folder_name];
cd(new_dir)

%% Including Alignment between BP and ECG
% Finometer Pro has one beat delay in BP
% Aligning BP to ECG
% By: Bishoy Samy
% Date: 10/30/17
% Added from cleandata

SR = round(1./mean(diff(t)));
if SR == 100
    BVcm_int = BVcm(ikp,:);    %added by BS, never indicated BVcm of interest only BV
    %      ECG=ECG(ikp);
elseif SR ==1000
end

t_test=t';

if length(ECG) ~= length(t_test)
    ECG = ECG(ikp);
end

ECG_new=ECG;
ECG_new=ECG_new./max(ECG_new);
ECG_new=((max(BP)-min(BP)).*ECG_new)+min(BP);

dumy_BP=BP'; ECG_new=ECG_new';
BP_delay_shift=0;
BP_shifter=0;

ECG_orig1 = ECG_new;

%  BPECG_align = questdlg('Does BP need to be aligned to ECG?', ...
%                          'BP and ECG Alignment','Yes', 'No', 'No');
%    switch BPECG_align,
%      case 'Yes',
%       BPECG_align=1;
%      case 'No',
%       BPECG_align=0;
%    end % switch


%%%%% Generate and plot data
% dx is the width of the axis 'window'
button_2=0;
dx=3;
ha = [0 t(end)];
basefilter1 = 0;

TIMEP_length = length(timeP);
TIMETEST_length = length(t_test);



while button_2 == 0
    
    figure(10); clf
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
    legend('BP','ECG'); xlabel('Time'); title('Check for ECG to BP Alignment'); %scrollplot
    
    set(gcf,'doublebuffer','on');
    % This avoids flickering when updating the axis
    xlim([ha]);%set(a,'xlim',ha);
    %set(a,'ylim',[min(ECG_new*mult_factor) max(ECG_new*mult_factor)]);
    %%%%% Generate constants for use in uicontrol initialization
    pos=get(a,'position');
    Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
    % This will create a slider which is just underneath the axis
    % but still leaves room for the axis labels above the slider
    xmax=max(t_test);
    %ha = get(gca,'XLim');
    S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %     pb3l = uicontrol(gcf,'Style','pushbutton','String','shift','units','normalized','Position',[.925 .33 .05 .03],'Callback','dumy_MCA=circshift(dumy_MCA,round(SR)); shifter = shifter+1; uiresume;');
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
    
    %XL = xlim(ha);
    
end

close(gcf);

BP=dumy_BP';


clear button_2 a dx ha h pb1au pb2au pb3au pb4au pbexit pbreset pbwin pbwin bb aa ECG_orig1 ha_c xmax basefilter1 t_test ECG_new dumy_BP
%% Plot All signals
% plotting in Percent
txtfontsize= 11;
t_name1 = ' ' ;
t_name2 = ' ' ;

figure(1)
clf

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
k =num_artery+3;
subplot(k,1,1)
plot(t,bp_f, 'b', lw, lw_thin), hold on;
plot(tRRi,MBP,'k', lw, lw_thick),
plot(tRRi,SBP,'r', lw, lw_thick),
plot(tRRi,DBP,'r', lw, lw_thick), hold off;
if isempty(ticktimes(comtick))==0
    for txt=1:length(markers)
        tetx_handle = text (t(markers_Index(txt)) ,bp_f(markers_Index(txt)),markers{txt});
        set(tetx_handle,'Colo','r','FontSize',txtfontsize,'Rotation',90)
    end
end
axis tight
yl = ylim;
yl(1) = yl(1)-yl(2)*.05;
yl(2) = yl(2)+yl(2)*.05;
ylim(yl);
%    xlabel('Time (s)', fs, fs_label)
ylabel({'Arterial';'Blood Pressure';'(mmHg)'}, fs, fs_label)
title(plot_title, fs, fs_title)


for figcont = 1:num_artery
    subplot(k,1,figcont+1)
    plot(t,bv_f(:,figcont), 'g', lw, lw_thin), hold on;
    plot(tRRi,MBV(:,figcont),'k', lw, lw_thick),
    plot(tRRi,SBV(:,figcont),'r', lw, lw_thick),
    plot(tRRi,DBV(:,figcont),'r', lw, lw_thick),
    line([0 t(end)], [100 100], 'Color','m'),hold off;
    axis tight
    if isempty(ticktimes(comtick))==0
        for txt=1:length(markers)
            tetx_handle = text (t(markers_Index(txt)) ,bv_f(markers_Index(txt),figcont),markers{txt});
            set(tetx_handle,'Colo','r','FontSize',txtfontsize,'Rotation',90)
        end
    end
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    xlabel('Time (s)', fs, fs_label);
    ylabel({cell2mat(labels (figcont));'(%)'}, fs, fs_label);
    
    
    
end


% thi spart to plot Cratiod not normalized
if exist('Carotid')
    
    
    subplot(k,1,figcont+1)
    plot(ti,Carotid_clean, 'g', lw, lw_thin), hold on;
    plot(tRRi,M_Carotid_Clean,'k', lw, lw_thick),
    plot(tRRi,S_Carotid_Clean,'r', lw, lw_thick),
    plot(tRRi,D_Carotid_Clean,'r', lw, lw_thick),,hold off;
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
% ends here



subplot(k,1,num_artery+2)
plot(tRRi,1000./RRi*60,'r', lw, lw_thick);
axis tight
yl = ylim;
yl(1) = yl(1)-yl(2)*.05;
yl(2) = yl(2)+yl(2)*.05;
ylim(yl);
ylabel('HR (bpm)', fs, fs_label);
%    xlabel('Time (s)', fs, fs_label);
subplot(k,1,num_artery+3)
plot(tco2,CO2, 'g', lw, lw_thin), hold on;
plot(tetCO2,etCO2,'r', lw, lw_thick), hold off;
axis tight
if isempty(ticktimes(comtick))==0
    for txt=1:length(markers)
        tetx_handle = text (tco2(markers_Index(txt)) ,CO2(markers_Index(txt)),markers{txt});
        set(tetx_handle,'Colo','r','FontSize',txtfontsize,'Rotation',90)
    end
end
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

% button = myquestdlg('Is the data correctly peak detected?','Continue Operation','Yes','No','Help','No');
% if strcmp(button,'Yes')
%     disp('Matlab are currently saving your files !!')
% elseif strcmp(button,'No')
%     disp('Ending Program')
%     %break
% elseif strcmp(button,'Help')
%     disp('Sorry, no help available')
%     %break
% end

%Trim Data to just include aligned ECG and BP data
% st1 = max([t(1) tECG(1)]);
% st2 = min([t(end) tECG(end)]);
% good = find(tECG>=st1 & tECG<=st2);
% ECG = ECG(good);
% clear good
% good = find(t>=st1 & t<=st2);
% t = t(good);
% BP = BP(good);
%
% markerChart = markerChart(good);
% clear good
%% Select channel with marker if there is one

numch = min(size(data));
choosing_marker = 0;
Marker_Channel_Use=zeros(numch,1);

while choosing_marker<1
    figure(33)
    clf
    
    for i=1:numch
        subplot(numch,1,i)
        plot(data(i,:))
        title(titles(i,:))
        axis tight
        possub = get(gca,'position');
        BU(i) = uicontrol(gcf,'Style','checkbox','String',titles(i,:),'BackgroundColor','c','FontSize',10,'units','normalized','Position',[.925 (possub(2)+possub(4)/3) .1 .05],'Callback', ' uiresume;');
        set(BU(i),'value',Marker_Channel_Use(i));
    end
    
    b4u = uicontrol(gcf,'Style','pushbutton','String',' Finish ','BackgroundColor','r','FontSize',15,'units','normalized','Position',[.025 .05 .05 .05],'Callback', 'choosing_marker =2 ; uiresume;');
    uiwait(gcf)
    
    
    for chkboxnum=1:(numch)
        Marker_Channel_Use(chkboxnum) =get(BU(chkboxnum),'value');
    end
end

MarkerChannel = find(Marker_Channel_Use==1);

if isempty(MarkerChannel)==1
    markerChart=NaN;
else
    markerChart=abs(data(MarkerChannel,:));
end


%% Do a quick alignment based on the marker channels
%Velocitydlg  =questdlg('What signal would you like to choose for Velocity Peak Detection?','Please choose the appropriate Velocity Signal' ,'MATLAB','Brachial_Analyzer','');
% Assign variables for marker alignments

%this will be skipped for files with no matlab velocity analysis
if mlv_exist == 1 %if matlab velocity matrix is empty, program will assume BA vel is present
    timeMat = timeBAvel; % Matlab Time marker
    VelMat = Velocity_BA; % Matlab Velocity marker Signal
    Velocitydlg='Brachial_Analyzer';
else
    Velocitydlg='MATLAB';
end

timemarker = timeMat; % Matlab Time marker
SignalMarker = VelMat; % Matlab Velocity marker Signal

%% align data based on marker if possible

if mlv_exist == 1     %checks to see if matlab velocity is empty, if it is, there is no marker
    markerChart=nan;
end

buttonqdlg  = questdlg('Do you have a marker?','Marker Option','Yes','No','Help', 'No');

if strcmp(buttonqdlg,'Yes')
    disp('Proceeding to Alignment')
    
    
    
    
    % Verify the presence of a marker signal
    if isnan(markerChart(1))
        buttonqdlg = 'No';
        
    else
        buttonqdlg = 'Yes';
    end
    
    if strcmp(buttonqdlg,'Yes')
        
        button = 0;
        
        XL(1) = min([min(t_1000HZ) min(timemarker)]);
        XL(2) = max([max(t_1000HZ) max(timemarker)]);
        
        figure('Name', 'Quickalignmentbasedonthemarkerchannels','NumberTitle','on','NextPlot', 'add')
        
        while button==0
            
            
            clf
            %     scrsz = get(0,'ScreenSize');
            %     set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
            plot(t_1000HZ,markerChart./max(markerChart),'b-')
            hold on
            %     plot(timeMat,Velocity_PeakDetect./max(Velocity_PeakDetect),'r-')
            plot(timemarker,SignalMarker./max(SignalMarker),'r-')
            
            hold off
            xlim(XL)
            ha = get(gcf,'CurrentAxes');
            %scrollplot;
            
            pb2u = uicontrol(gcf,'Style','pushbutton','String','Shift+','units','normalized','Position',[.925 .9 .05 .03],'Callback','timemarker = timemarker+0.01; uiresume;');
            pb2l = uicontrol(gcf,'Style','pushbutton','String','Shift-','units','normalized','Position',[.925 .8 .05 .03],'Callback','timemarker = timemarker-0.01; uiresume;');
            pb3u = uicontrol(gcf,'Style','pushbutton','String','Shift++','units','normalized','Position',[.925 .7 .05 .03],'Callback','timemarker = timemarker+1; uiresume;');
            pb3l = uicontrol(gcf,'Style','pushbutton','String','Shift--','units','normalized','Position',[.925 .6 .05 .03],'Callback','timemarker = timemarker-1; uiresume;');
            pb4u = uicontrol(gcf,'Style','pushbutton','String','Shift+++','units','normalized','Position',[.925 .5 .05 .03],'Callback','timemarker = timemarker+10; uiresume;');
            pb4l = uicontrol(gcf,'Style','pushbutton','String','Shift---','units','normalized','Position',[.925 .4 .05 .03],'Callback','timemarker = timemarker-10; uiresume;');
            
            pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
            
            uiwait(gcf)
            
            XL = xlim(ha);
            
        end
    end
    %    savefig('Quick alignment based on the marker channels')
    close(gcf)
elseif strcmp(buttonqdlg,'No')
    disp('Proceeding to Peak Detection')
    timeMarker = timeMat; % Matlab Time marker
    markerChart = nan(size(markerChart));
    if exist('marker')==0
        marker=nan;
    end
    marker = nan(size(marker));
    
elseif strcmp(buttonqdlg,'Help')
    disp('Quick, call an engineer! x2184')
end



if strcmp(Velocitydlg,'MATLAB')
    %Section 5:  Align diameter with Velocity signal if using velocity signal from matlab
    
    Vel = VelMat;
    diam = diamBA;
    
    if exist ('markerChart', 'var')
        marker = markerChart;
    else
        n = size (diam);
        marker = NaN(n);
        markerChart = NaN(n);
        
    end
    % Align diam with Vel signal
    button = 0;
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
    %     Vel = Vel(good);
    %     Velocity_PeakDetect = Vel;
    
    %     good = find(timemarker>=st1 & timemarker<=st2);
    timeMarker = timemarker(good);
    clear good
    
elseif strcmp(Velocitydlg,'Brachial_Analyzer')
    Velocity_PeakDetect = Velocity_BA';
    tdiam = timeBAvel;
    
    %Trim Data to just include aligned diam and Vel data
    st1 = max([t_1000HZ(1) timemarker(1)]);
    st2 = min([t_1000HZ(end) timemarker(end)]);
    good = find(tdiam>=st1 & tdiam<=st2);
    diam = diamBA;
    Velocity_PeakDetect  = Velocity_BA(good);
    diam = diam(good);
    timeMarker = tdiam(good);
    clear good
    
end


% %% Choose Signal for Peak Detection
% %     if strcmp(buttonqdlg,'No')
% Peak_detection_dlg  =questdlg('What signal would you like to choose for Peak Detection?','please choose peak detection signal' ,'BP','ECG','');
%
% if strcmp(Peak_detection_dlg,'BP')
%     Peak_detection_signal =BP;
% elseif strcmp(Peak_detection_dlg,'ECG')
%     if exist('ECG')
%         Peak_detection_signal= ECG;
%     else
%         Peak_detection_signal=BP;
%     end
%
% elseif strcmp(Peak_detection_dlg,'MCA2')
%     if exist('MCA2')
%         Peak_detection_signal= MCA2;
%     else
%         Peak_detection_signal=BP;
%     end
% end
%
% % Peak_detection_signal = interp1(t_1000HZ,Peak_detection_signal_1000HZ,t,'linear');
% tco2 = t;
%
% % Perform peak detection
%
% Signal = Peak_detection_signal;
% SR = round(1./mean(diff(t)));
%
% points = ceil(0.2*SR);
% thresh = max(Signal)*0.35;
% maxthresh = max(Signal)*1.1;
% minthresh = min(0,min(Signal));
% Value = 'value';
% Remove_Group = 0;
% button = 0;
% basefilter = 0;
% LPfilter = 0;
%
% iRRi = Peaksmsa(Signal,points,thresh);
% XL = [0 length(Signal)];
% XL2 = [0 max(t)];
% redetect = 0;
%
% HR = 1./diff(t(iRRi))*60;
% RR = diff(t(iRRi))*1000;
%
% Signal_Original = Signal;
% k=50;
% cutoff_freq = SR/4;
%
% minC    = nanmedian(Signal) - 2*nanstd(Signal);
% maxC    = nanmedian(Signal) + 2*nanstd(Signal);
% Signal_Norm = (Signal - minC) /(maxC - minC);
%
% minC    = nanmedian(BP) - 2*nanstd(BP);
% maxC    = nanmedian(BP) + 2*nanstd(BP);
% BPdetect = (BP - minC) /(maxC - minC);
%
% Signal_Norm = Signal;
% % thresh = max(timeMat)*0.35;
% %% Peak Detection of chosen signal
% % SR = round(1./mean(diff(t)));
% m = 1; %value preset for use in basefilter (Highpass) function.
% n = 1; %value preset for use in LPF function.
% figure('Name', 'PeakDetectionofchosensignal','NumberTitle','on','NextPlot', 'add')
% button = 0;
% while button==0
%
%     if redetect == 1
%         iRRi = Peaksmsa(Signal,points,thresh);
%         bad = find(diff(t(iRRi))<0.15);
%         iRRi(bad+1) = [];
%         HR = 1./diff(t(iRRi))*60;
%         redetect = 0;
%     end
%     clf
%
%     subplot(3,1,1:2)
%     plot(Signal);
%     hold on
%     hline = refline(0,thresh);
%     set(hline,'Color','g');
%     plot(iRRi,Signal(iRRi),'ro')
%     hold off
%     scrollplot;
%
%     xlim(XL)
%     ha = get(gcf,'CurrentAxes');
%
%     subplot(3,1,3)
%     plot(t(iRRi(2:end)),HR)
%     ylabel('LabChart RRI')
%     scrollplot;
%
%     xlim(XL2)
%     ha2 = get(gcf,'CurrentAxes');
%     ss = get(gcf,'Position');
%
%     pb1u = uicontrol(gcf,'Style','slider','Max',ceil(max(Signal)*1.2),'Min',floor(min(Signal)-max(Signal)*0.2),'Value',round(thresh),'SliderStep',[0.001 0.1],'units','normalized','Position',[.925 .5 .025 .4],'Callback','thresh = get(pb1u,Value); uiresume;');
%     pb1ut1 = uicontrol(gcf,'style','text','string',ceil(YL1(2)),'units','normalized','HorizontalAlignment','left','Position',[.96 .86 .025 .025]);
%     pb1utime = uicontrol(gcf,'style','text','string',floor(YL1(1)),'units','normalized','HorizontalAlignment','left','Position',[.96 .5 .025 .025]);
%     pb1ut3 = uicontrol(gcf,'style','text','string',round(thresh),'units','normalized','HorizontalAlignment','left','Position',[.96 .68 .025 .025]);
%
%     pb2u = uicontrol(gcf,'Style','pushbutton','String','Delay+','units','normalized','Position',[.925 .45 .05 .03],'Callback','points = points + points*.1; redetect = 1; uiresume;');
%     pb2l = uicontrol(gcf,'Style','pushbutton','String','Delay-','units','normalized','Position',[.925 .41 .05 .03],'Callback','points = points - points*.1; redetect = 1; uiresume;');
%     pb3u = uicontrol(gcf,'Style','pushbutton','String','Delay++','units','normalized','Position',[.925 .37 .05 .03],'Callback','points = points + points*.5; redetect = 1; uiresume;');
%     pb3l = uicontrol(gcf,'Style','pushbutton','String','Delay--','units','normalized','Position',[.925 .33 .05 .03],'Callback','points = points - points*.5; redetect = 1; uiresume;');
%
%     pb5u = uicontrol(gcf,'Style','pushbutton','String','Thresh+','units','normalized','Position',[.925 .29 .05 .03],'Callback','thresh = thresh + thresh*.05; redetect = 1; uiresume;');
%     pb5l = uicontrol(gcf,'Style','pushbutton','String','Thresh-','units','normalized','Position',[.925 .25 .05 .03],'Callback','thresh = thresh - thresh*.05; redetect = 1; uiresume;');
%
%     pb4u = uicontrol(gcf,'Style','pushbutton','String','LPF','units','normalized','Position',[.925 .2 .05 .03],'Callback','n = n + 1;LPfilter = 1; redetect = 1; uiresume;');
%     pb4l = uicontrol(gcf,'Style','pushbutton','String','BaseF','units','normalized','Position',[.925 .16 .05 .03],'Callback','m = m + 1; basefilter = 1; redetect = 1; uiresume;');
%     pb6l = uicontrol(gcf,'Style','pushbutton','String','ReDetect','units','normalized','Position',[.925 .12 .05 .03],'Callback','redetect = 1; uiresume;');
%
%     pb10 = uicontrol(gcf,'Style','pushbutton','String','inverse','units','normalized','Position',[.2 .025 .05 .03],'Callback','Signal = -Signal;Signal_Original = -1*Signal_Original ; uiresume;');
%     %         pb7 = uicontrol(gcf,'Style','pushbutton','String','Left Artery','units','normalized','Position',[.4 .025 .05 .03],'Callback','ECG = BV(:,1); thresh = max(ECG)*0.35; uiresume;');
%     % %       pb8 = uicontrol(gcf,'Style','pushbutton','String','Right Artery','units','normalized','Position',[.6 .025 .05 .03],'Callback','ECG = BV(:,2); thresh = max(ECG)*0.35; uiresume;');
%     %         pb9 = uicontrol(gcf,'Style','pushbutton','String','ECG','units','normalized','Position',[.8 .025 .05 .03],'Callback','ECG = BP;  thresh = max(BP)*0.35; uiresume;');
%     pb8 = uicontrol(gcf,'Style','pushbutton','String','ECG','units','normalized','Position',[.6 .025 .05 .03],'Callback','Signal = ECG; thresh = max(ECG)*0.35; uiresume;');
%     pb9 = uicontrol(gcf,'Style','pushbutton','String','BP','units','normalized','Position',[.8 .025 .05 .03],'Callback','Signal = BP;  thresh = max(BP)*0.35; uiresume;');
%
%     pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
%
%     uiwait(gcf)
%
%     ECGf = zeros(length(Signal),1);
%
%     if basefilter == 1
%         %            set(gcf,'Pointer','watch')
%         %         baseECG = medfilt1(ECG,round(SR));
%         %         ECGf=ECG - baseECG + mean(baseECG);
%         %       use a 5 sec moving median to remove baseline changes
%
%         ECGold = Signal;
%         %             for i = 1:length(ECG)
%         %                 if i+round(SR)*3 <= length(ECG)
%         %                     ECGf(i) = median(ECG(i:i+round(SR)*3));
%         %                 else
%         %                     ECGf(i) = median(ECG(length(ECG)-(round(SR)*3):length(ECG)));
%         %                 end
%         %             end
%         %             szECGf = size(ECGf);
%         %             if szECGf(1)>1
%         %                 ECGf = ECGf.';
%         %             end
%         %             ECG = ECG - ECGf;
%         basefilter = 0;
%         if m > 5
%             m = 5
%         end;
%         [b,a] = butter(m,5/SR,'high');
%         ECGf = filtfilt(b,a,Signal);
%         Signal = ECGf;
%         %             set(gcf,'Pointer','arrow')
%     end
%
%     if LPfilter == 1
%         set(gcf,'Pointer','watch')
%         if n > 8
%             n = 8
%         end;
%         [b,a] = butter(n,cutoff_freq/(SR/2),'low');
%         ECGf = filtfilt(b,a,Signal_Norm);
%         Signal = ECGf;
%         LPfilter = 0;
%         cutoff_freq=cutoff_freq*.9;
%         set(gcf,'Pointer','arrow')
%     end
%
%     XL = xlim(ha);
%     XL2 = xlim(ha2);
%     % XL = xlim;
%
% end
% savefig('Peak Detection of chosen signal')
% close(gcf)
% %% Add good beats or remove bad beats
% closefig = get(0,'CurrentFigure');
% close(closefig);
%
% iRRi = Peaksmsa(Signal,points,thresh);
% bad = find(diff(t(iRRi))<0.25);
% iRRi(bad+1) = [];
%
% button = 0;
% bad = [];
% iRRiold = iRRi;
%
% tRRi = t(iRRi);
% RRi = diff(tRRi)*1000;
% HR = 60000./RRi;
% timeP = t(iRRi)';
%
% XL = [0 max(t)];
% XL2 = XL;
%
% deff = max(Signal) - max(BPdetect);
% if abs(deff) < .25
%     Mul = 1;
% elseif abs(deff) > .25
%     Mul = max(Signal)/2;
% end
%
%
% com_signal=BPdetect;
%
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% from new clean data_nasa
% %
% Remove_Group =0;
% notes = 0;
% removed_peaks = iRRi;
% iRRi_original = iRRi;
% notes=0;
% Remove_Group = 0;
% button = 0;
%
% figure('Name', 'PeakDetectionofchosensignal','NumberTitle','on','NextPlot', 'add')
%
% tsig=t;
%
% while button == 0
%
%
%     if length(t)> length(Signal)
%         tsig  = t(1:length(Signal));
%     end
%     if length(com_signal)> length(tsig)
%         com_signal = com_signal(1:length(tsig));
%     end
%
%     clf
%     %     scrsz = get(0,'ScreenSize');
%     %     set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
%
%     subplot(3,1,1:2)
%     plot(tsig ,Signal)
%     hold on
%     plot(tsig,com_signal*Mul,'k-')
%     plot(tsig(iRRi),Signal(iRRi),'ro','linewidth',2 )
%     hold off
%     scrollplot
%
%     xlim(XL)
%     ha = get(gcf,'CurrentAxes');
%
%     subplot(3,1,3)
%     plot(t(iRRi(2:end)),HR)
%     ylabel('HR')
%     scrollplot;
%
%     xlim(XL2)
%     ha2 = get(gcf,'CurrentAxes');
%     ss = get(gcf,'Position');
%
%     pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove','Position',[ss(3)*.925 ss(4)*.7 ss(3)*.05 ss(4)*.05],'Callback', ...
%         '[x, y] = ginput(1); x = x(end); [mpt bad] = min(abs(t(iRRi)-x)); iRRiold = iRRi; iRRi(bad) = []; HR = 1./diff(t(iRRi)).*60; uiresume;');
%     pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','Position',[ss(3)*.925 ss(4)*.5 ss(3)*.05 ss(4)*.05],'Callback', ...
%         '[x, y] = ginput(1); x = round(x(end).*SR); rows2 = find(x > iRRi_original); iRRiold = iRRi; iRRi(end+1) = iRRi_original(rows2(end)); iRRi = sort(iRRi); iRRi = unique(iRRi); HR = 1./diff(t(iRRi)).*60; uiresume;');
%     %     pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','Position',[ss(3)*.925 ss(4)*.5 ss(3)*.05 ss(4)*.05],'Callback', ...
%     %         '[x, y] = ginput(1); x = round(x(end).*SR); [mpt good] = max(Signal(x-round(SR/5):x+round(SR/5))); iRRiold = iRRi; iRRi(end+1) = (good + x-round(SR/5)); iRRi = sort(iRRi); iRRi = unique(iRRi); HR = 1./diff(t(iRRi)).*60; uiresume;');
%     pb3u = uicontrol(gcf,'Style','pushbutton','String','Undo','Position',[ss(3)*.925 ss(4)*.3 ss(3)*.05 ss(4)*.05],'Callback', ...
%         'iRRi = iRRiold; HR = 1./diff(t(iRRi)).*60; uiresume;');
%     pb4u = uicontrol(gcf,'Style','pushbutton','String','AddClick','Position',[ss(3)*.925 ss(4)*.9 ss(3)*.05 ss(4)*.05],'Callback', ...
%         '[x, y] = ginput(1); rows2 = find(x < t); iRRiold = iRRi; iRRi(end+1) = rows2(1); iRRi = sort(iRRi); iRRi = unique(iRRi); HR = 1./diff(t(iRRi)).*60;  uiresume;');
%     %    pb4u = uicontrol(gcf,'Style','pushbutton','String','Zoom Out','Position',[ss(3)*.4 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*1.25; uiresume;');
%     %    pb4l = uicontrol(gcf,'Style','pushbutton','String','Zoom In','Position',[ss(3)*.2 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*0.75; uiresume;');
%     %    pb5u = uicontrol(gcf,'Style','pushbutton','String','Shift Right','Position',[ss(3)*.8 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL+diff(XL)*.5; uiresume;');
%     pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove Group','Position',[ss(3)*.025 ss(4)*.7 ss(3)*.1 ss(4)*.05],'Callback','Remove_Group = 1; uiresume;');
%     pb6u = uicontrol(gcf,'Style','pushbutton','String','Undo All','Position',[ss(3)*.025 ss(4)*.3 ss(3)*.1 ss(4)*.05],'Callback', ...
%         'iRRi = iRRi_original; HR = 1./diff(t(iRRi)).*60; uiresume;');
%
%     pb6u = uicontrol(gcf,'Style','pushbutton','String','Notes','FontSize',12,'Position',[ss(3)*.025 ss(4)*.55 ss(3)*.1 ss(4)*.05],'BackgroundColor','g','Callback', ...
%         'notes = 1; uiresume;');
%
%
%
%     if exist('MCA1')
%         pbA1 = uicontrol(gcf,'Style','pushbutton','String','MCA1','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','com_signal=BV(:,1)/mean(MCA1); uiresume;');
%     end
%     if exist('MCA2')
%         pbA2 = uicontrol(gcf,'Style','pushbutton','String','MCA1','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','com_signal=BV(:,2)/mean(MCA2); uiresume;');
%     end
%     if exist('MCA')
%         pbA = uicontrol(gcf,'Style','pushbutton','String','MCA','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','com_signal=MCA/mean(MCA); uiresume;');
%     end
%
%
%
%     pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','Position',[ss(3)*.925 ss(4)*.1 ss(3)*.05 ss(4)*.05],'Callback','button = 1; uiresume;');
%
%     uiwait(gcf)
%     if notes ==1
%         answer_notes = inputdlg('please enter your notes');
%         Notes_list{length(Notes_list)+1}= answer_notes;
%         notes =0;
%     end
%
%     if Remove_Group > 0
%         iRRiold = iRRi;
%         removed_peaks =iRRi;
%         k = waitforbuttonpress;
%         point1 = get(gca,'CurrentPoint');    % button down detected
%         finalRect = rbbox;                   % return figure units
%         pointime = get(gca,'CurrentPoint');    % button up detected
%         point1 = point1(1,1:2);              % extract x and y
%         pointime = pointime(1,1:2);
%         p1 = min(point1,pointime);             % calculate locations
%         offset = abs(point1-pointime);         % and dimensions
%         x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)]; % is the time
%         y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
%         hold on
%         axis manual
%         fig = plot(x,y,'r','linewidth',2) ;
%
%         Index_Start_Array=find(t(iRRi)<= min(x));
%         Index_end_Array  = find(t(iRRi)<= max(x));
%
%         if  (isempty(Index_Start_Array)) ||  (isempty(Index_end_Array) )
%
%             display('');
%
%         else
%             Index_Start = Index_Start_Array(end);
%             Index_end = Index_end_Array(end);
%
%             removed_peaks(Index_Start:Index_end)=[];
%             iRRi = removed_peaks;
%
%             HR = 1./diff(t(iRRi)).*60;
%             Remove_Group = 0;
%             Index_Start_Array=[];
%             Index_end_Array=[];
%             removed_peaks=[];
%         end
%     end
%
%
%     XL = xlim(ha);
%     XL2 = xlim(ha2);
%
% end
% savefig('Peak Detection of chosen signal')
% close(gcf)
% ikp = find(iRRi<=length(t));
% iRRi = iRRi(ikp);
% tRRi = t(iRRi);
% RRi = diff(tRRi)*1000;
%
% tRRi(end) = [];
% time = t;
%

%% Peak Detectection for Velocity_PeakDetect
% BP = Velocity_PeakDetect;
% timeMat = min(timeMatraw):timeMatdiff:max(timeMatraw);
if strcmp(Velocitydlg,'MATLAB')
    
    m = 1; %value preset for use in basefilter (Highpass) function.
    n = 1; %value preset for use in LPF function.
    timeMat = timeMarker;
    timeMatSound = timeMat; % Created new variable for the purpose of being able to select Velocity Sound Signal during loop.
    tRRi_corr = tRRi;
    % timeMat = timeMat;
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
    
    figure('Name','PeakDetectectionforVelocity_PeakDetect','NumberTitle','on','NextPlot', 'add')
    
    while button==0
        %     Velocity_PeakDetect = Velocity_PeakDetect;
        
        if redetect == 1
            iRRiV = Peaksmsa(Velocity_PeakDetect,points,thresh);
            bad = find(diff(timeMat(iRRiV))<0.15);
            iRRiV(bad+1) = [];
            HRV = 1./diff(timeMat(iRRiV))*60;
            XL = [0 max(timeMat)];
            XL2 = [0 max(timeMat)];
            redetect = 0;
        end
        
        clf('reset')
        subplot(3,1,1:2)
        plot(timeMat,Velocity_PeakDetect)
        ha = get(gcf,'CurrentAxes');
        hold on
        hline = refline(0,thresh);
        set(hline,'Color','g');
        plot(timeMat(iRRiV),Velocity_PeakDetect(iRRiV),'ro','linewidth',2 )
        hold off
        %scrollplot
        title ('Peak Detection of Velocity or Diameter Signal');
        
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
        
        %     pb10 = uicontrol(gcf,'Style','pushbutton','String','ECG inverse','units','normalized','Position',[.2 .025 .05 .03],'Callback','ECG = -ECG; uiresume;');
        %         pb7 = uicontrol(gcf,'Style','pushbutton','String','Left Artery','units','normalized','Position',[.4 .025 .05 .03],'Callback','ECG = BV(:,1); thresh = max(ECG)*0.35; uiresume;');
        %         pb8 = uicontrol(gcf,'Style','pushbutton','String','Right Artery','units','normalized','Position',[.6 .025 .05 .03],'Callback','ECG = BV(:,2); thresh = max(ECG)*0.35; uiresume;');
        %         pb9 = uicontrol(gcf,'Style','pushbutton','String','MCA','units','normalized','Position',[.8 .025 .05 .03],'Callback','ECG = MCA;  thresh = max(ECG)*0.35; uiresume;');
        
        pb9 = uicontrol(gcf,'Style','pushbutton','String','Diameter','units','normalized','Position',[.8 .025 .05 .03],'Callback','Velocity_PeakDetect = diamBA; timeMat = timeBAvel; thresh = max(diamBA)*0.35; redetect = 1; uiresume;');
        pb7 = uicontrol(gcf,'Style','pushbutton','String','Velocity_Sound','units','normalized','Position',[.4 .025 .05 .03],'Callback','Velocity_PeakDetect = VelMat; timeMat = timeMatSound; thresh = max(VelMat)*0.35;  redetect = 1;uiresume;');
        pb8 = uicontrol(gcf,'Style','pushbutton','String','Velocity_Brachial Analyzer','units','normalized','Position',[.6 .025 .05 .03],'Callback','Velocity_PeakDetect = Velocity_BA; timeMat = timeBAvel; thresh = max(Velocity_BA)*0.35; redetect = 1; uiresume;');
        
        pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
        
        uiwait(gcf)
        
        Velocityf = zeros(length(Velocity_PeakDetect),1);
        
        if basefilter == 1
            %            set(gcf,'Pointer','watch')
            %         baseECG = medfilt1(ECG,round(SR));
            %         ECGf=ECG - baseECG + mean(baseECG);
            %       use a 5 sec moving median to remove baseline changes
            
            %         Velocityold = VelocitSelect;
            %             for i = 1:length(ECG)
            %                 if i+round(SR)*3 <= length(ECG)
            %                     ECGf(i) = median(ECG(i:i+round(SR)*3));
            %                 else
            %                     ECGf(i) = median(ECG(length(ECG)-(round(SR)*3):length(ECG)));
            %                 end
            %             end
            %             szECGf = size(ECGf);
            %             if szECGf(1)>1
            %                 ECGf = ECGf.';
            %             end
            %             ECG = ECG - ECGf;
            basefilter = 0;
            if m > 5
                m = 5;
            end
            [b,a] = butter(m,5/SR,'high');
            Velocityf = filtfilt(b,a,Velocity_PeakDetect);
            Velocity_PeakDetect = Velocityf;
            %             set(gcf,'Pointer','arrow')
        end
        
        if LPfilter == 1
            set(gcf,'Pointer','watch')
            if n > 8
                n = 8
            end;
            [b,a] = butter(n,cutoff_freq/(SR/2),'low');
            Velocityf = filtfilt(b,a,Velocity_PeakDetect);
            Velocity_PeakDetect = Velocityf;
            LPfilter = 0;
            cutoff_freq=cutoff_freq*.9;
            set(gcf,'Pointer','arrow')
        end
        
        XL = xlim(ha);
        
    end
    %    savefig('Peak Detectection for Velocity_PeakDetect')
    close(gcf)
    % Section 9A: Add good beats, remove bad beats
    closefig = get(0,'CurrentFigure');
    close(closefig);
    
    % iRRiV = Peaksmsa(Velocity_PeakDetect,points,thresh);
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
    
    % deff = max(Velocity_PeakDetect) - max(MCAdetect);
    % if abs(deff) < .25
    %     Mul = 1;
    % elseif abs(deff) > .25
    %     Mul = max(Velocity_PeakDetect)/2;
    % end
    
    
    % com_signal=MCAdetect;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% feom ncleandata_nasa
    %
    Remove_Group =0;
    notes = 0;
    removed_peaks = iRRiV;
    iRRiV_original = iRRiV;
    notes=0;
    Remove_Group = 0;
    button = 0;
    SR = round(1./mean(diff(timeMat)));
    figure(  'Name','PeakDetectectionforVelocity_PeakDetect','NumberTitle','on','NextPlot', 'add')
    
    while button == 0
        
        
        %     if length(com_signal)> length(timeMat)
        %         com_signal = com_signal(1:length(timeMat));
        %     end
        
        
        clf
        %     scrsz = get(0,'ScreenSize');
        %     set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
        clf
        subplot(3,1,1:2)
        plot(timeMat,Velocity_PeakDetect)
        hold on
        %     plot(timeMat,com_signal*Mul,'k-')
        plot(timeMat(iRRiV),Velocity_PeakDetect(iRRiV),'ro','linewidth',2 )
        hold off
        %scrollplot
        xlim(XL)
        ha = get(gcf,'CurrentAxes');
        subplot(3,1,3)
        plot(timeMat(iRRiV(2:end)),HRV)
        ylabel('DVD RRI')
        %scrollplot;
        xlim(XL2)
        ha2 = get(gcf,'CurrentAxes');
        ss = get(gcf,'Position');
        
        pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove','Position',[ss(3)*.925 ss(4)*.7 ss(3)*.05 ss(4)*.05],'Callback', ...
            '[x, y] = ginput(1); x = x(end); [mpt bad] = min(abs(timeMat(iRRiV)-x)); iRRiVold = iRRiV; iRRiV(bad) = []; HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
        pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','Position',[ss(3)*.925 ss(4)*.5 ss(3)*.05 ss(4)*.05],'Callback', ...
            '[x, y] = ginput(1); x = round(x(end).*SR); rows2 = find(x > iRRiV_original); iRRiVold = iRRiV; iRRiV(end+1) = iRRiV_original(rows2(end)); iRRiV = sort(iRRiV); iRRiV = unique(iRRiV); HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
        % '[x, y] = ginput(1); x = round(x(end).*SR); [mpt good] = max(Velocity_PeakDetect(x-round(SR/5):x+round(SR/5))); iRRiVold = iRRiV; iRRiV(end+1) = (good + x-round(SR/5)); iRRiV = sort(iRRiV); iRRiV = unique(iRRiV); HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
        pb3u = uicontrol(gcf,'Style','pushbutton','String','Undo','Position',[ss(3)*.925 ss(4)*.3 ss(3)*.05 ss(4)*.05],'Callback', ...
            'iRRiV = iRRiVold; HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
        pb4u = uicontrol(gcf,'Style','pushbutton','String','AddClick','Position',[ss(3)*.925 ss(4)*.9 ss(3)*.05 ss(4)*.05],'Callback', ...
            '[x, y] = ginput(1); rows2 = find(x < timeMat); iRRiVold = iRRiV; iRRiV(end+1) = rows2(1); iRRiV = sort(iRRiV); iRRiV = unique(iRRiV); HRV = 1./diff(timeMat(iRRiV)).*60;  uiresume;');
        %         '[x, y] = ginput(1); x = round(x(end).*SR); iRRiVold = iRRiV; iRRiV(end+1) = x; iRRiV = sort(iRRiV); iRRiV = unique(iRRiV); HRV = 1./diff(timeMat(iRRiV)).*60;  uiresume;');
        %    pb4u = uicontrol(gcf,'Style','pushbutton','String','Zoom Out','Position',[ss(3)*.4 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*1.25; uiresume;');
        %    pb4l = uicontrol(gcf,'Style','pushbutton','String','Zoom In','Position',[ss(3)*.2 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*0.75; uiresume;');
        %    pb5u = uicontrol(gcf,'Style','pushbutton','String','Shift Right','Position',[ss(3)*.8 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL+diff(XL)*.5; uiresume;');
        pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove Group','Position',[ss(3)*.025 ss(4)*.7 ss(3)*.1 ss(4)*.05],'Callback','Remove_Group = 1; uiresume;');
        pb6u = uicontrol(gcf,'Style','pushbutton','String','Undo All','Position',[ss(3)*.025 ss(4)*.3 ss(3)*.1 ss(4)*.05],'Callback', ...
            'iRRiV = iRRiV_original; HRV = 1./diff(timeMat(iRRiV)).*60; uiresume;');
        
        pb6u = uicontrol(gcf,'Style','pushbutton','String','Notes','FontSize',12,'Position',[ss(3)*.025 ss(4)*.55 ss(3)*.1 ss(4)*.05],'BackgroundColor','g','Callback', ...
            'notes = 1; uiresume;');
        
        
        %     if exist('MCA1')
        %         pbA1 = uicontrol(gcf,'Style','pushbutton','String','MCA1','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','com_signal=BV(:,1)/mean(MCA1); uiresume;');
        %     end
        %     if exist('MCA2')
        %         pbA2 = uicontrol(gcf,'Style','pushbutton','String','MCA1','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','com_signal=BV(:,2)/mean(MCA2); uiresume;');
        %     end
        %     if exist('MCA')
        %         pbA = uicontrol(gcf,'Style','pushbutton','String','MCA','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','com_signal=MCA/mean(MCA); uiresume;');
        %     end
        
        
        
        pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','Position',[ss(3)*.925 ss(4)*.1 ss(3)*.05 ss(4)*.05],'Callback','button = 1; uiresume;');
        
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
    %    savefig('Peak Detectection for Velocity_PeakDetect')
    close(gcf)
    ikpv = find(iRRiV<length(timeMat));
    tRRiV = timeMat(iRRiV);
    RRiV = diff(tRRiV)*1000;
    
    %tRRiV(end) = [];
    % timeVel = t;
    
    % Section 10 : Plot R-R Intervals for Signal data and Matlab derived velocity data to overlap
    
    %Calculate HR from properly peak detected ECG and Velocity signals
    %HR = 1./diff(t(iRRi))*60;
    % HRV = 1./diff(timeMat(iRRiV)).*60;
    HRV = 1./diff(tRRiV).*60;
    
    button = 0;
    figure(  'Name','CalculateHRfromproperlypeakdetectedECGandVelocitysignals','NumberTitle','on','NextPlot', 'add')
    %    HRnorm = zscore(HR) ;
    %    HRVnorm = zscore(HRV);
    %     while button == 0
    %
    %         clf('reset')
    %         %     scrsz = get(0,'ScreenSize');
    %         %     set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
    %         plot(t(iRRi(1:end-1)), HRnorm, 'bo-')
    %         hold on
    %         plot(timeMarker(iRRiV(1:end-1)),HRVnorm, 'ro-')
    %         hold off
    %         scrollplot;
    %         % xlim(XL);
    %         title('Velocity and (BP or ECG) RRI Plot')
    %         xlabel('Time(Seconds)')
    %         legend('LabChart Signal', 'Velocity');
    %
    %         pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
    %
    %         uiwait(gcf)
    %     end
    %    savefig('Calculate HR from properly peak detected ECG and Velocity signals')
    close(gcf)
    %
    % button = questdlg('Would you like to align the data?','Continue Operation','Yes','No','Help','No');
    % if strcmp(button,'Yes')
    %     disp('Proceding to analysis')
    
    %     uiwait( errordlg('Data needs to be aligned, Please call an Engineer.'));
    
    
    % Section 12: Align Velocity and Signal R-R Intervals
    
    % if strcmp(button,'Yes')
    %
    SR = round(1./mean(diff(t)));
    XL = [0 max(t)];
    tHR = tRRi_old;
    tHRV = timeMarker;
    HRVShift = HRV;
    
    %  uiwait(msgbox('Align the signals as close as possible'));
    %% plot heart rates to align
    
    figure(  'Name','AlignVelocityandSignalR-RIntervals','NumberTitle','on','NextPlot', 'add')
    button = 0;
    
    %     HRnorm = zscore(HR) ;
    %     HRVShiftnorm = zscore(HRVShift);
    
    YL = [min([HR HRVShift]) max([HR HRVShift])];
    
    while button == 0
        
        
        clf('reset')
        %         scrsz = get(0,'ScreenSize');
        %         set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
        plot(tHR,HR, 'bo-')
        %        plot(HR, 'bo-')
        
        %   plot(t(iRRi(2:end)),HR, 'b', t(iRRiV(2:end)),HRV, 'r')
        %   plot((timehr(inthr)-timehr(inthr(1))),hr(inthr),'bo-')
        hold on
        plot(tHRV(iRRiV(1:end-1)),HRVShift, 'ro-')
        %        plot(HRVShift, 'ro-')
        
        %scrollplot
        %   plot((time24rr(intrr)-time24rr(intrr(1))),60./rr(intrr),'ro-')
        hold off
        xlim(XL);
        ylim(YL);
        title('Velocity and (BP or ECG)RRI Plot')
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
        %         shifter
        
    end
    %     graphfile = [infile2(1:end-4),'_Align_f',num2str(gcf)];
    %     saveas(gcf,graphfile,'jpg');
    %     savefig(graphfile)
    close(gcf)
    
    %% Section 13: Further Align Velocity(or Diameter) and Signal based on previous shift
    % decided to comment out since don't think we need a second alignment
    
    button = 0;
    delay_shift=0;
    % Velocity_Shift = VelMat';
    %    Velocity_Shift = Velocity_PeakDetect;
    
    %     %zscore ECG signal
    %     maxECG = max(ECG);
    %     minECG = min(ECG);
    %     ECGw = ECG-minECG;
    %     ECGnorm = ECGw./median(ECGw);
    
    %     %zscore BP signal
    %     maxBP = max(BP);
    %     minBP = min(BP);
    %     BPw = BP-minBP;
    %     BPnorm = BPw./median(BPw);
    %     markernorm = markerChart./median(markerChart);
    
    %    diam_Shift = diamBA;
    
    %     XL = [min([min(tHR) min(tHRV)]) max([max(tHR) max(tHRV)])];
    %
    %     BPnorm = zscore(BP) ;
    %     markernorm = zscore(markerChart) ;
    %
    %     mul = max(BPnorm) / max(markernorm);
    %
    %
    
    BP_Shiftnorm = zscore(bp_f);
    Velocity_Shiftnorm = zscore(Velocity_PeakDetect);
    figure(  'Name','FurtherAlignVelocity(orDiameter)andSignalbasedonpreviousshift','NumberTitle','on','NextPlot', 'add')
    
    while button == 0
        
        
        clf
        switch Velocitydlg
        end
        plot(t, BP_Shiftnorm, 'b-')
        hold on
        %         plot(tHRV,Velocity_Shift./median(Velocity_Shift), 'r-')
        plot(tHRV,Velocity_Shiftnorm, 'r-')
        hold off
        %    axis tight
        %scrollplot;
        xlim(XL);
        title('Velocity and BP Plot')
        xlabel('Time(Seconds)')
        legend('LabChart Signal', 'marker','Velocity');
        
        ha = get(gcf,'CurrentAxes');
        %    YL1 = ylim;
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
    %     graphfile = [infile2(1:end-4),'_FurtherAlign_f',num2str(gcf)];
    %     saveas(gcf,graphfile,'jpg');
    %     savefig(graphfile)
    close(gcf)
    timeAligned = tHRV;
    % /////////////////////////////////////////////////////////////////////////
    %check for debugging purposes
    
    time_dur_vel=tHRV(end)-tHRV(1);
    
    time_dur_diam=timeBA(end)-timeBA(1);
    % /////////////////////////////////////////////////////////////////////////
    
    ShiftTime = nanmean(tHRV(1)-t(1));
    
    if ba_dia_exist==0                      %added by BS 1/13/18 since if BA diam is chosen needs to be aligned same as matlab velocity, can be further aligned later.
        %shift diameter time the same amount
        timeBA_orig = timeBA;
        timeBA = timeBA + ShiftTime;
    else
        %shift diameter time the same amount
        timeBA_orig = timeBA;
        timeBA = timeBA + ShiftTime;
    end
    
    
    button = 0;
    
    Diam_Shiftnorm = zscore(diamBA);
    figure(  'Name','FurtherAlignDiameterandSignalbasedonpreviousshift','NumberTitle','on','NextPlot', 'add')
    
    % Added 9/5/2019 //////////////////////////////////////////////////////////
    Diam_Shiftnorm = interp1(timeBA,Diam_Shiftnorm,tHRV);
    
    % /////////////////////////////////////////////////////////////////////////
    
    
    while button == 0
        
        
        clf
        switch Velocitydlg
        end
        plot(t, BP_Shiftnorm, 'b-')
        hold on
        %         plot(tHRV,Velocity_Shift./median(Velocity_Shift), 'r-')
        
        % Added 9/5/2019 //////////////////////////////////////////////////////////
        %         plot(timeBA,Diam_Shiftnorm, 'r-')
        plot(tHRV, Diam_Shiftnorm,'r')
        % /////////////////////////////////////////////////////////////////////////
        hold off
        %    axis tight
        %scrollplot;
        xlim(XL);
        title('Diameter and BP Plot')
        xlabel('Time(Seconds)')
        legend('LabChart Signal', 'marker','Velocity');
        
        ha = get(gcf,'CurrentAxes');
        %    YL1 = ylim;
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
    
    %     disp('Ending Program')
    % elseif strcmp(button,'Help')
    %     disp('Quick, call an engineer! x2178')
    % SR = round(1./mean(diff(t)));
    % XL = [0 max(t)];
    % tHR = t;
    % tHRV = timeMarker;
    % HRVShift = HRV;
elseif strcmp(Velocitydlg,'Brachial_Analyzer')
    %     SR = round(1./mean(diff(t)));
    %     XL = [0 max(t)];
    tHR = t;
    tHRV = timeMarker;
    iRRiV = Peaksmsa(Velocity_PeakDetect,points,thresh);
    
    %     HRVShift = HRV;
    disp('Continue...')
    
    %Correct Velocity angle from BA
    % Brachial Analyzer velocity correction script
    % For: Concussion analysis (velocity avg was calculated using brachial analyzer *wrong angle in video, hence needs correction)
    % Written by: Bishoy Samy
    % Date: 12/14/2015
    
    
    subset_avgflow = Velocity_PeakDetect;
    subset_time =   timeMarker;
    
    prompt = {'Enter angle from P300:','Enter correct angle:'};
    dlg_title = 'Input angles';
    num_lines = 1;
    defaultans = {'00','00'};
    angles = inputdlg(prompt,dlg_title,num_lines,defaultans);
    angles=str2double(angles);
    correct_angle= angles(2);
    p300_angle= angles(1);
    
    corrected = subset_avgflow.*(cosd(p300_angle)/cosd(correct_angle));
    %corrected2 = subset_avgflow./cosd(abs(correct_angle)-abs(p300_angle));
    
    % corrected_velmax = subset_vel_max.*(cosd(p300_angle)/cosd(correct_angle));
    
    % figure; plot(subset_time,subset_avgflow,'r');xlabel('Time (s)');ylabel('Avg Flow (cm/s)');hold on; plot(subset_time,corrected,'g'); ...
    %     %plot(subset_time,subset_vel_max,'ro'),hold on; plot(subset_time,corrected_velmax,'go'); ...        %corrected max velocities from excel sheet
    %     legend(['Original_' num2str(p300_angle)],['Corrected_' num2str(correct_angle)],-1);hold on; title('Mean velocity - Original vs corrected');%hold on; plot(subset_time,corrected2,'g');
    % scrollplot;
    
    Velocity_PeakDetect = corrected;
    
end
%% zscore all signals - do not zscore signals (Dr Serrador)
%button = questdlg('Are both Signals aligned?', 'Continue Operation','Yes','No','Help');

button = 'Yes';

% if strcmp(button,'Yes')
%     disp('Pick points for alignment')
%     disp('Ending alignment')
% elseif strcmp(button,'No')
%     disp('Ending alignment')
% elseif strcmp(button,'Help')
%     disp('Quick! call the engineers')
%     pause
% end

if strcmp(button,'Yes')
    
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
        bv_fo_cm = BVcm_int(intwf1,:);      %added by BS 2/16/18 in order to include tcd cm/s
        %     end
    else
        bv_f = BVcm_int;
        bv_fo = BVcm_int;
        bv_fo_cm = BVcm_int;
    end
    clear i
    
    intwf2 = find(time>=startbin & time<=stopbin);
    tECGo = time(intwf2);
    ECGo = ECG;
    %     ECGo = ECG(intwf2);
    
    intRRi_temp = find(iRRi>=intwf2(1) & iRRi<=intwf2(end));
    iRRio = iRRi(intRRi_temp);
    iRRio = iRRio - intwf2(1)+1;
    
    %trim peak detection times to area of interest
    %     intRR1 = find(iRRi>=intwf2(1) & iRRi<=intwf2(end));
    %     iRRio = iRRi(intRR1)-intwf2(1)+1;
    
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
    
    
    %correct location of peaks by removing number of points removed from
    %iRRi value
    %     iRRio = iRRi - int(1);
    %     goodiRRi = find(iRRio>0 & iRRio<=length(int));
    %     iRRio = iRRio(goodiRRi);
    
    % Replace Velocity signal with actual velocity signal to move forward
    
    %     if eq(length(Velocity_PeakDetect),length(VelMat))
    %         Velocity_Shift = VelMat;
    %
    %     elseif eq(length(Velocity_PeakDetect),length(Velocity_BA))
    %         Velocity_Shift = Velocity_BA;
    %
    %
    %     end
    Velocity_Shift = Velocity_PeakDetect;
    
    int2 = find(tHRV>=startbin & tHRV<=stopbin);
    tHRVo = tHRV(int2);
    MatVelo = Velocity_Shift(int2);
    
    int5 = find(timeBA>=startbin & timeBA<=stopbin);
    timeBAo = timeBA(int5);
    DiamBAo = diamBA(int5);
    
    int3 = find(tRRi>=startbin & tRRi<=stopbin);
    tRRio = tRRi(int3);
    %     HRo = HR(int3);
    iRRio = iRRi(int3);
    
    tRRiV = timeMat(iRRiV);
    int4 = find(tRRiV>=startbin & tRRiV<=stopbin);
    iRRiVo = iRRiV(int4);
    tRRiVo = tRRiV(int4);
    
    %     %adjust times so all are relative to startbin
    %     tHRo = tHRo-startbin;
    %     tHRVo = tHRVo - startbin;
    %     tRRio = tRRio - startbin;
    %     tRRiVo = tRRiVo - startbin;
    %     to = to - startbin;
    %     tco2o = tco2o - startbin;
    %     tetCO2 = tetCO2 - startbin;
    
    endaligntime = max([max(tHRo) max(tHRVo) max(tRRio) max(tRRiVo) max(to) max(tco2o) max(tetCO2)]);
    
    
end

if strcmp(button,'No')
    uiwait( errordlg('Signals are not aligned, Please try again.'));
end

% %% Section 15: Remove velocity spike if there
%
% MatVelo = Velocity_Shift(int2);
% % SpikeRemove(MatVelo,tHRVo);
%
%
% mkrbutton = questdlg('Do you want to remove marker spike?', 'Continue Operation','Yes','No','Help');
%
% if strcmp(mkrbutton,'Yes')
%
%     button = 0;
%     XL = [tHRVo(1) tHRVo(end)];
%
%     figure('Position',screensize(1), 'Name','Remove velocity spike if there','NumberTitle','on','NextPlot', 'add')
%     %   clear spikeval;
%     while button == 0
%
%
%         clf
%         %         scrsz = get(0,'ScreenSize');
%         %         set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
%         plot(tHRVo,MatVelo,'b-')
%
%         if exist('spikeval','var')
%             hold on
%             plot(tHRVo(spikeval),MatVelo(spikeval),'r-')
%             hold off
%         end
%
%         %       scrollplot
%         xlim(XL);
%         title('Velocity')
%         xlabel('Time(Seconds)')
%
%         ha = get(gcf,'CurrentAxes');
%         YL1 = ylim;
%
%         ss = get(gcf,'Position'); scrollplot;
%         %         pb1u = uicontrol(gcf,'Style','slider','Max',round(max(tHRVo)),'Min',[0],'Value',round(XL(2)),'SliderStep',[0.01 0.1],'units','normalized','Position',[.925 .6 .025 .3],'Callback','XL(2) = get(pb1u,Value); xlim(XL); uiresume;');
%         %         pb1ut1 = uicontrol(gcf,'style','text','string',ceil(XL(2)),'units','normalized','HorizontalAlignment','left','Position',[.96 .9 .025 .025]);
%         %         pb1ut2 = uicontrol(gcf,'style','text','string',floor(XL(1)),'units','normalized','HorizontalAlignment','left','Position',[.96 .6 .025 .025]);
%         %
%         %         pb1u2 = uicontrol(gcf,'Style','slider','Max',round(max(tHRVo)),'Min',[0],'Value',round(XL(1)),'SliderStep',[0.01 0.1],'units','normalized','Position',[.925 .1 .025 .3],'Callback','XL(1) = get(pb1u2,Value); xlim(XL); uiresume;');
%         %         pb1ut3 = uicontrol(gcf,'style','text','string',ceil(XL(2)),'units','normalized','HorizontalAlignment','left','Position',[.96 .4 .025 .025]);
%         %         pb1ut4 = uicontrol(gcf,'style','text','string',floor(XL(1)),'units','normalized','HorizontalAlignment','left','Position',[.96 .1 .025 .025]);
%
%         pb1au = uicontrol(gcf,'Style','pushbutton','String','Select','units','normalized','Position',[.925 .5 .05 .05],'Callback','rect = getrect; uiresume;');
%
%         pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.95 .05 .05 .05],'Callback','button = 1; uiresume;');
%
%         uiwait(gcf)
%
%         if exist('rect','var')
%             spikeval = find(tHRVo>=rect(1) & tHRVo<=(rect(1)+rect(3)));
%             clear rect;
%         end
%
%         XL = xlim(ha);
%
%         %replace spike with linear interpolation
%
%
%     end
%
%     MatVelo_clean = MatVelo;
%     MatVelo_clean(spikeval) = [];
%     tHRVo_clean = tHRVo;
%     tHRVo_clean(spikeval) = [];
%     MatVelo_int = interp1(tHRVo_clean,MatVelo_clean,tHRVo,'linear');
%     MatVelo_orig = MatVelo;
%     MatVelo = MatVelo_int;
%
% elseif strcmp(mkrbutton,'No')
%     disp('Pick points for alignment')
% elseif strcmp(mkrbutton,'Help')
%     disp('Quick! call the engineers')
%     pause
% end

%  elseif strcmp(button,'No')
%      disp('Contnuing Analysis')
% end

%% Plotting all the channels

figure('Name','PlottingalltheSignals','NumberTitle','on','NextPlot', 'add')
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
    if isempty(ticktimes(comtick))==0
        for txt=1:length(markers)
            tetx_handle = text (t(markers_Index(txt)) ,bp_f(markers_Index(txt)),markers{txt});
            set(tetx_handle,'Colo','r','FontSize',txtfontsize,'Rotation',90)
        end
    end
    axis tight
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    %    xlabel('Time (s)', fs, fs_label)
    ylabel({'Arterial';'Blood Pressure';'(mmHg)'}, fs, fs_label)
    title(plot_title, fs, fs_title)
    
    for figcont = 1:num_artery
        subplot(k,1,figcont+1)
        plot(to,bv_fo(:,figcont), 'g', lw, lw_thin), hold on;
        plot(tRRio,MBVo(:,figcont),'k', lw, lw_thick),
        plot(tRRio,SBVo(:,figcont),'r', lw, lw_thick),
        plot(tRRio,DBVo(:,figcont),'r', lw, lw_thick),
        line([to(1) to(end)], [100 100], 'Color','m'),hold off;
        axis tight
        if isempty(ticktimes(comtick))==0
            for txt=1:length(markers)
                tetx_handle = text (t(markers_Index(txt)) ,bv_f(markers_Index(txt),figcont),markers{txt});
                set(tetx_handle,'Colo','r','FontSize',txtfontsize,'Rotation',90)
            end
        end
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
    %    xlabel('Time (s)', fs, fs_label);
    
    subplot(k,1,num_artery+3)
    plot(timeBAo,DiamBAo,'r', lw, lw_thick);
    axis tight
    meanDiamBAo= mean(DiamBAo((round(length(DiamBAo)/4)):((round(length(DiamBAo)/4))*3)));
    stdDiamBAo = std(DiamBAo((round(length(DiamBAo)/4)):((round(length(DiamBAo)/4))*3)));
    yl(1) = meanDiamBAo-stdDiamBAo*4;
    yl(2) = meanDiamBAo+stdDiamBAo*4;
    ylim(yl); xlim([tHRVo(1) tHRVo(end)]);
    ylabel({'Diameter';'(cm)'}, fs, fs_label)
    %    xlabel('Time (s)', fs, fs_label);
    
    subplot(k,1,num_artery+4)
    plot(tHRo,HRo,'r', lw, lw_thick);
    axis tight
    yl = ylim;
    yl(1) = yl(1)-yl(2)*.05;
    yl(2) = yl(2)+yl(2)*.05;
    ylim(yl);
    ylabel({'HR';'(bpm)'}, fs, fs_label)
    %    xlabel('Time (s)', fs, fs_label);
    
    subplot(k,1,num_artery+5)
    plot(tco2o,CO2o, 'g', lw, lw_thin), hold on;
    plot(tetCO2o,etCO2o,'r', lw, lw_thick), hold off;
    axis tight
    if isempty(ticktimes(comtick))==0
        for txt=1:length(markers)
            tetx_handle = text (tco2(markers_Index(txt)) ,CO2(markers_Index(txt)),markers{txt});
            set(tetx_handle,'Colo','r','FontSize',txtfontsize,'Rotation',90)
        end
    end
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
%graphfile = [infile2(1:end-4),'_AllSignals_f',num2str(gcf)];
%saveas(gcf,graphfile,'jpg');
%savefig(graphfile)
%close(gcf)

%% Select region of interest

figure(  'Name','Selectregionofinterest','NumberTitle','on','NextPlot', 'add')
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
    %         ylabel(Y_Axis{j});
    %         Y_Axis({Artery_to_use})
    %scrollplot
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
    %         ylabel(Y_Axis{j});
    %         Y_Axis({Artery_to_use})
    %scrollplot
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
    %     pb5u = uicontrol(gcf,'Style','pushbutton','String','Diameter LPF','units','normalized','Position',[.05 .15 .05 .05],'Callback','n = n + 1; LPfilter = 1; redetect = 1; uiresume;');
    
    uiwait(gcf)
    
    %     if LPfilter == 1
    %         set(gcf,'Pointer','watch')
    % %         [b,a] = butter(8,cutoff_freq/n,'low');
    %         DIAMoriginal = newdata_Select(3,:);
    % %         Diamfilt = filtfilt(b,a,DIAMoriginal);
    %         Diamfilt = filtfilt(lpFilt,DIAMoriginal);
    %         DIAM = Diamfilt;
    %         LPfilter = 0;
    %         cutoff_freq=cutoff_freq*.9;
    %         set(gcf,'Pointer','arrow')
    %     end
    
end
close(gcf)

block_data_start = start;
block_data_stop = stop;
startime = start(1);
endtime = stop(1);

% SR = round(1./mean(diff(t)));
% points = ceil(0.2*SR);
% Signal2 = newdata_Select(1,:);
% thresh = max(Signal2)*0.15;
% iRRii = Peaksmsa(Signal2,points,thresh);

%correct peak detection points again based on section
% good2 = find(t(iRRii)>=startime & t(iRRii)<=endtime);
% iRRii_old = iRRii;
% iRRii = iRRii(good2);
% iRRii = iRRii - (good(1)-1);

%% save old interpolated data _old variable and replace with new overlapped
%interopolated data
%ECGi_old = ECGi;
BPi_old = BPi;
BVi_old = BVi;
CO2i_old = CO2i;
clear BPi BVi CO2i

%resample everything to 100 Hz for cleaning and further analysis
tint = startime:0.01:endtime;
if length(tECGo)~=length(ECGo);             %added if statement since if data was collected at higher sampling caused error
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


%iRRii = round(iRRio./(SR/100));

% %% choosing which signal the user would Like to clean
% Artery_to_remove=[];
% choose_artries=0;
% plotting_color= ['r' 'b' 'k' 'r' 'k'  ];
% %     Artery_to_use= ones(num_artery,1);
% n=2; %Added because we only need to clean Velocity and Diameter
% Artery_to_use= ones(n,1);
%
% %Create new matrix for data to be selected
% p = length(tint);
% newdata_Select = zeros(n,p);
% %newdata_Select(1,:)= BPi;  %Assign 1st channel to BP
% newdata_Select(1,:)= MatVeli; %Assign 1st channel to Velocity
% newdata_Select(2,:)= DiamBAi; %Assign 2nd channel to Diameter
%
% % %Assign the appropriate titles to the selected signals
% % titles_Select = cell(n,1);
% %
% % %Select the right title for the BP signal going forward
% % titles_str = cellstr(titles);
% % for i = 1: length(titles_str)
% %     if strcmp (titles_str(i,:),'BP')
% %         titles_Select(1) = titles_str(i);
% %     end
% % end
%
% titles_Select(1) = {'Velocity'};
% titles_Select(2) = {'Diameter'};
%
% %Assign the appropriate units to the selected signals
% Units_Select = cell(n,1);
% Units_Select(1) = cellstr('cm/sec');
% Units_Select(2) = cellstr('mm');
%
% t=tint;
%
% figure(  'Name','choosingwhichsignaltheuserwouldLiketoclean','NumberTitle','on','NextPlot', 'add')
%
% num_artery = n;
% while choose_artries <1
%
%     for nm=1:num_artery
%         subplot(num_artery,1,nm)
%         plot(t,newdata_Select(nm,:),plotting_color(nm))
%
%
%         ylabel(titles_Select(nm) )
%         title('Please check Only the signal that you would like to Clean ');
%
%         BU(nm) = uicontrol(gcf,'Style','checkbox','String',titles_Select(nm),'BackgroundColor','c','FontSize',10,'units','normalized','Position',[.925 (1.1/num_artery)*(num_artery-nm) .1 .05],'Callback', ' uiresume;');
%
%         set(BU(nm),'value',Artery_to_use(nm));
%     end
%
%     b4u = uicontrol(gcf,'Style','pushbutton','String',' Finish ','BackgroundColor','r','FontSize',15,'units','normalized','Position',[.025 .05 .05 .05],'Callback', 'choose_artries =2 ; uiresume;');
%     uiwait(gcf)
%     for chkboxnum=1:(num_artery)
%         Artery_to_use(chkboxnum) =get(BU(chkboxnum),'value');
%     end
% end
% savefig('choosing which signal the user would Like to clean')
% close(gcf)
%
% a=1;
%
% for all = 1:length(Artery_to_use);
%
%     if Artery_to_use(all) < 1
%         Artery_to_remove (a) = all;
%         a=a+1;
%     end
% end
%
% num_artery=num_artery-length(Artery_to_remove);
% newdata_Select(Artery_to_remove,:)=[];
% titles_cell(Artery_to_remove)=[];

%% start cleaning of velocity and diameter

% Section 19: Plot parameters
%numchannels = ch;
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

n=1;

while button==0
    
    clf
    
    plot(BVfilt);
    %scrollplot;
    
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
    %scrollplot;
    
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
    figure(99)
    subplot(2,2,1)
    plot(tint(startptplot:endptplot),MatVeli(startptplot:endptplot))
    hold on
    plot(tint(SVeloPt(i)),MatVeli(SVeloPt(i)),'go')
    plot(tint(DVeloPt(i)),MatVeli(DVeloPt(i)),'ro')
    hold off
    xlim([tint(startptplot) tint(endptplot)])
    
    %        pause
    subplot(2,2,3)
    plot(tint(startptplot:endptplot),DiamBAi(startptplot:endptplot))
    hold on
    plot(tint(SDiamPto(i)),DiamBAi(SDiamPto(i)),'go')
    plot(tint(DDiamPto(i)),DiamBAi(DDiamPto(i)),'ro')
    hold off
    xlim([tint(startptplot) tint(endptplot)])
    
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
        
        subplot(2,2,4)
        hold on
        plot(DiamBAiPlot(i,:))
        plot((SDiamPto(i)-(intbeat(1)-30)),DiamBAiPlot(i-1,(SDiamPto(i)-(intbeat(1)-30))),'go')
        plot((DDiamPto(i)-(intbeat(1)-30)),DiamBAiPlot(i-1,(DDiamPto(i)-(intbeat(1)-30))),'ro')
    end
    
    
    if exist('BV1')==1
        [LIA,LOCB] = ismember(SBV(i,1),BV1(:,1));
        NDX_Max(i,m) = LOCB;
        
        [LIA,LOCB] = ismember(DBV(i,1),BV1(:,1));
        NDX_Min(i,m) = LOCB;
    end
    
    %     MBP(i) = nanmean(bp_f(iRRi(i):iRRi(i+1)));
    %     SBP(i) = nanmax(bp_f(iRRi(i):iRRi(i+1)));
    %     DBP(i) = nanmin(bp_f(iRRi(i):iRRi(i+1)));
    %     MBV(i) = nanmean(bv_f(iRRi(i):iRRi(i+1)));
    %     SBV(i) = nanmax(bv_f(iRRi(i):iRRi(i+1)));
    %     DBV(i) = nanmin(bv_f(iRRi(i):iRRi(i+1)));
    %     MBD(i) = nanmean(bd_f(iRRi(i):iRRi(i+1)));
    %     SBD(i) = nanmax(bd_f(iRRi(i):iRRi(i+1)));
    %     DBD(i) = nanmin(bd_f(iRRi(i):iRRi(i+1)));
    %        MBV(i,2) = nanmean(bv_f(iRRi(i):iRRi(i+1),2));
    %        SBV(i,2) = nanmax(bv_f(iRRi(i):iRRi(i+1),2));
    %        DBV(i,2) = nanmin(bv_f(iRRi(i):iRRi(i+1),2));
    %        tRRim(i) = nanmean(t(iRRi(i):iRRi(i+1)));
    %        RRim(i) = (t(iRRi(i+1))-t(iRRi(i)))*1000;
    
    %     %calculate critical closing pressure
    %     if max(isnan((BV1(iRRi(i):iRRi(i+1),m))))==0
    %         if max(isnan((BV1(iRRi(i):iRRi(i+1),m)))==0
    %             %           figure(10)
    %             %           plot(bp_f(iRRi(i):iRRi(i+1)),bv_f(iRRi(i):iRRi(i+1),1))
    %             %  slope = polyfit(bp_f(iRRi(i):iRRi(i+1)),bv_f(iRRi(i):iRRi(i+1),m),m);
    %
    %             %                 slope = polyfit(bp_f(iRRi(i):iRRi(i+1),m),bv_f(iRRi(i):iRRi(i+1),m),bp_f(iRRi(i):iRRi(i+1),m));
    %
    %             slope = polyfit((BV1(iRRi(i):iRRi(i+1),m)),(BV1(iRRi(i):iRRi(i+1),m)),1);
    %
    %             CCP_ICH(i,m) = -slope(2)/slope(1);
    %             %           xlim([0 max(bp_f(iRRi(i):iRRi(i+1))*1.1)])
    %             %           ylim([0 max(bv_f(iRRi(i):iRRi(i+1),1)*1.1)])
    %             %           refline(slope)
    %
    %             %           pause(0.1)
    %         end
    %
    %     end
end

MatVeliPlotMedian = nanmedian(MatVeliPlot);
DiamBAiPlotMedian = nanmedian(DiamBAiPlot);
SingleBeatFlow = MatVeliPlotMedian.*(pi()*(DiamBAiPlotMedian/2).^2).*60;

figure(98)
subplot(3,1,1)
plot(MatVeliPlotMedian(1:150))
axis tight
subplot(3,1,2)
plot(DiamBAiPlotMedian(1:150))
axis tight
subplot(3,1,3)
plot(SingleBeatFlow(1:150))
axis tight

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
% tRRi_old = tRRi;
% tRRi = tRRim;
% time = t;
%
% RRi_old = RRi;
% RRi = RRim;

figure(  'Name','Plot Velocity,andDiameterSys,Dias,andmean','NumberTitle','on','NextPlot', 'add')

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

tVDBeat(:,1) = tRRio;
tVDBeat(:,2) = tRRio;


while button == 0
    for nm=1:2
        subplot(2,1,nm)
        plot(tint,BVel1(:,nm), 'b', lw, lw_thin), hold on;
        plot(tVDBeat(:,nm),SBVel(:,nm),'ro', lw, lw_thick),
        plot(tVDBeat(:,nm),MBVel(:,nm),'ko', lw, lw_thick),
        plot(tVDBeat(:,nm),DBVel(:,nm),'go', lw, lw_thick), hold off;
        axis tight
        yl = ylim;
        yl(1) = yl(1)-yl(2)*.05;
        yl(2) = yl(2)+yl(2)*.05;
        ylim(yl);
        xlabel('Time (s)', fs, fs_label)
        
        %ylabel(titles_Select(nm) )
        %         title('Please check Only the signal that you would like to Clean ');
        %scrollplot
    end
    
    b4u = uicontrol(gcf,'Style','pushbutton','String',' Continue ','BackgroundColor','g','FontSize',15,'units','normalized','Position',[.025 .05 .1 .05],'Callback', 'button = 1 ; uiresume;');
    uiwait(gcf)
    
end
%graphfile = [infile2(1:end-4),'_DiamVel_f',num2str(gcf)];
%saveas(gcf,graphfile,'jpg');
%savefig(graphfile)
%close(gcf)

%% Drop Filter Procedure
% Preallocate vector for removed beats
tVel_bad= nan(length(tRRio),2);
MBVel_tempf_bad = nan(length(tRRio),2);
SBVel_tempf_bad = nan(length(tRRio),2);
DBVel_tempf_bad = nan(length(tRRio),2);



for m = 1:2
    
    figure(  'Name','DropFilterProcedure','NumberTitle','on','NextPlot', 'add')
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
    %ylabel(titles_Select(m) )
    %     title (labels{m})
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
            
            %             figure(1)
            clf
            %        scrsz = get(0,'ScreenSize');
            %        set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
            %    ss = get(gcf,'Position');
            
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
            xlabel('Time (s)', fs, fs_label)
            %ylabel(strcat(titles_Select(m),'_',Units_Select(m)), fs, fs_label)
            %             title (labels{m})
            %scrollplot
            xlim(XL)
            %             axis tight
            
            subplot(2,1,2)
            plot(tVDBeat(1:length(PP),m),PP, 'bo')
            hold on
            hlineU = refline(0,LPPthreshold);
            set(hlineU,'Color','r')
            hlineL = refline(0,UPPthreshold);
            set(hlineL,'Color','g')
            hold off
            %scrollplot
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
            % %    pb4l = uicontrol(gcf,'Style','pushbutton','String','Zoom In','Position',[ss(3)*.2 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*0.75; uiresume;');
            % %    pb5u = uicontrol(gcf,'Style','pushbutton','String','Shift Right','Position',[ss(3)*.8 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL+diff(XL)*.5; uiresume;');
            % %    pb5l = uicontrol(gcf,'Style','pushbutton','String','Shift Left','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL-diff(XL)*.5; uiresume;');
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
            
            %             figure(1)
            %    scrsz = get(0,'ScreenSize');
            %        scrsz = get(0,'ScreenSize');
            %        set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
            clf
            
            subplot(3,1,1:2)
            plot(tint,bv_f_temp, 'b', lw, lw_thin), hold on;
            plot(tbp,MBV_tempf,'ko', lw, lw_thick),
            plot(tbp,SBV_tempf,'ro', lw, lw_thick),
            plot(tbp,DBV_tempf,'ro', lw, lw_thick), hold off;
            %scrollplot
            xlim(XL)
            ha = get(gcf,'CurrentAxes');
            subplot(3,1,3)
            plot(tbp,PPf,'ro-')
            ylabel('HR')
            %scrollplot;
            xlim(XL2)
            ha2 = get(gcf,'CurrentAxes');
            ss = get(gcf,'Position');
            
            pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove','units','normalized','Position',[.925 .8 .05 .05],'Callback', ...
                '[x,y,butnum] = ginput(1); [mpt bad] = min(abs(tbp-x)); tbpold = tbp; MBV_tempfold = MBV_tempf; SBV_tempfold = SBV_tempf; DBV_tempfold = DBV_tempf; PPfold = PPf; tbp(bad) = []; MBV_tempf(bad) = []; SBV_tempf(bad) = []; DBV_tempf(bad) = []; PPf(bad)=[]; uiresume;');
            %                 pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','units','normalized','Position',[.925 .6 .05 .05],'Callback', ...
            %                     '[x, y] = ginput(1); x = x(end); [minpt addpt] = min(abs(tbp-x)); tbpold = tbp; MBV_tempfold = MBV_tempf; SBV_tempfold = SBV_tempf; DBV_tempfold = DBV_tempf; PPfold = PPf; tbp(end+1) = tbp(addpt); MBV_tempf(end+1) = MBV_temp(addpt); SBV_tempf(end+1) = SBV_temp(addpt); DBV_tempf(end+1) = DBV_temp(addpt); [tbp ix] = sort(tbp); MBV_tempf = MBV_tempf(ix); SBV_tempf = SBV_tempf(ix); DBV_tempf = DBV_tempf(ix);  PPf = SBV_tempf-DBV_tempf; uiresume;');
            
            pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','units','normalized','Position',[.925 .6 .05 .05],'Callback', ...
                '[x, y] = ginput(1); x = x(end); [minpt addpt] = min(abs(tRRio-x)); tbpold = tbp; MBV_tempfold = MBV_tempf; SBV_tempfold = SBV_tempf; DBV_tempfold = DBV_tempf; PPfold = PPf; tbp(end+1) = tRRio(addpt); MBV_tempf(end+1) = MBV_temp(addpt); SBV_tempf(end+1) = SBV_temp(addpt); DBV_tempf(end+1) = DBV_temp(addpt); [tbp ix] = sort(tbp); MBV_tempf = MBV_tempf(ix); SBV_tempf = SBV_tempf(ix); DBV_tempf = DBV_tempf(ix);  PPf = SBV_tempf-DBV_tempf; uiresume;');
            
            pb3u = uicontrol(gcf,'Style','pushbutton','String','Undo','units','normalized','Position',[.925 .4 .05 .05],'Callback', ...
                'tbp = tbpold; MBV_tempf = MBV_tempfold; SBV_tempf = SBV_tempfold; DBV_tempf = DBV_tempfold; PPf = PPfold; uiresume;');
            %    pb4u = uicontrol(gcf,'Style','pushbutton','String','Zoom Out','Position',[ss(3)*.4 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*1.25; uiresume;');
            %    pb4l = uicontrol(gcf,'Style','pushbutton','String','Zoom In','Position',[ss(3)*.2 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*0.75; uiresume;');
            %    pb5u = uicontrol(gcf,'Style','pushbutton','String','Shift Right','Position',[ss(3)*.8 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL+diff(XL)*.5; uiresume;');
            %    pb5l = uicontrol(gcf,'Style','pushbutton','String','Shift Left','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL-diff(XL)*.5; uiresume;');
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
        
        % Return only unique values, this is to prevent errors during interpolation
        %        tbp = unique(tbp);
        %        MBV_tempf =  unique(MBV_tempf);
        %        SBV_tempf =  unique(SBV_tempf);
        %        DBV_tempf =  unique(DBV_tempf);
        
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
        %
    end
    
    close(gcf)
end

% %% Calculate CO2 Peaks
% tco2 = time;
% SR = 1/mean(diff(tco2));
% points = ceil(1.5*SR);
% % points = ceil(SR);
% thresh = round(max(CO2)/2);
%
% %Low pass filter the CO2 signal
% button = 0;
% XL2 = [tco2(end) tco2(end)];
% XL3 = XL2;
% basefilter1 = 0;
% k=1;
%
% figure(  'Name','CalculateCO2Peaks','NumberTitle','on','NextPlot', 'add')
%
% while button==0
%
%     instruct1 = ('Use filter to smooth CO2 signal. Press LPF button to smooth signal and improve peak detection. You can press the button several times to improve smoothing. Smooth until the majority of peaks are correctly detected. Once you are happy, press Continue');
%
%     %     scrsz = get(0,'ScreenSize');
%     %     set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
%     clf
%
%     iCO2i = Minsmsa(CO2,points,thresh);
%     %     iCO2i = Minsmsa(CO2,100,10);
%
%     tCO2i = tco2(iCO2i);
%     etCO2 = [];
%     tetCO2 = [];
%     ietCO2 = [];
%
%     for i = 1:(length(iCO2i)-1)
%         [etCO2(i) temp] = max(CO2(iCO2i(i):iCO2i(i+1)));
%         %         temp = find((CO2(iCO2i(i):iCO2i(i+1))) == etCO2(i));
%         %         temp = temp+iCO2i(i);
%         ietCO2(i) = temp(1)+iCO2i(i);
%         tetCO2(i) = tco2(ietCO2(i));
%         %         if etCO2(i) < 20
%         %             etCO2(i) = NaN;
%         %         end
%     end
%     tco2 = tco2-tco2(1);
%     subplot(3,1,1:2)
%     plot(tco2,CO2);
%     hold on
%     plot(tco2(ietCO2),CO2(ietCO2),'ro')
%     hold off
%     ha2 = get(gcf,'CurrentAxes');
%     ylabel('CO2 (mmHg)')
%     scrollplot;
%     %     xlim(XL2);
%
%     subplot(3,1,3)
%     plot(tco2(ietCO2),CO2(ietCO2),'bo-');
%     hold on
%     ha2 = get(gcf,'CurrentAxes');
%     ylabel('CO2 (mmHg)')
%     hold off
%     %     xlim(XL2);
%     scrollplot;
%
%     pb3u = uicontrol(gcf,'Style','pushbutton','String','LPF','units','normalized','Position',[.925 .76 .05 .05],'Callback','basefilter1 = 1; uiresume;');
%     %    pb3l = uicontrol(gcf,'Style','pushbutton','String','LPF','units','normalized','Position',[.925 .275 .05 .05],'Callback','basefilter2 = 1; uiresume;');
%     pb7l = uicontrol(gcf,'Style','pushbutton','String','Instructions','units','normalized','Position',[.925 .3 .05 .05],'Callback','questdlg(instruct1,''Instructions'',''Continue'',''Continue''); uiresume;');
%     pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.9 .01 .1 .05],'Callback','button = 1; uiresume;');
%
%     uiwait(gcf)
%
%     if basefilter1 == 1
%         CO2old = CO2;
%         npt = round(SR/10);
%         CO2 = fastsmooth(CO2,npt*k);
%         basefilter1 = 0;
%         k=k+1;
%     end
%
%     XL2 = xlim(ha2);
%
% end
%
% %determine best threshold for detection
%
% Thresh = thresh;
% Points = points;
%
% ietCO2_old = ietCO2;
%
% button = 0;
%
% % figure(  'Name','Calculate CO2 Peaks','NumberTitle','on','NextPlot', 'add')
%
% while button == 0
%
%     instruct1 = ('This is a simple threshold cutoff. Any end tidal values below the line will be removed. Once you are happy, press Continue');
%
%
%     clf
%     %     scrsz = get(0,'ScreenSize');
%     %     set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
%     %    ss = get(gcf,'Position');
%
%     good = find(CO2(ietCO2)>thresh);
%     ietCO2_plot = ietCO2(good);
%
%     clf
%     plot(tco2,CO2, 'b', lw, lw_thin), hold on;
%     plot(tco2(ietCO2_plot),CO2(ietCO2_plot),'ro')
%     refline(0,thresh);
%     axis tight
%     yl = ylim;
%     yl(1) = yl(1)-yl(2)*.05;
%     yl(2) = yl(2)+yl(2)*.05;
%     ylim(yl);
%     xlabel('Time (s)', fs, fs_label)
%     ylabel('CO2  (mmHg)', fs, fs_label)
%     scrollplot
%
%     pb1u = uicontrol(gcf,'Style','pushbutton','String','Thresh+','units','normalized','Position',[.925 .9 .05 .05],'Callback', ...
%         'thresh=thresh+Thresh*0.05; uiresume;');
%     pb2u = uicontrol(gcf,'Style','pushbutton','String','Thresh-','units','normalized','Position',[.925 .8 .05 .05],'Callback', ...
%         'thresh=thresh-Thresh*0.05; uiresume;');
%     pb3u = uicontrol(gcf,'Style','pushbutton','String','Thresh++','units','normalized','Position',[.925 .7 .05 .05],'Callback', ...
%         'thresh=thresh+Thresh*0.2; uiresume;');
%     pb4u = uicontrol(gcf,'Style','pushbutton','String','Thresh--','units','normalized','Position',[.925 .6 .05 .05],'Callback', ...
%         'thresh=thresh-Thresh*0.2; uiresume;');
%     %     pb5u = uicontrol(gcf,'Style','pushbutton','String','Points+','units','normalized','Position',[.925 .5 .05 .05],'Callback', ...
%     %         'points=points+Points*0.1; uiresume;');
%     %     pb6u = uicontrol(gcf,'Style','pushbutton','String','Points-','units','normalized','Position',[.925 .4 .05 .05],'Callback', ...
%     %         'points=points-Points*0.1; uiresume;');
%     %     pb7u = uicontrol(gcf,'Style','pushbutton','String','Points++','units','normalized','Position',[.925 .3 .05 .05],'Callback', ...
%     %         'points=points+Points*0.5; uiresume;');
%     %     pb8u = uicontrol(gcf,'Style','pushbutton','String','Points--','units','normalized','Position',[.925 .2 .05 .05],'Callback', ...
%     %         'points=points-Points*0.5; uiresume;');
%     % %    pb4l = uicontrol(gcf,'Style','pushbutton','String','Zoom In','Position',[ss(3)*.2 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*0.75; uiresume;');
%     % %    pb5u = uicontrol(gcf,'Style','pushbutton','String','Shift Right','Position',[ss(3)*.8 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL+diff(XL)*.5; uiresume;');
%     % %    pb5l = uicontrol(gcf,'Style','pushbutton','String','Shift Left','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL-diff(XL)*.5; uiresume;');
%     pb7l = uicontrol(gcf,'Style','pushbutton','String','Instructions','units','normalized','Position',[.925 .125 .05 .05],'Callback','questdlg(instruct1,''Instructions'',''Continue'',''Continue''); uiresume;');
%     pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
%
%     uiwait(gcf)
%
%
% end
%
% ietCO2 = ietCO2(good);
% %
% %
% button = 0;
%
% % figure(  'Name','Calculate CO2 Peaks','NumberTitle','on','NextPlot', 'add')
%
% while button == 0
%
%     instruct1 = ['Add or remove peak detections points till CO2 signal is correctly peak detected. Once you are happy, press Continue'];
%
%
%     %     scrsz = get(0,'ScreenSize');
%     %     set(gcf,'Position',[1 scrsz(4)*.05 scrsz(3) scrsz(4)*.85])
%     clf
%     subplot(3,1,1:2)
%     plot(tco2,CO2, 'b', lw, lw_thin), hold on;
%     plot(tco2(ietCO2),CO2(ietCO2),'ro'), hold off;
%     scrollplot
%     %     axis tight
%     %     yl = ylim;
%     %     xlim(XL)
%     ha = get(gcf,'CurrentAxes');
%     subplot(3,1,3)
%     plot(tco2(ietCO2),CO2(ietCO2))
%     ylabel('ETCO2')
%     scrollplot;
%     ha2 = get(gcf,'CurrentAxes');
%
%     ss = get(gcf,'Position');
%
%     pb1u = uicontrol(gcf,'Style','pushbutton','String','Remove','units','normalized','Position',[.925 .7 .05 .05],'Callback', ...
%         '[x, y] = ginput(1); x = x(end); [mpt bad] = min(abs(tco2(ietCO2)-x)); ietCO2old = ietCO2; ietCO2(bad) = []; uiresume;');
%     pb2u = uicontrol(gcf,'Style','pushbutton','String','Add','units','normalized','Position',[.925 .5 .05 .05],'Callback', ...
%         '[x, y] = ginput(1); x = round(x(end).*SR); [mpt good] = max(CO2(x-round(SR):x+round(SR))); ietCO2old = ietCO2; ietCO2(end+1) = (good + x-round(SR)); ietCO2 = sort(ietCO2); ietCO2 = unique(ietCO2); uiresume;');
%     pb3u = uicontrol(gcf,'Style','pushbutton','String','Undo','units','normalized','Position',[.925 .3 .05 .05],'Callback', ...
%         'ietCO2 = ietCO2old; uiresume;');
%     %    pb4u = uicontrol(gcf,'Style','pushbutton','String','Zoom Out','Position',[ss(3)*.4 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*1.25; uiresume;');
%     %    pb4l = uicontrol(gcf,'Style','pushbutton','String','Zoom In','Position',[ss(3)*.2 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL*0.75; uiresume;');
%     %    pb5u = uicontrol(gcf,'Style','pushbutton','String','Shift Right','Position',[ss(3)*.8 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL+diff(XL)*.5; uiresume;');
%     %    pb5l = uicontrol(gcf,'Style','pushbutton','String','Shift Left','Position',[ss(3)*.6 1 ss(3)*.1 ss(4)*.05],'Callback','XL = xlim; XL = XL-diff(XL)*.5; uiresume;');
%     pb7l = uicontrol(gcf,'Style','pushbutton','String','Instructions','units','normalized','Position',[.925 .2 .05 .05],'Callback','questdlg(instruct1,''Instructions'',''Continue'',''Continue''); uiresume;');
%     pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .1 .05 .05],'Callback','button = 1; uiresume;');
%
%     uiwait(gcf)
%
%     ha = get(gcf,'CurrentAxes');
%     XL = xlim(ha2);
%
% end
% savefig('Calculate CO2 Peaks')
% close(gcf)
%
% tetCO2 = tco2(ietCO2);
% etCO2 = CO2(ietCO2);
% good = find(~isnan(etCO2));
% ietCO2 = ietCO2(good);
% etCO2 = etCO2(good);
% tetCO2 = tetCO2(good);
%

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

%Carotid Stiffness Calculations
%Butlin M (2007). "Structural and function effects on large artery stiffness: an in-vivo experimental investigation". Graduate School of Biomedical Engineering, University of New South Wales.
% PWV = sqrt((E_inc * h)/ (2*r*rho))
% Ratio of ln(systolic/diastolic pressures) to (relative change in diameter) Ultrasound*
% b = ln(Ps/Pd) / ((Ds - Dd )/ Dd)

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
%
% startime = startime - t(1);
% endtime = endtime - t(1);

% ikp = find(iRRiV<length(t));
% [s d m irri]= SDMcalc (CBF,t,'flow',1,iRRiV(ikp));c


% t = tRRi;

% plot(tRRiV),iRRiV(1:end-1)
button = 0;

figure(  'Name','PlotofDiameter,Velocity,andFlow','NumberTitle','on','NextPlot', 'add')

% clf('reset')
switch Velocitydlg
    case   'MATLAB'
        
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
        
        button = 0;
        
        figure(  'Name','PlotofMeanDiameter,Velocity,Flow','NumberTitle','on','NextPlot', 'add')
        
        % clf('reset')
        
        
        
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
        
    case  'Brachial_Analyzer'
        while button == 0
            hold on;
            % Plot of BP
            subplot(3,2,1)
            plot(tint,BpFlow,'k');hold on;
            %             plot(tRRi,SBV(:,1),'bo');
            plot(tRRi,MBV(:,1),'ko');
            plot(tbp_bad(:,1),MBVel_tempf_bad(:,1),'ro');
            
            %             plot(tRRi,DBV(:,1),'bo');
            hold off
            axis tight
            title(titles_cell(1))
            ylabel('mmHg')
            %scrollplot
            
            % Plot of Diameter
            subplot(3,2,2)
            plot(t,DiamFlow_cm,'k');hold on;
            %             plot(tRRi,SBV(:,3),'bo');
            plot(tRRi,MBV(:,3),'ko');
            plot(tbp_bad(:,3),MBVel_tempf_bad(:,3),'ro');
            
            %             plot(tRRi,DBV(:,3),'bo');
            hold off
            axis tight
            ylim([min(DBV(:,3))*.9 max(SBV(:,3))*1.1])
            title(titles_cell(7))
            ylabel('cm')
            %scrollplot
            
            % Plot of Cerebral Velocity
            subplot(3,2,3)
            plot(time,CBV, lw, lw_thin),hold on
            
            %             plot(tRRi,SBV(:,2),'ro');
            plot(tRRi,MBV(:,2),'ko');
            plot(tbp_bad(:,2),MBVel_tempf_bad(:,2),'ro');
            
            %             plot(tRRi,DBV(:,2),'ro');
            hold off
            axis tight
            %     ylim([min(DBV)/10*.9 max(SBV)/10*1.1])
            title(titles_cell(6))
            ylabel('cm/sec')
            %scrollplot
            
            % Plot of Cerebral Blood Flow
            subplot(3,2,4)
            %             plot(time,CBF,'k')
            
            %             plot(tRRi,SBD_flow,'ro');hold on
            plot(tRRi,MBD_flow,'ko');
            %             plot(tRRi,DBD_flow,'ro');
            hold off
            axis tight
            title('Cerebral Blood Flow')
            ylabel('Flow mL/min')
            %scrollplot
            
            % Plot of Blood Flow per Hundred Grams of Brain Tissue
            subplot(3,2,5)
            
            %             plot(time,CBF_Hgbt)
            %             plot(tRRi,SBD_flow_Hgbt,'ro');hold on
            plot(tRRi,MBD_flow_Hgbt,'ko');
            %             plot(tRRi,DBD_flow_Hgbt,'ro');
            hold off
            title('Blood Flow per Hundred Grams of Brain Tissue ')
            ylabel('Flow mL/g.min')
            %scrollplot
            hold off
            
            pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
            uiwait(gcf)
            
        end
        
        %         graphfile2 = [infile2(1:end-4),'_FlowSignals_f',num2str(gcf)];
        %         saveas(gcf,graphfile2,'jpg');
        %         savefig(graphfile2)
        %         close(gcf)
        
        
        
        button = 0;
        
        figure(  'Name','PlotofMeanDiameter,Velocity,Flow','NumberTitle','on','NextPlot', 'add')
        
        % clf('reset')
        
        while button == 0
            
            subplot(411)
            plot(tRRi,MBV(:,3),'-b*');
            title('Mean Diameter')
            ylabel('Cm')
            %scrollplot
            
            
            subplot(412)
            plot(tRRi,MBV(:,2),'-k*');
            title('Mean Blood Velocity')
            ylabel('Cm/s')
            %scrollplot
            
            subplot(413)
            plot(tRRi,MBD_flow,'-r*')
            title('Mean Blood flow')
            ylabel('Flow mL/min ')
            %scrollplot
            
            %We cannont determine systolics and diastolics accurately from Brachial analyzer, therefore distensibility and stiffness calulcations is not reported
            %             subplot(324)
            %             plot(tRRi,DD,'-g*')
            %             title('Distensibility')
            %             ylabel('')
            %             scrollplot
            
            %             subplot(325)
            %             plot(tRRi,Beta,'-b*')
            %             hold on
            %             title('Stiffness')
            %             ylabel('')
            %             scrollplot
            
            subplot(414)
            plot(tRRi,CVR,'-m*')
            title('Resistance')
            ylabel('')
            %scrollplot
            
            
            pbexit = uicontrol(gcf,'Style','pushbutton','String','Continue','units','normalized','Position',[.925 .05 .05 .05],'Callback','button = 1; uiresume;');
            uiwait(gcf)
        end
end
%graphfile = [infile2(1:end-4),'_MeanDiamVelFlow_f',num2str(gcf)];
%saveas(gcf,graphfile,'jpg');
%savefig(graphfile)
%close(gcf)

%% Save data
Question = 'Output File Already Exists, Overwrite?';
Qtitle= 'Data Aligned Analysis';


% checking if file originally exist
check_exist = exist(folder_name, 'file');
if check_exist >0
    ButtonName = questdlg(Question,Qtitle);
    switch ButtonName
        case 'Yes'
            
            folder = mkdir (name,folder_name);
            new_dir=[name '\' folder_name];
            cd(new_dir)
            outfile_clean = [folder_name '_AlignedClean.mat'];
            save([new_dir '\' outfile_clean])
            
            
        case 'No'
            new_name = inputdlg('Rename File', 'Qtitle', 1,{folder_name});
            new_name=cell2str (new_name);
            new_name(1)=[];
            new_name(end)=[];
            folder = mkdir (name , new_name);
            new_dir=[name '\' new_name];
            cd(new_dir)
            outfile_clean = [folder_name '_AlignedClean.mat'];
            save([new_dir '\' outfile_clean])
            
            
        case 'Cancel'
            STATUS = 0;
            return;
    end
    
    
    
elseif check_exist ==0
    folder = mkdir (name,folder_name);
    new_dir=[name '\' folder_name];
    cd(new_dir)
    % Saving Matlab File
    Matlab_Saving_bar_massage = ('Saving Data AS ****AlignedClean.mat');
    Matlab_Saving_bar = waitbar(.5,Matlab_Saving_bar_massage);
    hw=findobj(Matlab_Saving_bar,'Type','Patch');
    set(hw,'FaceColor','b','EdgeColor','b')
    outfile_clean = [folder_name 'AlignedClean.mat'];
    save([pathname folder_name '\' outfile_clean])
    waitbar(1,Matlab_Saving_bar,Matlab_Saving_bar_massage);
    close(Matlab_Saving_bar)
end


%% Saving Outputs to Excel

% %This section replaces the calculated systolics and diastolics with nans
% %for cases in which Brachial analyzer velocity is selected
% if ~strcmp(Velocitydlg,'Brachial_Analyzer')
%     SBV_Vel = (SBV(:,2));
%     DBV_Vel = (DBV(:,2));
% elseif strcmp(Velocitydlg,'Brachial_Analyzer')
%     SBV_Vel = nan(size(SBV(:,2)));
%     DBV_Vel = nan(size(DBV(:,2)));
%
%     SBD_flow = nan(size(SBD_flow));
%     DBD_flow = nan(size(DBD_flow));
%
%     SBD_flow_Hgbt = nan(size(SBD_flow_Hgbt));
%     DBD_flow_Hgbt = nan(size(DBD_flow_Hgbt));
%
%     Sys_Vel_Mean = nan(size(Sys_Vel_Mean));
%     Dias_Vel_Mean = nan(size(Dias_Vel_Mean));
%
%     Sys_Flow_Mean = nan(size(Sys_CBF_Mean));
%     Dias_Flow_Mean = nan(size(Dias_CBF_Mean));
%
%     SBD_flow_Hgbt = nan(size(SBD_flow_Hgbt));
%     DBD_flow_Hgbt = nan(size(DBD_flow_Hgbt));
%
%     DD = nan(size(DD));
%     Beta = nan(size(Beta));
% end
% % save([new_dir '\' name])
% h = waitbar(0,'Please wait...Saving Beat-by-Beat Data as **.xls');
% subject_ID = {infile(1:6)};
% book_name = cell2mat(strcat(subject_ID, Procedure, '_Brain_Blood_Flow'));
% % A = {'Time','BDIAMM(cm)','CBV','CBF','CBF_HGBT', 'Distensibility','Stiffness','CVR','tRRim'};
% % A = {'Time','Sys_Diam(cm)','Mean_Diam(cm)','Dias_Diam(cm)','Sys_Vel(cm/s)','Mean_Vel(cm/s)','Dias_Vel(cm/s)','Sys_CBF(mL/min)','Mean_CBF(mL/min)','Dias_CBF(mL/min))','CBF_HGBT(Flow mL/g.min)', 'Distensibility','Stiffness','CVR'};
% A = {'Time','Sys_BP(mmHg)','Mean_BP(mmHg)','Dias_BP(mmHg)','Sys_Diam(cm)','Mean_Diam(cm)','Dias_Diam(cm)','Sys_Vel(cm/s)','Mean_Vel(cm/s)','Dias_Vel(cm/s)','Sys_CBF(mL/min)','Mean_CBF(mL/min)','Dias_CBF(mL/min)','Sys_CBF_Hgbt(mL/g.min)','Mean_CBF_Hgbt(mL/g.min)','Dias_CBF_Hgbt(mL/g.min)','Distensibility','Stiffness','CVR'};
%
% A_mean = {'subject_ID','Artery(Common=1,Internal=2)','Position(Seated=1,Supine = 2)','Sys_BP_Mean','Mean_BP_Mean','Dias_BP_Mean','Sys_Vel_Mean','Mean_Vel_Mean','Dias_Vel_Mean','Sys_Diam_Mean','Mean_Diam_Mean','Dias_Diam_Mean','Sys_Flow_Mean','Mean_Flow_Mean','Dias_Flow_Mean','Sys_CBF_Hgbt_Mean','Mean_CBF_Hgbt_Mean','Dias_CBF_Hgbt_Mean',' Distensibility_Mean',' Stiffness_mean','CVR_mean'};
% letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];
%
% %Write Excel Sheet Cerebral bood flow values
%
% time = t;
% output = {tRRi,SBV(:,1),MBV(:,1),DBV(:,1),SBV(:,3),MBV(:,3),DBV(:,3),SBV_Vel,MBV(:,2),DBV_Vel,SBD_flow,MBD_flow,DBD_flow, SBD_flow_Hgbt,MBD_flow_Hgbt,DBD_flow_Hgbt,DD,Beta,CVR};
%
%
% B = length (output);
% sheet = 1;
%
% xlswrite(book_name, A,sheet,'A1');
%
% for c = 1:B
%     output2 =  reshape(cell2mat(output(1,c)),[],1);
%     range=[letters(c) num2str(2)];
%     %     range2=[letters(c) num2str(2)];
%     %     xlswrite(book_name, output2,sheet,'CerebralBloodFlow',range) %Need to
%     %     figure out how to rename sheet without creating a new sheet
%     xlswrite(book_name, output2,sheet,range)
%
%     waitbar(c / B)
% end
%
% %Write Excel Sheet Average values
% sheet2 = 2; %Define sheet to use
%
% output_Mean = {Artery,Position,Sys_BP_Mean,Mean_BP_Mean,Dias_BP_Mean,Sys_Vel_Mean,Mean_Vel_Mean,Dias_Vel_Mean,Sys_Diam_Mean,Mean_Diam_Mean,Dias_Diam_Mean,Sys_CBF_Mean,Mean_CBF_Mean,Dias_CBF_Mean,Sys_CBF_Hgbt_Mean,Mean_CBF_Hgbt_Mean,Dias_CBF_Hgbt_Mean,mean_Distensibility, Stiffness_mean, CVR_mean};
% xlswrite(book_name,A_mean,sheet2,'A1');%Write headers starting in cell A1
% xlswrite(book_name,subject_ID,sheet2,'A2');%Write Subject ID in A2
% xlswrite(book_name,output_Mean,sheet2,'B2');%Write Mean values starting in cell B2
%
% close(h)
%
% % Compile each run into a single Excel Sheet located \\vhaeasfpc4\Projects\WRIISC\STUDY - Head Injury\DATA_Project IDs
% book_name2 = ['HeadInjury_Report'];
%
% %Define directory
% dir1 = '\\vhaeasfpc4\Projects\WRIISC\STUDY - Head Injury\DATA_Project IDs';
% if  exist(dir1,'dir')%Save to specified directory if it exists
%     % Define all Letters to use
%     letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];
%
%     % Write the header file to excel Report_OutputFile.xlsx
%     xlswrite([dir1 '\' book_name2], A_mean,sheet,'A1');
%
%     % Request the numeric data, text, and a copy of the unprocessed (raw)
%     % data from the first worksheet:
%     [ndata, text, alldata] = xlsread([dir1 '\' book_name2],'Sheet1'); % read the existing data for appending this vHIT
%     if isempty(ndata)
%         xlswrite([dir1 '\' book_name2],subject_ID ,sheet,'A2');
%         xlswrite([dir1 '\' book_name2], output_Mean,sheet,'B2');
%
%     elseif ~isempty(ndata)
%         nData = (size(alldata));
%         % last row with data in the file
%         nRows = nData(1);
%         % last row with data in the file
%         nColmn = nData(2);
%
%         % plus 1 to write subject name and data in the next line
%         nameRange=['A' num2str(nRows + 1)];
%         dataRange=['B' num2str(nRows + 1)];
%         xlswrite([dir1 '\' book_name2],subject_ID ,sheet,nameRange);
%         xlswrite([dir1 '\' book_name2], output_Mean,sheet,dataRange);
%     end
%
% end
%
%
%
% %% Save figures
% % Convert .fig file to jpg file.
%
%
% % button = questdlg('Do you want to save the graphs?','Save Graphs','Yes','No','Help','No');
% % if strcmp(button,'Yes')
% % filterspec = '*.fig';
% % Title = 'Pick Figures you would like to save';
% % [infile,pathname] = uigetfile(filterspec,Title);
% %
% % if infile == 0
% %     %        errordlg('File not found','File Error');
% %     return;
% % end
% Fig_saving_bar_massage= ('Saving Figure as **.jpg');
% Fig_saving_bar = waitbar(0,Fig_saving_bar_massage);
% hw=findobj(Fig_saving_bar,'Type','Patch');
% set(hw,'FaceColor','g','EdgeColor','g')
% d = dir(fullfile(new_dir,'*.fig'));
% %     d=dir(dir(name,'*.fig')'*.fig');
% for k=1:length(d)
%     fname=d(k).name;
%     %     [pathstr, name, ext] = fileparts(fname);
%     %     movefile(fname, fullfile(pathstr, [name '.txt']))
%
%     % ...
%     %     if fname==0
%     %         %break;
%     %     end
%     %     if isequal(fname,0) || isequal(pathstr,0)
%     %         disp('User pressed cancel')
%     %     else
%     %         disp(['User selected ', fullfile(pathstr, fname)])
%     %     end
%     openfig(fname);
%     %             figure(k)
%     %             set(gcf,'name',cell2str(initials))
%     %             set(0,'CurrentFigure')
%     %             orient landscape
%     %        set(gcf,'PaperUnits','inches', 'PaperType','Letter' ,'PaperOrientation','landscape')
%     %        print
%
%     %             graphfile = [pathname,folder_name,'\',infile(1:end-5),'f',num2str(i)];
%     graphfile = [pathname,folder_name,'\',fname(1:end-4)];
%     % graphfile = [pathname,fname(1:end-4)];
%     saveas(gcf,graphfile,'jpg');
%     %             close(gcf)
%     waitbar(.5*k,Fig_saving_bar,Fig_saving_bar_massage);
%
% end
% close all;
% close(Fig_saving_bar);
% % infile1 = [pathname infile];
%
%
% %     %
% %     if  strcmp(Velocitydlg,'MATLAB')
% %
% %         for  i = [4,7,10,12,13,14,16,17];
% %             waitbar(.5*i,Fig_saving_bar,Fig_saving_bar_massage);
% %             figure(i)
% %             set(gcf,'name',cell2str(initials))
% %             set(0,'CurrentFigure')
% %             orient landscape
% %             %        set(gcf,'PaperUnits','inches', 'PaperType','Letter' ,'PaperOrientation','landscape')
% %             %        print
% %
% %             graphfile = [pathname,folder_name,'\',infile(1:end-5),'f',num2str(i)];
% %             saveas(gcf,graphfile,'jpg')
% %         end
% %     elseif   strcmp(Velocitydlg,'Brachial_Analyzer')
% %         for  i = [4,7,10,12,13,14] ;
% %             waitbar(.5*i,Fig_saving_bar,Fig_saving_bar_massage);
% %             figure(i)
% %             set(gcf,'name',cell2str(initials))
% %             set(0,'CurrentFigure')
% %             orient landscape
% %             %        set(gcf,'PaperUnits','inches', 'PaperType','Letter' ,'PaperOrientation','landscape')
% %             %        print
% %
% %             graphfile = [pathname,folder_name,'\',infile(1:end-5),'f',num2str(i)];
% %             saveas(gcf,graphfile,'jpg')
% %         end
% %
% %     end
% %     close(Fig_saving_bar)
%
% % end
% %% Saving Text Files
% Notes_saving_bar_massage= ('Saving Figure as **.jpg');
% Notes_saving_bar = waitbar(0,Notes_saving_bar_massage);
% hw=findobj(Notes_saving_bar,'Type','Patch');
% set(hw,'FaceColor','g','EdgeColor','g')
%
% fid=fopen(folder_name,'wt');
% if ~isempty(Notes_list)
%     for i=1:length(Notes_list);
%         waitbar((1/length(Notes_list))*i,Notes_saving_bar,Notes_saving_bar_massage);
%         fprintf(fid,'%s',mat2str(cell2mat(Notes_list{i})));
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
% end
% close(Notes_saving_bar)
%
% %% Write table to file with imformation pertaining to analysis
%
% if mauid_exist == 1;
%     diamsoft='Brachial Analyzer';
% else
%     diamsoft='MAUI';
% end
%
% % Create a table
% switch diamsoft
%     case 'Brachial Analyzer'
%         LastName = {'Report_Date';'ScriptName';'Operator_ID';'Subject_ID';'Pat_ID_type';'Study_ID';'Analysis_Date'; 'Filename'};
%         Report_Date = T{2};
%         ScriptName = mfilename;
%         Operator_ID = cell2mat(initials);
%         if isempty (Operator_ID)
%             Operator_ID = nan;
%         end
%         Subject_ID = T{3};
%         Pat_ID_type = T{4};
%         Study_ID = T{5};
%         % Analysis_Date = datetime('now'); Works in Matlab 2015
%         Analysis_Date =  datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')
%         DateString = datestr(Analysis_Date);
%         Filename = graphfile;
%         % Filename = pathname;
%         Analysis_Report = [{Report_Date};{ScriptName};{Operator_ID};{Subject_ID};{Pat_ID_type};{Study_ID};{DateString};{Filename}];
%
%
%         T = table(Analysis_Report,'RowNames',LastName)
%
%         % Write the table, T, to a comma delimited text file, called myPatientData.dat, and display the file contents.
%         AnalysisReport = [num2str(Subject_ID) '_AnalysisFileSummary_' cell2mat(initials) '.txt']
%
%         writetable(T,AnalysisReport,'WriteRowNames',true)
%         clc;
%         % type(AnalysisReport)
%         T
%         close all;
%
%     case 'MAUI'
%         LastName = {'ScriptName';'Analysis_Date'};
%         % Report_Date = T{1,2};
%         ScriptName = mfilename;
%         Operator_ID = cell2mat(initials);
%         if isempty (Operator_ID)
%             Operator_ID = nan;
%         end
%         % Subject_ID = T{2,2};
%         % Pat_ID_type = T{3,2};
%         % Study_ID = T{4,2};
%         % Analysis_Date = datetime('now'); Works in Matlab 2015
%         Analysis_Date =  datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM');
%         DateString = datestr(Analysis_Date);
%         Filename = graphfile;
%         % Filename = pathname;
%         % Analysis_Report = [Report_Date;ScriptName;Operator_ID;Subject_ID;Pat_ID_type;Study_ID;DateString; Filename];
%         Analysis_Report = {ScriptName;DateString};
%         % Write the table, T, to a comma delimited text file, called myPatientData.dat, and display the file contents.
%         Subject_ID = T{3};
%         AnalysisReport = [num2str(Subject_ID) '_AnalysisFileSummary_' cell2mat(initials) '.txt']
%
%         T = table(Analysis_Report,'RowNames',LastName);
%
%         % Write the table, T, to a comma delimited text file, called myPatientData.dat, and display the file contents.
%         [cell2mat(subject_ID) '_AnalysisFileSummary_' cell2mat(initials) '.txt']
%
%         writetable(T,AnalysisReport,'WriteRowNames',true)
%         clc;
%         % type(AnalysisReport)
%         T;
%         close all;
% end
