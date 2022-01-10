clear;
close all;

%% define basic parameters
sound_speed = 1540; % m/s speed in the air should be the speed in blood
                                                                                                %% Aaslid says that the sound_speed = 1540m/s
%% parameters from video file
%Hz sensor sound frequency f_source = 4.3 * 10^6;

%%Create UI to select sensor sound frequency buttons and input angle
f_source = 0;
answer = questdlg('What is the sensor sound frequency? f_source = ',...
                    'Sensor Frequency', '4 Mhz', '4.3 Mhz', '4.6 Mhz', '4.3 Mhz');            
switch answer
    case '4 Mhz'
        f_source = 4 * 10^6;
    case '4.3 Mhz'
        f_source = 4.3 * 10^6;
    case '4.6 Mhz'
        f_source = 4.6 * 10^6;
end

angle_input = inputdlg('Enter the angle in degrees, angle = ','Angle Input');
angle = str2double(cell2mat(angle_input));

sound_input = inputdlg('Enter the speed of sound in blood, m/s = ','Sound Speed');
sound_speed = str2double(cell2mat(angle_input));

%% load audio file from .mp4, .m4a, .mp3, or .wav
[file, path] = uigetfile('*.mp4;*.m4a;*.mp3;*.wav');
file_path = [path, file];
[y,Fs] = audioread(file_path);
sound_data = y(:,1);

%% low pass filter orignal sound file
sound_data = lowpass(sound_data,5,Fs);

%% peform short time fft to extract frequence change over time
% STFT plot
figure()
stft(sound_data,Fs,'Window',kaiser(256,5),'OverlapLength',220,'FFTLength',512);
%win = hamming(256,'periodic');
[s,f,t] = stft(sound_data,Fs,'Window',kaiser(256 * 4, 5),'OverlapLength',220 * 4,'FFTLength',512 * 4);


% build frequence change data
 f_change_overtime = [];
for i = 1 : size(s, 2)
    % max magniture at given time
    mag_at_time = 20*log10(abs(s(:,i)));
    % find max index
    [M, index] = max(mag_at_time);
    % find corresponding f
    f_at_time = f(index);
    % build f_change_overtime
    f_change_overtime = [f_change_overtime, abs(f_at_time)];
end
% frequency change over time
figure()
hold on;
plot(f_change_overtime);
% use hilbert transform
f_change_overtime_envelope = envelope(f_change_overtime, 10, 'peak');
plot(f_change_overtime_envelope);
title('Frequence Over Time');

%% intrapolate on frequence change over time to get speed over time
% original length
x_original = 1 : length(sound_data)/length(f_change_overtime) : length(sound_data);
% length after interpolation
xq = 1 : length(sound_data);
f_change_overtime_interpolated = interp1(x_original, f_change_overtime_envelope, xq, 'pchip');
figure()
plot(f_change_overtime_interpolated);
title('Frequence Over Time Interpolated');

%% calcalte speed based on frequence
% formula v_blood = (f_o / f_sound - 1) * v_sound
speed_overtime = (f_change_overtime_envelope / f_source - 1) * sound_speed;
figure()
plot(speed_overtime);
title('Speed Over Time Interpolated');

%% compare with previous results
previous_res = load('E:\Min and Steph\subjects_data\matlab codes\sample speed data\CON178_Common_Seated_Velocity_SI.mat');
figure()
plot(previous_res.V);
title('Previous Speed Over Time Interpolated');

%% plot together
figure()
subplot(2,1,1)
plot(speed_overtime);
title('Min')
subplot(2,1,2)
plot(previous_res.V);
title('Previous')
% same scale
figure()
hold on
normalizedData = mat2gray(speed_overtime);
plot(normalizedData);
plot(previous_res.V / max(previous_res.V));
title('Min vs Previous on same scale')


%% if we interpolate both

% new speed
% intrapolate on frequence change over time to get speed over time
x_original = 1 : length(sound_data)/length(f_change_overtime) : length(sound_data);
% length after interpolation
xq = 1 : length(sound_data);
f_change_overtime_interpolated = interp1(x_original, f_change_overtime_envelope, xq, 'pchip');

% old speed
x_original = 1 : length(sound_data)/length(previous_res.V) : length(sound_data);
% length after interpolation
xq = 1 : length(sound_data);
previous_V_interpolated = interp1(x_original, previous_res.V, xq, 'pchip');

figure()
hold on
normalizedData = mat2gray(f_change_overtime_interpolated);
plot(normalizedData);
plot(previous_V_interpolated / max(previous_V_interpolated));
title('Min vs Previous on same scale intrapolated')
legend('Min', 'previous');
