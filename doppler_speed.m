%clear;
close all;

%% define basic parameters
sound_speed = 1500; % m/s speed in the air should be the speed in blood
f_source = 2000000; %Hz sensor sound frequency, don't know

%% load audio file from .mp4 or .avi file
[file, path] = uigetfile('*');
file_path = [path, file];
[y,Fs] = audioread(file_path);

%% peform short time fft to extract frequence change over time
sound_data = y(:,1);
% STFT plot
figure()
stft(sound_data,Fs,'Window',kaiser(256,5),'OverlapLength',220,'FFTLength',512);
win = hamming(256,'periodic');
[s,f,t] = stft(sound_data,Fs,'Window',kaiser(256,5),'OverlapLength',220,'FFTLength',512);
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
f_change_overtime_envelope = envelope(f_change_overtime, 60, 'peak');
plot(f_change_overtime_envelope);

%% intrapolate on frequence change over time to get speed over time
% original length
x_original = 1 : length(sound_data)/length(f_change_overtime) : length(sound_data);
% length after interpolation
xq = 1 : length(sound_data);
f_change_overtime_interpolated = interp1(x_original, f_change_overtime_envelope, xq, 'pchip');
figure()
plot(f_change_overtime_interpolated);

%% calcalte speed based on frequence
% formula (f_o / f_sound - 1) * v_sound
speed_overtime = (f_change_overtime_interpolated / f_source - 1) * sound_speed;
figure()
plot(speed_overtime);
figure(5)
hold on
plot(speed_overtime);
plot(V);



