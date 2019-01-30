% Example of how to use the 

clearvars    %MATLAB2016
close all
clc

% Spectrogram parameters
fft_size = 512;
overlap = 80; % overlapping

% % load the observation
% [x,fs,tx] = data_testbase_loading(name);
[x,fs] = audioread('MPBW.wav');
x= x-mean(x) ; x = x/max(abs(x));
% Filtering of frequencies below 5 Hz
[b,a]=butter(10,5/(fs/2),'high'); x = filter(b,a,x);

% Spectrogram
[stft,f,t,p] = spectrogram(x,hann(fft_size),round((overlap/100)*fft_size),fft_size,fs);

% Detection parameters
delta_t = 3;% 3(s)
delta_f = 0.8;%(Hz)
signal_mini_duration_s = 5; % 2(s)

dT_temp = 1/fs;
dT_spec = t(2)-t(1);
signal_mini_duration_samp_temp = round(signal_mini_duration_s/dT_temp); 
signal_mini_duration_samp_spec = round( signal_mini_duration_s/dT_spec );


% Run the ridges detector
addpath Ridge_detector
addpath Ridge_detector/util
addpath Ridge_detector/detector
[my_tracks] = lea_run_shyam_detector(p,f,t,signal_mini_duration_samp_spec,delta_t/(t(2)-t(1)),8);

figure
% subplot of the spectrogram and the ridge detector results
ax1 = subplot(2,1,1);
imagesc(t,f,10*log10(p))
set(ax1,'clim',[-85 -20])
axis xy
hold on
for i = 1:length(my_tracks)
    time_tf =  my_tracks{i}(1,:);
    freq_tf = my_tracks{i}(2,:);
    plot(time_tf,freq_tf,'k')
end
ylim([10 40])
ylabel('Frequency (Hz)')

%% Reconstruction of the sounds
% If stft is already computed
[x_reconstructed, tx] = sound_reconstruction(stft, f, t, fs, fft_size, overlap, my_tracks,[0 15]);

% If stft is not already computed
% [x_reconstructed, tx] = sound_reconstruction_stft(x, fs, fft_size, overlap, my_tracks, [21 24]);

% subplot of the input waveform and the reconstructed signal
subplot(2,1,2);
hold on
plot((0:length(x)-1)/fs,x)
plot(tx,x_reconstructed)
xlim([0 max(tx)])
grid on 
xlabel('Time (s)')
ylabel('Amplitude')

