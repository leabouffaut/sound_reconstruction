% [x_reconstructed, tx] = sound_reconstruction_stft(x, fs, window, overlap, my_tracks,f_limit)
%
% Author: Lea Bouffaut, PhD candidate at the Naval Academy Research
% Institute, Brest, France.
% Date: 01-30-2019
%
% This function reconstructs waveforms from tonal sounds detected by the 
% ridge detector which is based on an image-processing technique applied to
% spectrograms of acoustic signals.
%
% THE STFT IS COMPUTED IN THIS FUNCTION, IF IT IS ALREADY
% COMPUTED ELSEWHERE, CHECK THE FUNCTOION "sound_reconstruction.m".
%
% The acoustic signal short time fourier tramsform (stft) is derived to 
% provide a time-frequency representation conservative of the complex
% signal, hence the phase.
% 
% A binary detection matrix is generated from the detected time-frequency 
% coordinates provided by the ridge detector. It is applied to the stft of 
% the acoustic input signal as a binary filter.
%
% Once filtered, inverse short-time Fourier transform is applied, to 
% reconstruct only the detections. The output signal can be displayed as a
% waveform, or listen to.
%
%
%
% INPUTS
%   - x, input acoustic signal vector,
%   - fs, sampling frequency (Hz),
%   - window, size of the stft window (samples), 
%   - overlap, overlap between two consecutive windows in the stft (%),
%   - my_tracks, ridge detector detection output as a Nx1 cell, were each
%     cell has 3 lines: 
%           * my_tracks{n}(1,:) is the detection time vector of the nth detected tonal,
%           * my_tracks{n}(2,:) is the detection frequency vector of the nth detected tonal,
%           * my_tracks{n}(3,:) is the detection amplitude vector of the nth detected tonal,
%   - f_limit (optional), is the limitation of the desired reconstructed frequency range
%           * f_limit = [fmin fmax]; (Hz).
%
% OUTPUTS
%   - x_reconstructed, the reconstructed acoustic signal,
%   - tx, the associated time (s),



function [x_reconstructed, tx] = sound_reconstruction_stft(x, fs, window, overlap, my_tracks,f_limit)
% Test the number of inputs
if nargin == 5,
    f_limit = [0 fs/2];
elseif nargin < 5
    disp('Input missing')
end

% Compute the spectrogram with the indicated parameters
[stft,f,t,~] = spectrogram(x,hann(window),round((overlap/100)*window),window,fs);

% Initialization of the detection matrix
% This matrix is a binary matrix of the same size as the stft matrix.
% we want to have zeros everywhere except for the detections
detection_matrix = zeros(length(f),length(t)); % set the matrix to zeros

% Go through ALL the detections
for i = 1:length(my_tracks) 
     detected_tonal_time =  my_tracks{i}(1,:); % time vector of the ith detected tonal
     detected_tonal_freq = my_tracks{i}(2,:); % frequency vector of the ith detected tonal
     
    if (mean(detected_tonal_freq)>f_limit(1)) && (mean(detected_tonal_freq)<f_limit(2)), % a limited bandwidth can be set using f_limit        
        % Find the start and end indices of the detected tonal in stft time
        % axis
        start_time_tf_index = find(t >= detected_tonal_time(1),1);
        end_time_tf_index = find(t >= detected_tonal_time(end),1);
        for tt = 0:end_time_tf_index-start_time_tf_index-1 % Go through time
            ff = find(f >= detected_tonal_freq(tt+1),1); % find the frequency bin of the corresponding to the detection
            detection_matrix(ff,start_time_tf_index+tt) = 1; % mark the coresponding (t,f) point as 1
        end
    end 
end

% Multiplication between the binary detection matrix and the stft
stft_filt = stft.*detection_matrix;
% ISTFT to retrieve the signal
[x_reconstructed, tx] = istft(stft_filt, hann(window), hann(window),round(((100-overlap)/100)*window), window, fs);




function [x, t] = istft(stft, awin, swin, hop, nfft, fs)
% Author: Ph.D. Eng. Hristo Zhivomirov        12/26/13 
% function: [x, t] = istft(stft, awin, swin, hop, nfft, fs)
% 
% Input:
% stft - STFT-matrix (only unique points, time
%        across columns, frequency across rows)
% awin - analysis window function
% swin - synthesis window function
% hop - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
%
% Output:
% x - signal in the time domain
% t - time vector, s
% signal length estimation and preallocation
L = size(stft, 2);          % determine the number of signal frames
wlen = length(swin);        % determine the length of the synthesis window
xlen = wlen + (L-1)*hop;    % estimate the length of the signal vector
x = zeros(1, xlen);         % preallocate the signal vector
% reconstruction of the whole spectrum
if rem(nfft, 2)             
    % odd nfft excludes Nyquist point
    X = [stft; conj(flipud(stft(2:end, :)))];
else                        
    % even nfft includes Nyquist point
    X = [stft; conj(flipud(stft(2:end-1, :)))];
end
% columnwise IFFT on the STFT-matrix
xw = real(ifft(X));
xw = xw(1:wlen, :);
% Weighted-OLA
for l = 1:L
    x(1+(l-1)*hop : wlen+(l-1)*hop) = x(1+(l-1)*hop : wlen+(l-1)*hop) + ...
                                      (xw(:, l).*swin)';
end
% scaling of the signal
W0 = sum(awin.*swin);                  
x = x.*hop/W0;                      
% generation of the time vector
t = (0:xlen-1)/fs;                 

