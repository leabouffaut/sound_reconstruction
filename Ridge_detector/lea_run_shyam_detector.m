function [my_tracks] = lea_run_shyam_detector(p,f,t,min_contour_len,max_contour_inactivity,SNR_thresh,min_intensity)

if nargin == 5,
    SNR_thresh = 4; % choose in the table in dt_SpectrogramRidgeDetector.m 
    min_intensity = -Inf;
end
if nargin == 6,
    min_intensity = -Inf;
end
dT = (t(2)-t(1)); % Time delta for the spectrogram
dF = f(2)-f(1); %Frequency delta for the spectrogram

h_tracker = dt_SpectrogramRidgeTracker(dT, dF, f);
% min_contour_len = 10; % mini nb of frames for a detection
% max_contour_inactivity = 4; %compensate the Lloyds mirror effect

% % Easy way :
%  TestSpectrogramContourTracker(p(:,1:1000), dT, f, -Inf, 10, 4)

% Set a global variable called my_tracks where all the detection data are
% going to be put in the form of 1 listed item per track. Each track has 3
% rows: time, freq and power.
global my_tracks ;
my_tracks = []; 

% Set parameters for the class
% h_tracker.XXX = XXX function of the class described in
% dt_Spectrogram_Ridge_Detector

h_tracker.Set_threshold_value(SNR_thresh);
if ~isinf(min_intensity)
    h_tracker.Set_min_intensity(min_intensity);
end
h_tracker.SetMinContourLength(min_contour_len);
h_tracker.SetMaxContourInactivity(max_contour_inactivity);
h_tracker.SetTrackingStartTime(t(1));
h_tracker.SetCallback(@GatherTracksCB, 0);

% Analysis starts
h_tracker.ProcessFrames(p);


% time_shyam = my_tracks{:}(1,:);
% f0_shyam = my_tracks{:}(2,:);
% ampl_shyam = my_tracks{:}(3,:);
 h_tracker.Flush();