% Session 4 - Exercise 1 - Single-channel FxNLMS ANC
% 
% Main points:
% (1) Generate the noisy microphone signal at the left ear
% (2) Implement the FxNLMS ANC.
% (3) Compute input/output SNRs and SNR improvement.
% (4) Implement the filter for left and right ears. Hear the synthesised signal

clear all;
% close all

% Load RIRs
load ComputedRIRs.mat
% Set length
sigLenSec = 5; 
RIRtruncSamples = 500;
fs_resample = 8000;
speech_filename='audio_files/speech1.wav';
noise_filename='audio_files/White_noise1.wav';

%%
% Plot the RIRs of the noise
figure(1); clf;
plot(RIR_noise(1:400,1));

% Read in the noise source (resample if necessary)
[noise_raw,fs_noise]=audioread(noise_filename);
noise_raw=noise_raw(1:fs_noise*sigLenSec);
noise=resample(noise_raw,fs_resample,fs_noise);

% Plot the noisy signal
figure(2); clf;
plot(noise)

%% FxLMS

M = 400;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length

mu = 0.5;   % Step size

W = zeros(L,1); % Initialize adaptive filter
    
tic
for n = 1:sigLenSample

    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation)
    X_Hmat = 

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y
    y =
    
    % STEP 3: Compute the error signal e(n)
    e(n) =

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf
    xf = 
    
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w
    w = 
end
toc


%%
% Calculate the noise suppression


%%
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation

figure(2); hold on;
