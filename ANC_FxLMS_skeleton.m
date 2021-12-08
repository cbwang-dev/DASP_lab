% Session 4 - Exercise 1 - Single-channel FxNLMS ANC
% 
% Main points:
% (1) Generate the noisy microphone signal at the left ear
% (2) Implement the FxNLMS ANC.
% (3) Compute input/output SNRs and SNR improvement.
% (4) Implement the filter for left and right ears. Hear the synthesised signal

clearvars;
close all;
% Load RIRs
load Computed_RIRs.mat
M = 400;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length
p = RIR_noise(:,1);
h = RIR_sources(1:M,1,1);
% Set length
sigLenSec = 5; 
RIRtruncSamples = 500;
fs_resample = 8000;
sigLenSample = sigLenSec*fs_resample;
speech_filename='audio_files/speech1.wav';
noise_filename='audio_files/White_noise1.wav';

%%
% Plot the RIRs of the noise
% figure(1); clf;
% plot(RIR_noise(1:400,1));

% Read in the noise source (resample if necessary)
[noise_raw,fs_noise]=audioread(noise_filename);
noise_raw=noise_raw(1:fs_noise*sigLenSec);
noise=resample(noise_raw,fs_resample,fs_noise);

% Plot the noisy signal
% figure(2); clf;
% plot(noise)

filt_noise = conv(noise, p);
filt_noise = filt_noise(1:size(noise,1));
%% FxLMS



mu = 0.5;   % Step size
delta = 5*10^-5;

w = zeros(L,1); % Initialize adaptive filter
x = noise;
d = filt_noise;
y = zeros(sigLenSample, 1);
y_h = zeros(sigLenSample, 1);
e = zeros(sigLenSample, 1);
xf = zeros(sigLenSample, 1);
tic
for n = L:sigLenSample

    % STEP 1 : Arrange the previous L−1 samples of x(n) up until the
    % current sample, n, in an [Lx1] vector
    len_x_samples = max(M,L);
    samples_x = x(n-len_x_samples+1:n, 1);

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y
    y(n, 1) = w' * flip(samples_x);
    samples_y = y(n-M+1:n, 1);
    y_h(n, 1) = h'*flip(samples_y);
    
    % STEP 3: Compute the error signal e(n)
    e(n, 1) =  filt_noise(n,1) + y_h(n,1);

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf
    xf(n) = h' * flip(samples_x);
    
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w
    w = w - (mu / (norm(xf(n,1))^2 + delta)) * xf(n,1)*e(n,1);
end
toc


%%
% Calculate the noise suppression
C = 10*log10((e.^2) ./ (d.^2));

%%
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation

figure();
plot(noise);
hold on;
plot(e);
title("error e")
hold off;
