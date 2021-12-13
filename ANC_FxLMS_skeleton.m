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
sigLenSec = 2; 
fs_resample = 8000;
sigLenSample = sigLenSec*fs_resample;
speech_filename='audio_files/speech1.wav';
noise_filename='audio_files/White_noise1.wav';

% Read in the noise source (resample if necessary)
[noise_raw,fs_noise]=audioread(noise_filename);
noise_raw=noise_raw(1:fs_noise*sigLenSec);
noise=resample(noise_raw,fs_resample,fs_noise);
d = conv(noise, p);
d = d(1:size(noise,1));

% soundsc(noise,fs_resample);pause;
% soundsc(filt_noise,fs_resample);pause;

%%
% Plot the RIRs of the noise
% figure(1); clf;
% plot(RIR_noise(1:400,1));

% Plot the noisy signal
% figure(2); clf;
% plot(noise)


%% FxLMS
mu = 0.4;
delta = 5*10^-5;

x = noise;
xf = conv(x, h);

w = zeros(L, 1);
e = zeros(sigLenSample,1);
y = zeros(sigLenSample, 1);

tic
for n = L:sigLenSample

    x_samples = x(n:-1:n-L+1); % get last L samples of X in reverse order
  
    y(n,1) = w' * x_samples; % get y[n] by multiplying reversed x samples with w
    y_samples = y(n:-1:n-L+1,1); % get last L samples of y in reverse order
    y_h = h' * y_samples; % get y_h[n] by multiplying reversed y samples with h

    e(n) = y_h + d(n); % error is the sum of y_h[n] and d[n]

    xw = xf(n-L+1:n); % get the last L samples of xf (x filtered with h)
    w = w + ( mu / ((norm(xw))^2 + delta) ) * xw*e(n); % calculate new w

end
toc


%%
% Calculate the noise suppression
E_e = mean( e(round(0.1*sigLenSample):end,1) ); % approximate E(e), by taking the mean of 90% last values of e
E_d = mean( d(round(0.1*sigLenSample):end,1) ); % approximate E(d), by taking the mean of 90% last values of d

C = 10*log10((E_e.^2) ./ (E_d.^2));

%%
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation
% 
figure(1);
plot(d, 'DisplayName', "noise");
hold on;
plot(e, 'DisplayName', "error");
hold off;
legend;
title('Comparison of noise and error')
