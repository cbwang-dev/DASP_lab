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
M = 150;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length
p = RIR_noise(1:1000,1);
h = RIR_sources(1:M,1,1); % choose first loud speaker


% Set length
sigLenSec = 10; 
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
figure(1);
subplot(2,1,1);
plot(p);
title('p - RIR primary path');

subplot(2,1,2);
plot(h);
title('h - RIR secundary path');


%% FxLMS
mu = 0.5;
delta = 5*10^-5;

x = noise;
xf = conv(x, h);

w = zeros(L, 1);
e = zeros(sigLenSample,1);
y = zeros(sigLenSample, 1);

tic
for n = L+M:sigLenSample
    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation)
    x_old=flip(x(n-L-M+2:n));
    X_Hmat=hankel(x_old(1:M),x_old(M:end));
    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y
    y=w'*X_Hmat'; % [1 M]
    % STEP 3: Compute the error signal e(n)
    e(n)=d(n)+h'*y'; 
    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf
    xf = (h' * X_Hmat)';
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w
    w = w - (mu/(norm(xf)^2+delta))*xf*e(n);
end
toc


%%
% Calculate the noise suppression
E_e = mean( e(round(0.3*sigLenSample):end,1) ); % approximate E(e), by taking the mean of 90% last values of e
E_d = mean( d(round(0.3*sigLenSample):end,1) ); % approximate E(d), by taking the mean of 90% last values of d

C = 10*log10((E_e.^2) ./ (E_d.^2))

%%
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation
% 
figure(2);
% plot(e+xf, 'DisplayName', 'Result')
plot(d, 'DisplayName', "noise");
hold on;
plot(e, 'DisplayName', "error");
hold off;
legend;
title('Comparison of noise and error')
