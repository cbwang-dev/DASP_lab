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
M = 50;  % Length of secondary path (RIR)
L = M;  % Adaptive filter length
p = RIR_noise(1:100,:);
h = RIR_sources(1:M,:,1); % choose first loud speaker


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
filt_noise = fftfilt(p,noise);

% soundsc(noise,fs_resample);pause;
% soundsc(filt_noise,fs_resample);pause;

%%
% Plot the RIRs of the noise
% figure(1);
% subplot(2,1,1);
% plot(p(:,1));
% title('p - RIR primary path - left');
% subplot(2,1,2);
% plot(h(:,1));
% title('h - RIR secundary path - left');


%% FxLMS
mu = 0.5;
delta = 5*10^-5;

x = noise;

% w = zeros(L, 1);
e = zeros(sigLenSample,2);
% y = zeros(sigLenSample, 1);
for i=1:2
    w = zeros(L, 1);
    y = zeros(sigLenSample, 1);
    tic
    for n = L+M:sigLenSample
        % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
        % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
        % use for the filtering operation)
        x_old_flip=flip(x(n-L-M+2:n));
        X_Hmat=hankel(x_old_flip(1:M),x_old_flip(M:end));
        % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
        % generate loudspeaker signals y(n). Store the samples of y(n) into the
        % vector y
        y=w'*X_Hmat'; % [1 M]
        % STEP 3: Compute the error signal e(n)
        e(n,i)=filt_noise(n)+h(:,i)'*y'; 
        % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
        % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
        % the vector xf
        xf = (h(:,i)' * X_Hmat)';
        % STEP 5: Update the filter w(n). Store the samples of w(n) into
        % the vector w
        w = w - (mu/(norm(xf)^2+delta))*xf*e(n,i);
    end
    toc
end

%%
% Calculate the noise suppression
E_e = mean( e(round(0.3*sigLenSample):end,1) ); % approximate E(e), by taking the mean of 90% last values of e
E_d = mean( filt_noise(round(0.3*sigLenSample):end,1) ); % approximate E(d), by taking the mean of 90% last values of d

C = 10*log10((E_e.^2) ./ (E_d.^2))

%%
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation
% 
figure(2);
% plot(e+xf, 'DisplayName', 'Result')
subplot(2,1,1);
plot(filt_noise(:,1), 'DisplayName', "noise - left");
hold on;
plot(e(:,1), 'DisplayName', "error - left");
hold off;
legend;
title('Comparison of noise and error - left')
subplot(2,1,2);
plot(filt_noise(:,2), 'DisplayName', "noise - right");
hold on;
plot(e(:,2), 'DisplayName', "error - right");
hold off;
legend;
title('Comparison of noise and error - right')
