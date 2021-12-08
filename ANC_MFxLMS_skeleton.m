% Session 4 - Exercise 2 - Multi-channel FxNLMS ANC
%           - Exercise 3 - Multi-channel FxNLMS ANC in 3D audio scenario
% 
%
% Main points:
% (1) Generate the noisy microphone signals
% (2) Implement the ANC.
% (3) Compute input/output SNRs and SNR improvement
% (4) Implement the filter for left and right ears. Hear the synthesised signal

clearvars;
close all;
% Load RIRs
load Computed_RIRs.mat
M = 400;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length
p = RIR_noise(:,:);
h = RIR_sources(1:M,:,:);
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

filt_noise = conv2(noise, p);


% Run the cross cancel from session 1 to obtain the auralized speech

%% Plot the RIRs of the noise
figure(1); clf; 

%%

% Read in the noise source (resample if necessary)

% Boolean flag to decide whether to add or not thee auralized speec
% addBinSig = 1;

% Plot the noisy signals (1 for each subplot)
figure(2); clf;


%% MFxLMS

M = 400;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length
J = 5;
mu = 0.5;   % Step size
y_h = zeros(sigLenSample, 2, J);
e = zeros(sigLenSample, 2);
xf = zeros(sigLenSample, 2, J);
w = zeros(L,J); % Initialize adaptive filter
y = zeros(L, J);
x = noise;
tic
for n=L:sigLenSample
   
    % STEP 1 : Arrange the previous L−1 samples of x(n) up until the
    % current sample, n, in an [Lx1] vector
    samples_x = x(n-L+1:n, :);

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y
    y(n, :) = flip(samples_x)' * w;
    samples_y = y(n-M+1:n, :);
    y_h(n, 1, :) = flip(samples_y)' * reshape(h(:, 1, :), M, J); 
    y_h(n, 2, :) = flip(samples_y)' * reshape(h(:, 2, :), M, J); 

    % STEP 3: Compute the error signal e(n)
    e(n, :) =  filt_noise(n,:) + sum(y_h(n,:,:) ,3);

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf
    xf(n, 1, :) = flip(samples_x)' * h(:, 1, :);
    xf(n, 2, :) = flip(samples_x)' * h(:, 2, :);

    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w
    for j = 1:J
        w(:, j) = w(:, j) - (mu / (norm(xf(n,:,j), 'fro')^2 + delta)) * xf(n,:, j)*e(n,:)';
    end
end
toc

%%
% Calculate the noise suppression

%%
% In the existing plot of the noisy signals, superimpose the corresponding
% error signals to appreciate the noise cancellation

figure();
plot(noise);
hold on;
plot(e(:,1));
plot(e(:,2));
title("error e")
hold off;
