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
J = size(RIR_sources,3);
p1 = RIR_noise(:,1);
p2 = RIR_noise(:,2);
h1 = RIR_sources(1:L, 1, :);
h1 = reshape(h1,L, J);
h2 = RIR_sources(1:L, 2, :);
h2 = reshape(h2,L, J);


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
d1 = conv(noise, p1);
d1 = d1(1:size(noise,1));
d2 = conv(noise, p2);
d2 = d2(1:size(noise,1));
d = [d1 d2];

% soundsc(noise,fs_resample);pause;
% soundsc(filt_noise,fs_resample);pause;


% Run the cross cancel from session 1 to obtain the auralized speech

%% Plot the RIRs of the noise
% figure(1); clf; 

%%

% Read in the noise source (resample if necessary)

% Boolean flag to decide whether to add or not thee auralized speec
% addBinSig = 1;

% Plot the noisy signals (1 for each subplot)
% figure(2); clf;


%% MFxLMS


mu = 0.5;   % Step size
delta = 10^-5; % regularization factor

y = zeros(L, J); % vector of x samples filetered by w
y_h = zeros(sigLenSample, 2, J); % vector of y samples filtered by h

e1 = zeros(sigLenSample, 1); % channel1 error
e2 = zeros(sigLenSample, 1); % channel2 error
e = zeros(sigLenSample, 2); % dual channel error

x = noise;
xf1 = zeros(sigLenSample, J); % xbar of channel 1
xf2 = zeros(sigLenSample, J); % xbar of channel 2

for i=1:J % convolve xf1 and xf2 with h
    xf_tmp = conv(x, h1(:,i));
    xf1(:,i) = xf_tmp(1:sigLenSample,1);
    
    xf_tmp = conv(x, h2(:,i));
    xf2(:,i) = xf_tmp(1:sigLenSample,1);
end

w = zeros(L,J); % Initialize adaptive filter
tic

%%
for n=L:sigLenSample
   
   
    samples_x = x(n:-1:n-L+1, 1);

    y(n, :) = w' * samples_x;
    samples_y = y(n:-1:n-M+1, :);

    y_h1 = sum(h1 .* samples_y); 
    % multiply reversed y samples from each channel with respective h coefficients element wise 
    % and sum them column wise in order to receive vector y_h1

    y_h2 = sum(h2 .* samples_y); 
    % analogous to y_h1

    e1(n,1) = d1(n,1) + sum(y_h1, 'all'); % sum all y_h1 values and add them to e1 to get the error
    e2(n,1) = d2(n,1) + sum(y_h2, 'all'); % same as e1
    e(n, :) = [e1(n,1) e2(n,1)];
   
    for j = 1:J
        Xj = [xf1(n-L+1:n,1) xf2(n-L+1:n,1)]; % construct Xj matrix from L last samples of xf1 and xf2
        w(:, j) = w(:, j) - (mu / (norm(Xj, 'fro')^2 + delta)) * Xj*e(n,:)'; % update w
    end
end
toc

%%
% Calculate the noise suppression
E_e = mean( e(round(0.1*sigLenSample):end,:) ); % approximate E(e), by taking the mean of 90% last values of e
E_d = mean( d(round(0.1*sigLenSample):end,:) ); % approximate E(d), by taking the mean of 90% last values of d

C1 = 10*log10((E_e(:,1).^2) ./ (E_d(:,1).^2));
C2 = 10*log10((E_e(:,2).^2) ./ (E_d(:,2).^2));

%%
% In the existing plot of the noisy signals, superimpose the corresponding
% error signals to appreciate the noise cancellation

figure(2);

subplot(2,1,1)
plot(d(:,1), 'DisplayName', "noise");
hold on;
plot(e(:,1), 'DisplayName', "error");
title("d-e channel 1")
hold off;
legend;

subplot(2,1,2)
plot(d(:,2), 'DisplayName', "noise");
hold on;
plot(e(:,2), 'DisplayName', "error");
title("d-e channel 2")
hold off;
