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
load ("Computed_RIRs.mat")
M = 400;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length
J = size(RIR_sources,3);
p1 = RIR_noise(:,1);
p2 = RIR_noise(:,2);
h1 = RIR_sources(1:M, 1, :);
h1 = reshape(h1,M, J);
h2 = RIR_sources(1:M, 2, :);
h2 = reshape(h2,M, J);


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



% Run the cross cancel from session 1 to obtain the auralized speech
% Boolean flag to decide whether to add or not the auralized speec
addBinSig = 1;

if addBinSig
    load ("session1.mat", "g");
    % read in the speech signal
    [speech_raw,fs_speech]=audioread(speech_filename);
    speech_raw=speech_raw(1:fs_speech*sigLenSec);
    speech=resample(speech_raw,fs_resample,fs_speech);

    g = reshape(g, [], J);
    Lg = size(g, 1);

    binaural_sig = zeros(sigLenSample, 2);

end
%% Plot the RIRs of the noise
% figure(1); clf; 

%%


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

w = zeros(L,J); % Initialize adaptive filter
tic

%%
for n = L+M:sigLenSample
    if mod(n, 1000) == 0
        fprintf("%d out of %d\n",n,sigLenSample)
    end
    x_old=flip(x(n-L-M+2:n));
    X_Hmat=hankel(x_old(1:M),x_old(M:end));
   
    y=w'*X_Hmat'; % [M J]
    y1 = y;
    y2 = y;
 
    y_h1 = sum(h1'.*y1, 'all');
    y_h2 = sum(h2'.*y1, 'all');

    if addBinSig
        speech_old = flip(speech(n-Lg-M+2:n, :));
        Speech_Hmat=hankel(speech_old(1:Lg), speech_old(Lg:end));

        y_speech1 = (g' * Speech_Hmat)';
        y_speech2 = y_speech1;
         
        binaural_sig(n, 1) = sum(h1'.* y_speech1', 'all');
        binaural_sig(n, 2) = sum(h2'.* y_speech2', 'all');
        
    end

   
    e1(n)=d1(n) + sum(h1'.*y, 'all');
    e2(n)=d2(n) + sum(h2'.*y, 'all');

    e(n, :) = [e1(n) e2(n)];  
   
    xf1 = (h1' * X_Hmat)';
    xf2 = (h2' * X_Hmat)';
   
    for j=1:J
        xf = [xf1(:,j) xf2(:,j)];
        w(:,j) = w(:,j) - (mu/(norm(xf, 'fro')^2+delta))*xf*e(n,:)';
    end
end
toc

%%
if addBinSig
    T = min(5, sigLenSec);
    sig = binaural_sig(1: fs_resample*T)';
    n = d(1: fs_resample*T, 1);
    SNR = snr(sig(:,1), n);
    k = 10^(SNR/20);
    binaural_sig = binaural_sig / k;
end



%%
% Calculate the noise suppression
E_e = mean( e(round(0.3*sigLenSample):end,:) ); % approximate E(e), by taking the mean of 90% last values of e
E_d = mean( d(round(0.3*sigLenSample):end,:) ); % approximate E(d), by taking the mean of 90% last values of d

C1 = 10*log10((E_e(:,1).^2) ./ (E_d(:,1).^2));
C2 = 10*log10((E_e(:,2).^2) ./ (E_d(:,2).^2));

%%
% In the existing plot of the noisy signals, superimpose the corresponding
% error signals to appreciate the noise cancellation

figure(2);

subplot(2,1,1)
plot(d(:,1), 'DisplayName', "noise ");
hold on;
plot(e(:,1), 'DisplayName', "filtered signal");
title("noise and filtered signal channel 1")
hold off;
legend;

subplot(2,1,2)
plot(d(:,2), 'DisplayName', "noise");
hold on;
plot(e(:,2), 'DisplayName', "filtered signal");
title("noise and filtered signal channel 2")
hold off;

if addBinSig
    figure(3);
    
    subplot(2,1,1)   
    plot(d(:,1)+binaural_sig(:,1), 'DisplayName', "noise+speech");
    hold on;
    plot(e(:,1) + binaural_sig(:,1), 'DisplayName', "filtered signal");

    title("noise and filtered signal channel 1")
    hold off;
    legend;
    
    subplot(2,1,2)
    plot(d(:,2)+binaural_sig(:,2), 'DisplayName', "noise+speech");
    hold on;
    plot(e(:,2)+binaural_sig(:,2), 'DisplayName', "filtered signal");
    title("noise and filtered signal channel 2")
    hold off;
end

%%
listen = 0;
if listen
    disp("Playing speech:"); soundsc(speech, fs_resample); pause;
    disp("Playing noise:"); soundsc(d, fs_resample); pause;
    disp("Playing received signal:"); soundsc(d+binaural_sig, fs_resample); pause;

end