clear;
close all;
%% hyperparameter and audio definition
[source, fs] =audioread('audio_files/speech1.wav');
siglength = 5;
source = source(1:siglength*fs);
nfft = 512;
noverlap = 2; %default
analysis  = sqrt(hann(nfft,'periodic'));
synthesis = sqrt(hann(nfft,'periodic')); 
J = 5;

load session1
g=reshape(g,[size(g,1)/J,J]); %shape [265,5], modify if nfft<265
source_duplicated = repmat(source,1,J);

[X_wola,f] = WOLA_analysis(source_duplicated,fs,analysis,nfft,noverlap,g);
x_WOLA = WOLA_synthesis(X_wola,synthesis,nfft,noverlap);

% sound test
% disp("playing WOLA filtered g");soundsc(mean(x_WOLA,2),fs);pause;
% disp("playing WOLA filtered g speaker 1");soundsc(x_WOLA(:,1),fs);pause;
% disp("playing WOLA filtered g speaker 5");soundsc(x_WOLA(:,5),fs);pause;
% load OLA_filter_g
% disp("playing OLA filtered g");soundsc(OLA_filter_g,fs);pause;

% the dimension of filters should be smaller than nfft
load('Computed_RIRs.mat')
RIR_sources=RIR_sources(1:nfft/2,:,:);
x_left=zeros(size(x_WOLA,1),J);
x_right=zeros(size(x_WOLA,1),J);

for j = 1:J
    x_left (:,j)=OLA(x_WOLA(:,j),RIR_sources(:,1,j),nfft);
    x_right(:,j)=OLA(x_WOLA(:,j),RIR_sources(:,2,j),nfft);
end
WOLA_binaural = [mean(x_left,2) mean(x_right,2)];

% Look!
figure(3);load SOE_binaural_sound_fs;
subplot(3,1,1);plot(WOLA_binaural(:,1));hold on;
plot(WOLA_binaural(:,2));title("WOLA based binaural");
legend("left", "right");
subplot(3,1,2);plot(estimated_binaural(:,1));hold on;
plot(estimated_binaural(:,2));title("SOE binaural");
legend("left", "right");
load OLA_binaural_sound;subplot(3,1,3);
plot(OLA_binaural(:,1));hold on;
plot(OLA_binaural(:,2));title("OLA based binaural");
legend("left", "right");


% Listen!
% load SOE_binaural_sound_fs;
% disp("playing SOE estimated binaural (from session 1)");soundsc(estimated_binaural,fs);pause;
% disp("playing WOLA estimated binaural");soundsc(WOLA_binaural,fs);pause;
% load OLA_binaural_sound;disp("playing OLA based binaural");soundsc(OLA_binaural,fs);
