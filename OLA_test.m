clear;
close all;
siglength=5;
nfft=512;
J=5;
dirac = zeros(1,50);
dirac(1) = 1;
[source,fs] =audioread('audio_files/speech1.wav');
source = source(1:siglength*fs);
dirac_speech = OLA(source,dirac,nfft);

figure(1);
subplot(2,1,1);plot(source,'color','red');legend('speech');
title("original speech signal");
subplot(2,1,2);plot(dirac_speech,'color','blue');legend('dirac speech');
title("filtered speech signal");
synth_error=norm(dirac_speech-source);
fprintf("Synthesis error is %d.(OLA Dirac scenario)\n",synth_error);

%% binaural test
% plain convolution inherited from SOE.m
load session1 % g H
Lh=400;
conv_HRTF_total=H*g;
conv_HRTF_left =conv_HRTF_total(1:Lh);
conv_HRTF_right=conv_HRTF_total(Lh+1:end);
conv_left=conv(source,conv_HRTF_left);
conv_right=conv(source,conv_HRTF_right);
conv_binaural=[conv_left,conv_right];

% OLA implementation
load Computed_RIRs.mat
% the dimension of filters should be smaller than nfft
RIR_sources=RIR_sources(1:nfft/2,:,:);
g=reshape(g,[size(g,1)/J,J]); %shape [265,5], modify if nfft<265

% test for filtering with g, for comparison with WOLA
OLA_g=zeros(size(source,1),J);
for j=1:J
    OLA_g(:,j)=OLA(source,g(:,j),nfft);
end
% OLA_filter_g=mean(OLA_g,2);
% disp("playing OLA filtered g");soundsc(OLA_filter_g,fs);pause;
% disp("playing OLA filtered g speaker 1");soundsc(OLA_g(:,1),fs);pause;
% disp("playing OLA filtered g speaker 5");soundsc(OLA_g(:,5),fs);pause;
% save("OLA_filter_g.mat","OLA_filter_g");

OLA_binaural_left=zeros(size(source,1),J);
OLA_binaural_right=zeros(size(source,1),J);
for j=1:J
    filtered_g_source=OLA(source,g(:,j),nfft);
    OLA_binaural_left (:,j)=OLA(filtered_g_source,RIR_sources(:,1,j),nfft);
    OLA_binaural_right(:,j)=OLA(filtered_g_source,RIR_sources(:,2,j),nfft);
end
OLA_binaural=[mean(OLA_binaural_left,2) mean(OLA_binaural_right,2)];

% Now listen:
% disp("playing SOE estimated binaural (from session 1)");load SOE_binaural_sound_fs;soundsc(estimated_binaural,fs);pause;
% disp("playing convolution based binaural");soundsc(conv_binaural,fs);pause;
% disp("playing OLA based binaural");soundsc(OLA_binaural,fs);

% Now look:
figure(2);
subplot(2,1,1);plot(conv_left);;hold on;
plot(conv_right);
title("convolution based binaural");
legend("left", "right");
subplot(2,1,2);plot(OLA_binaural(:,1));hold on;
plot(OLA_binaural(:,2));title("OLA based binaural");
legend("left", "right");
save("OLA_binaural_sound.mat","OLA_binaural")