%% load RIR results
load('Computed_RIRs.mat')

%% define filenames of speech and noise
speechfilename{1} = 'audio_files/speech1.wav';
speechfilename{2} = 'audio_files/speech2.wav';
noisefilename = 'audio_files/Babble_noise1.wav';

%% load target audios
[s1,Fs1] = audioread(speechfilename{1}); % for first audio file
[s2,Fs2] = audioread(speechfilename{2}); % for second audio file
[n1,Fn1] = audioread(noisefilename); % for noise signal

%% set desired length (in second) of recorded microphone signals
desired_len = 5; % in seconds
s1 = s1(1:desired_len * Fs1,:,:);
s2 = s2(1:desired_len * Fs2,:,:);
n1 = n1(1:desired_len * Fn1,:,:);

%% resample
s1 = resample(s1,fs_RIR,Fs1);
s2 = resample(s2,fs_RIR,Fs2);
n1 = resample(n1,fs_RIR,Fn1);

res_s1 = fftfilt(RIR_sources(:,:,1), s1);
res_s2 = fftfilt(RIR_sources(:,:,2), s2);
res_n1 = fftfilt(RIR_noise, n1);

% mic = res_s1;
% mic = res_s1 + res_s2;
mic = res_s1 + res_s2 + res_n1;

%% save the microphone signals
save('mic.mat','mic','fs_RIR')

%% play the microphone signals that are recorded by the microphone array 
soundsc(mic(:,1), fs_RIR);
pause;
soundsc(mic(:,2), fs_RIR);

figure;
plot(mic(:,1),'r')
hold on
plot(mic(:,2),'g')
hold on;
mic_diff = mic(:,1) - mic(:,2);
plot(mic_diff,'y')