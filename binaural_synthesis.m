load HRTF.mat;
%% first part
% Read out the wav file speech1.wav, and resample it to 8 kHz. 
% Truncate the signal to a length of 10 s and store it in the variable x.
speechfilename = 'audio_files/speech1.wav';
[x,Fs] = audioread(speechfilename);
resample_frequency = 8000;
x = resample(x,resample_frequency,Fs);
desired_len = 10; % in seconds
x = x(1:desired_len * resample_frequency,:,:);
% soundsc(x,resample_frequency);

% create several binaural signals
binaural_sig1_1=[x x];
binaural_sig2_1=[x 0.5*x];
x_delay3=[0;0;0;x(1:size(x,1)-3)];
binaural_sig3_1=[x x_delay3];
binaural_sig4_1=[fftfilt(HRTF(:,1),x) fftfilt(HRTF(:,2),x)];

%% second part
% Read out the wav file speech1.wav, and resample it to 8 kHz. 
% Truncate the signal to a length of 10 s and store it in the variable x.

speechfilename = 'audio_files/speech2.wav';
[x,Fs] = audioread(speechfilename);
resample_frequency = 8000;
x = resample(x,resample_frequency,Fs);
desired_len = 10; % in seconds
x = x(1:desired_len * resample_frequency,:,:);
% soundsc(x,resample_frequency);

%% create several binaural signals
binaural_sig1_2=[x x];
binaural_sig2_2=[0.5*x x];
x_delay3=[0;0;0;x(1:size(x,1)-3)];
binaural_sig3_2=[x_delay3 x];
binaural_sig4_2=[fftfilt(HRTF(:,2),x) fftfilt(HRTF(:,1),x)];

%% combination
binaural_sig1=binaural_sig1_1+binaural_sig1_2;
binaural_sig2=binaural_sig2_1+binaural_sig2_2;
binaural_sig3=binaural_sig3_1+binaural_sig3_2;
binaural_sig4=binaural_sig4_1+binaural_sig4_2;

%% sound test
% soundsc(binaural_sig1,resample_frequency);pause;
% soundsc(binaural_sig2,resample_frequency);pause;
% soundsc(binaural_sig3,resample_frequency);pause;
% soundsc(binaural_sig4,resample_frequency);pause;