% Lab 1 for Digital Audio Signal Processing Lab Sessions
% Exercise 1-4: 3D audio
% 
% In this lab, we derive a set of filters g that can be used, along with
% the measured RIRs H, to produce the proper psychocoustic rendition of 3D
% audio
%
%

clear;
% close all

% Load ATFs
load Computed_RIRs

% Load measured HRTFs
load HRTF 

% Define the signal length
siglength = 10;

% Load the speech signal and resample
source_filename{1} = 'audio_files/speech1.wav';
[source,fs] = audioread(source_filename{1}); 
source = source(1:fs*siglength);


% Noise flag for the noise perturbing SOE
noiseFlag = 0;
% Noise flag for sweetspot tests
sweetspotFlag = 0;

% Number of loudspeakers
J = size(RIR_sources,3);

% Define the lengths of RIRs and g filters
Lh = 400;              % Length of impulse responses of listeniLg room
Lg = (2*Lh-2)/(J-2);   % Length of filters g_j

% Truncate impulse response to reduce computational complexity
RIR_sources = RIR_sources(1:Lh,:,:);

% Calculate delay for SOE
Delta=ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);
% RIR_sources = delayseq(RIR_sources,Delta); % delayed sequence of RIR

% Define the Toeplitz matrices for left and right ear (left side of SOE)
HL=[]; 
for j=1:J
    first_column=[RIR_sources(:,1,j)' zeros(1,Lg-1)];
    first_row=[RIR_sources(1,1,j) zeros(1,Lg-1)];
    part_toeplitz=toeplitz(first_column,first_row);
    HL = [HL part_toeplitz];
end
HR=[];  
for j=1:J
    first_column=[RIR_sources(:,2,j)' zeros(1,Lg-1)];
    first_row=[RIR_sources(1,2,j) zeros(1,Lg-1)];
    part_toeplitz=toeplitz(first_column,first_row);
    HR = [HR part_toeplitz];
end

% Define the HRTFs for left and right ear (right side of SOE) from the
% loaded HRTF
xL = [1 zeros(1,Lh+Lg-2)]'; % Left ear
xR = [1 zeros(1,Lh+Lg-2)]'; % Right ear

% Construct H (from HL and HR) and x (from xL and xR) and remove all-zero 
% rows in H, and the corresponding elements in x
H = [HL;HR];
x = [xL;xR];

% Solve the SOE
if ~noiseFlag
    g=H\x;
else
    % With noise
end

% Plot estimated and real HRTFs
figure(1); clf;

% Calculate synthesis error

% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x


