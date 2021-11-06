% Lab 1 for Digital Audio Signal Processing Lab Sessions
% Exercise 1-4: 3D audio
% 
% In this lab, we derive a set of filters g that can be used, along with
% the measured RIRs H, to produce the proper psychocoustic rendition of 3D
% audio
%
% Use MySA_GUI load /GUI_setups/1_4.mat first before experiment!!!

clear;
% close all

% Load ATFs
load Computed_RIRs

% Load measured HRTFs
load HRTF 

% Define the signal length
siglength = 10;

% Load the speech signal and resample
flag_resample_8k = 1;
source_filename{1} = 'audio_files/speech1.wav';
[source,fs] = audioread(source_filename{1}); 
source = source(1:fs*siglength); %truncate
if flag_resample_8k == 1 %resample if necessary
    source = resample(source,8000,fs);
end
  
        
% Noise flag for the noise perturbing SOE
noiseFlag = 0;
% Noise flag for sweetspot tests
sweetspotFlag = 1;

% Number of loudspeakers
J = size(RIR_sources,3);

% Define the lengths of RIRs and g filters
Lh = 400; % Length of impulse responses of listening room
% Lh = 1500; % Length of impulse responses of listening room
Lg = ceil((2*Lh-2)/(J-2));   % Length of filters g_j

% Truncate impulse response to reduce computational complexity
RIR_sources = RIR_sources(1:Lh,:,:);

% Calculate delay for SOE
Delta=ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);

% Define the Toeplitz matrices for left and right ear (left side of SOE)
HL=[]; 
for j=1:J
    part_toeplitz=toeplitz(RIR_sources(:,1,j),zeros(Lg-1,1));
    HL=[HL part_toeplitz];
end
HR=[];  
for j=1:J
    part_toeplitz=toeplitz(RIR_sources(:,2,j),zeros(Lg-1,1));
    HR=[HR part_toeplitz];
end

% Define the HRTFs for left and right ear (right side of SOE) from the
% loaded HRTF
xL = HRTF(:,1); % Left ear
xR = HRTF(:,2); % Right ear
% delayed version of xL and xR 
% Question: need to modify siglength regarding the loss introduced by delay?
% This version, xL(Delta+1,end) info lost.
xL_delayed = delayseq(xL,Delta);
xR_delayed = delayseq(xR,Delta);
xL_delayed = xL_delayed(1:Lh); %also need to truncate
xR_delayed = xR_delayed(1:Lh); %also need to truncate

% Construct H (from HL and HR) and x (from xL and xR) and remove all-zero 
% rows in H, and the corresponding elements in x
H = [HL;HR];
% x = [xL;xR];
x = [xL_delayed;xR_delayed];

index_nonzero = sum(abs(H),2)>0;
H_nonzero = H(index_nonzero,:);
x_nonzero = x(index_nonzero,:);

% Solve the SOE
if ~noiseFlag
    g=H_nonzero\x_nonzero;
else
    stddev=0.05*std(H_nonzero(:,1)); 
    H_nonzero=H_nonzero+wgn(size(H_nonzero,1),size(H_nonzero,2),stddev);
    g=H_nonzero\x_nonzero;
end

% Plot estimated and real HRTFs
figure(1);title("Estimated and Real HRTFs without noise");hold on;
subplot(3,1,1);plot(1:size(x_nonzero,1),x_nonzero,'color','red');
legend("x");
title("Real HRTFs wothout noise");
subplot(3,1,2);plot(1:size(x_nonzero,1),H_nonzero*g,'color','blue');
legend("H*g");
title("Estimated HRTFs wothout noise");
subplot(3,1,3);plot(1:size(x_nonzero,1),H_nonzero*g-x_nonzero);
legend("H*g-x")
title("Difference between Real and Estimated HRTFs");

% Calculate synthesis error
synth_error=norm(H_nonzero*g-x_nonzero);
fprintf("Synthesis error is %d.(1.4.5 scenario)\n",synth_error);

% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x
% binaural speech synthesized with H and g
estimated_HRTF_total=H*g;
estimated_HRTF_left =estimated_HRTF_total(1:Lh);
estimated_HRTF_right=estimated_HRTF_total(Lh+1:end);
estimated_left=conv(source,estimated_HRTF_left);
estimated_right=conv(source,estimated_HRTF_right);
estimated_binaural=[estimated_left,estimated_right];
% binaural speech synthesized with x (or, HRTF)
measured_left=conv(source,xL_delayed(1:Lh));
measured_right=conv(source,xR_delayed(1:Lh));
measured_binaural=[measured_left,measured_right];
% Calculate synthesis error
synth_error=norm(estimated_binaural-measured_binaural);
fprintf("Synthesis error is %d.(1.4.6 scenario)\n",synth_error);
% audio test
% soundsc(estimated_binaural,8000);pause;soundsc(measured_binaural,8000);

% Answer to 1.4.8 TODO

% Answer to 1.4.10 (sweet spot) By manipulating m_pos and apply in 
% create_rirs.m
if sweetspotFlag == 1
    synth_errors=[];
    for i = 1:10
%         m_pos(:,1)=m_pos(:,1)+ones(2,1)*0.01;
        m_pos(:,2)=m_pos(:,2)+ones(2,1)*0.01;
        [RIR_sources,~]=create_rirs(m_pos,s_pos,v_pos,room_dim,rev_time,fs_RIR,Lh);
        HL=[]; 
        for j=1:J
            part_toeplitz=toeplitz(RIR_sources(:,1,j),zeros(Lg-1,1));
            HL=[HL part_toeplitz];
        end
        HR=[];  
        for j=1:J
            part_toeplitz=toeplitz(RIR_sources(:,2,j),zeros(Lg-1,1));
            HR=[HR part_toeplitz];
        end
        H = [HL;HR];
        estimated_HRTF_total=H*g;
        estimated_HRTF_left =estimated_HRTF_total(1:Lh);
        estimated_HRTF_right=estimated_HRTF_total(Lh+1:end);
        estimated_left=conv(source,estimated_HRTF_left);
        estimated_right=conv(source,estimated_HRTF_right);
        estimated_binaural_new=[estimated_left,estimated_right];
        synth_errors=[synth_errors norm(estimated_binaural-estimated_binaural_new)];
    end
    figure
    plot(1:10,synth_errors)
    title('Synthesis Error with the Movement of Both Mics')
    xlabel('cm')
    % TODO: test mic movement in another direction
end
