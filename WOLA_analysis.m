% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2021
%
% The following is the skeleton code for the analysis stage of the WOLA method, which you need to
% complete


function [X,f] = WOLA_analysis(x,fs,window,nfft,noverlap)
%WOLA_analysis  short-time fourier transform
% INPUT:
%   x           : input time signal(s) (samples x channels)
%   fs          : sampling rate
%   window      : window function
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%   g           : transfer function (samples x channels)
%
% OUTPUT:
%   X           : STFT matrix (bins x frames x channels)
%   f           : frequency vector for bins


% use only half FFT spectrum
N_half = nfft / 2 + 1;

% get frequency vector
f = 0:(fs / 2) / (N_half - 1):fs / 2;

% init
L = floor((length(x) - nfft + (nfft / noverlap)) / (nfft / noverlap)); % number of frames
M = size(x,2); % number of channels
X = zeros(N_half, L, M);
hop = nfft - nfft / noverlap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code for binaural %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G=fft(g,nfft); %512 5
% G_half=G(1:N_half,:); %Use only half FFT spectrum

for m = 0:M-1 % channel index
    for l = 0:L-1 % Frame index

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete (3 - 5 lines) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        index_start=(l*(nfft-nfft/noverlap))+1;
        index_end=index_start+nfft-1;
        x_window=x(index_start:index_end,m+1).*window;
        X_fft=fft(x_window,nfft);
%         X(:,l+1,m+1)=X_fft(1:N_half).*G_half(:,m+1);
        X(:,l+1,m+1)=X_fft(1:N_half);
        %temp=X_fft.*G(:,m+1); %equivalent with line 51
        %X(:,l+1,m+1)=temp(1:N_half); %equivalent too
        
    end
end

end
