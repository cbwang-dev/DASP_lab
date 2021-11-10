% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2021
%
% The following is the skeleton code for the synthesis stage of the WOLA method, which you need to
% complete


function x_synthesized = WOLA_synthesis(X,window,nfft,noverlap)
%WOLA_synthesis inverse short-time fourier transform.
%
% INPUT:
%   X             : input matrix (bins x frames x channels)
%   window        : window function
%   nfft          : FFT size
%   noverlap      : frame overlap; default: 2 (50%)
%
% OUTPUT:
%   x_synthesized : output time signal(s)


L = size(X,2); %frames
M = size(X,3); %channels

% ## Perform IFFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1 - 3 lines) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[X;conj(flip(X(1:nfft/2,:,:)))];
x=ifft(X,nfft);

% ## Apply synthesis window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1 - 3 lines) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x.*window;
x = permute(x,[1 3 2]); %(bins x channels x frames)

% ## Obtain re-synthesised signals

x_synthesized = zeros(L*nfft/noverlap+(nfft-nfft/noverlap),M);

for m = 0:L-1 % frame index

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section of code to complete (1 - 3 lines) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    index_start=m*nfft/noverlap+1;
    index_end=index_start+nfft-1;
    x_synthesized(index_start:index_end,:)=x_synthesized(index_start:index_end,:)+x(:,:,m+1);
    
end

end