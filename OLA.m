% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2021
%
% The following is the skeleton code for the OLA function, which you need to
% complete.

function y = OLA(x,h,nfft)
%
% Overlap and add method for computing convolution between the filter, h,
% and signal x.
% INPUT:
%   x           : input time-domain signal(s) (samples x 1)
%   h           : filter (samples x 1)
%   nfft        : FFT size
%
% OUTPUT:
%   y           : convolved output signal

narginchk(2,3)

Lh = length(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1-2) lines %
% Calculate the appropriate value of Lx, the frame length of x,
% given your nfft and Lh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = nfft - Lh +1;

H = fft(h, nfft);   % DFT of h
H = H(:);           % make into a column vector (speed)
nx = length(x);     % total length of x
y = zeros(nx,1);

istart = 1;

while istart <= nx
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section of code to complete (5 - 10 lines) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iend=min(istart+Lx-1,nx);
    X_partial=fft(x(istart:iend),nfft);
    y_partial=ifft(X_partial.*H, nfft);
    iend=min(istart+nfft-1,nx);
    if iend==nx 
        % to get the size of output consistant so that the output shape
        % is determined, less bugs
        y(istart:iend)=y(istart:iend)+y_partial(1:iend-istart+1);
    else % normal situation
        y(istart:iend)=y(istart:iend)+y_partial;
    end
    istart=istart+Lx;
    
end


end
