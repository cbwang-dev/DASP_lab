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
Lx = nfft - Lh + 1;

H = fft(h, nfft);   % DFT of h
H = H(:);           % make into a column vector (speed)
nx = length(x);     % total length of x
y = zeros(nx,1);

istart = 1;
y_s_old = zeros(nfft, 1);
while istart <= nx-Lx

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section of code to complete (5 - 10 lines) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    s = x(istart: istart + Lx -1, :);

    S = fft(s, nfft);
    Y_s = S.*H;
    y_s = ifft(Y_s, nfft); 


    y(istart: istart+Lh-2, :) = y_s(1:Lh-1, :) + y_s_old(end-Lh+2:end, :); % overlap
    
    core_len = nfft - 2*Lh + 2;
    y(istart+Lh-1: istart + Lh + core_len-2, :) = y_s(Lh : Lh + core_len-1); % core block

    y_s_old = y_s;
    
    istart = istart + Lx;

end


end
