function S = STFT(signal, win, hopSize, F, Fs)
%   Q5:
%   implementation of STFT using the method we developed on Q2.

signal = signal(:).';  %set signal as a row vector.

if (length(win) == 1)
    win  =  hamming(win);  % calculate the window
end

if (length(F) == 1)
    omega = (2*pi/F)*(0:F-1);  % calculate the frequencies
else
    omega = F;
end

N = length(signal);  % length of the otiginal signal
N_downsampled = floor(N/hopSize);  % length of the downsampled signal
X = zeros([length(omega), N]);
X_downsample = zeros([hopSize, N_downsampled]);
S = zeros([length(omega), N_downsampled]);

for k = 1:length(omega)  % calculating for each frequency:
    X(k,:) = signal.*exp(-1i*omega(k)*(1:N));  % multiply with exponent
    for m = 1:hopSize  % polyphase implementation
          first_index = 1+hopSize-m;  % first index of the shifted signal
          last_index = 1+hopSize*N_downsampled-m;  % last index of the shifted signal
          X_downsample(m,:) = X(k,first_index:hopSize:last_index)*win(m);  % shift, downsample and filter
    end
    S(k,:) = sum(X_downsample);  % sum over the polyphase of the given frequency
end

S = fftshift(S,1);  % shift zero frequency to the center of the spectrum
end
