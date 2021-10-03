%% Main code
clc;clear all;close all;

%% Q.6
[touchtone_1, Fs] = audioread('touchtone_1.wav');
% sound(touchtone_1, Fs)
%% Q.7
N = length(touchtone_1);
T_vec = 0:1/Fs:(N-1)/Fs;    % time vector 
F_vec = linspace(-Fs/2,Fs/2,N); % frequency vector

X = fftshift(fft(touchtone_1));
figure;
plot(F_vec, abs(X));
title('Q.7: Fourier transform of the entire Signal');
xlabel('frequency[Hz]');ylabel('X(f)');
%% Q.8a

win = 512;
hopSize = 512;
F = 64;
S = STFT(touchtone_1, win, hopSize, F, Fs);

% plotting the spectrogram

S_size = size(S);
T_vec = linspace(0,(N-1)/Fs, S_size(2));
F_vec = linspace(-Fs/2, Fs/2, S_size(1));

figure;
imagesc(T_vec, F_vec, abs(S));
%show only right half of spectrum (symmetric)
axis ([0 (N-1)/Fs 0 Fs/2]);
axis xy;
xlabel('time[sec]');
ylabel('frequency[Hz]');
title(['Q.8a: STFT with win = ', num2str(win), ', hopsize = ', num2str(hopSize), ', F = ', num2str(F)]);

% figure;
% spectrogram(touchtone_1,512,[],64,'yaxis')



%% Q.8b

win = 256;
hopSize = 256;
F = 256;
S = STFT(touchtone_1, win, hopSize, F, Fs);

% plotting the spectrogram

S_size = size(S);
T_vec = linspace(0,(N-1)/Fs, S_size(2));
F_vec = linspace(-Fs/2, Fs/2, S_size(1));

figure;
imagesc(T_vec, F_vec, abs(S));
%show only right half of spectrum (symmetric)
axis ([0 (N-1)/Fs 0 Fs/2]);
axis xy;
xlabel('time[sec]');
ylabel('frequency[Hz]');
title(['Q.8b: STFT with win = ', num2str(win), ', hopsize = ', num2str(hopSize), ', F = ', num2str(F)]);

% figure;
% spectrogram(touchtone_1,256,[],256,'yaxis')

%% Q.8c

win = 1024;
hopSize = 1024;
F = 256;
S = STFT(touchtone_1, win, hopSize, F, Fs);

% plotting the spectrogram

S_size = size(S);
T_vec = linspace(0,(N-1)/Fs, S_size(2));
F_vec = linspace(-Fs/2, Fs/2, S_size(1));

figure;
imagesc(T_vec, F_vec, abs(S));
%show only right half of spectrum (symmetric)
axis ([0 (N-1)/Fs 0 Fs/2]);
axis xy;
xlabel('time[sec]');
ylabel('frequency[Hz]');
title(['Q.8c: STFT with win = ', num2str(win), ', hopsize = ', num2str(hopSize), ', F = ', num2str(F)]);

% figure;
% spectrogram(touchtone_1,1024,[],256,'yaxis')

%% Q.8d

win = 2048;
hopSize = 2048;
F = 256;
S = STFT(touchtone_1, win, hopSize, F, Fs);

% plotting the spectrogram

S_size = size(S);
T_vec = linspace(0,(N-1)/Fs, S_size(2));
F_vec = linspace(-Fs/2, Fs/2, S_size(1));

figure;
imagesc(T_vec, F_vec, abs(S));
%show only right half of spectrum (symmetric)
axis ([0 (N-1)/Fs 0 Fs/2]);
axis xy;
xlabel('time[sec]');
ylabel('frequency[Hz]');
title(['Q.8d: STFT with win = ', num2str(win), ', hopsize = ', num2str(hopSize), ', F = ', num2str(F)]);

% figure;
% spectrogram(touchtone_1,2048,[],256,'yaxis')

%% Q.10

digits_1 = decode(touchtone_1);

%% Q.11

[touchtone_2, Fs] = audioread('touchtone_2.wav');
load('touchtone_2_sequence');
sequence = double(sequence);

N = length(touchtone_2);
sigma = [0.05, 0.1, 0.25, 0.5, 1, 2];
error_rate = zeros(1, length(sigma));

for i=1:length(sigma)
    signal = touchtone_2 + sigma(i)*randn(N,1);
    win = 256;
    hopSize = 256;
    F = 128;
    S = STFT(signal, win, hopSize, F, Fs);

    % plotting the spectrogram
    S_size = size(S);
    T_vec = linspace(0,(N-1)/Fs, S_size(2));
    F_vec = linspace(-Fs/2, Fs/2, S_size(1));

    figure;
    imagesc(T_vec, F_vec, abs(S));
    %show only right half of spectrum (symmetric)
    axis ([0 (N-1)/Fs 0 Fs/2]);
    axis xy;
    xlabel('time[sec]');
    ylabel('frequency[Hz]');
    title(['STFT of touchtone_2 with AGW with sigma = ', num2str(sigma(i)), ', hopsize = ', num2str(hopSize), ', F = ', num2str(F)]);
end

for i = 1:length(sigma)
    signal = touchtone_2 + sigma(i)*randn(N,1);
    digits_2 = decode(signal);
    error_rate(i) = sum((digits_2~=sequence))/length(sequence);
end

figure;
plot(sigma,error_rate); axis tight;
title('Q.11: error rate as a function of sigma');
xlabel('sigma');ylabel('error rate');

%% Q.12

[touchtone_3, Fs] = audioread('touchtone_3.wav');

digits_3 = decode(touchtone_3);
