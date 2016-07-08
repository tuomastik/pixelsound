% This script converts image into sound and visualizes the produced sound
% with a spectrogram, where the original image is seen.

% The conversion is done by converting each pixel in the image to a sine
% wave, where frequency corresponds to the row of the pixel and amplitude
% corresponds to the intensity of the pixel.

clear all
close all
clc

%% Settings

path_to_im = fullfile('..', 'data', 'im.jpg');
max_freq = 20500; % Hz
fs = max_freq * 2; % Sampling frequency (Nyquist–Shannon theorem)

%% Open image

im = imread(path_to_im); % Open image
if length(size(im)) == 3
    im = rgb2gray(im); % Convert to grayscale if not already
end
im = double(im);
% im = imresize(im, 0.5);
im = mat2gray(im); % Normalize to [0...1]
figure(1), imshow(im, []), title('Original image');
[rows, cols] = size(im);

%% Convert each pixel into sine wave

col_dur_sec = rows/max_freq; % Duration of column in seconds
col_dur_sam = ceil(col_dur_sec*fs); % Duration of column in samples
t = 0:1/fs:col_dur_sec; % Time steps to sample sine wave
y = zeros(1, cols*col_dur_sam+1); % Initialize output audio signal
freqs = max_freq-(1:rows)/rows*max_freq; % Precalculate frequencies

for col = 1:cols
    
    start_sample = (col - 1) * col_dur_sam + 1;
    end_sample = start_sample + col_dur_sam;
    
    for row = 1 : rows
        
        amplitude = im(row, col)^2; % Pixel intensity
        
        y(1, start_sample:end_sample) = y(1, start_sample:end_sample) + ...
                                        amplitude * sin(2*pi*freqs(row)*t);
        
    end
end

% y = y./abs(max(y)); % Normalize 
soundsc(y, fs); % Listen
% clear sound; % Stop playback
n_bits = 32;
out_file = fullfile('..', 'data', 'im.wav');
audiowrite(out_file, y, fs, 'BitsPerSample', n_bits) % Save audio

%% Draw spectrograms with different window functions

windows = {'barthannwin', 'bartlett', 'blackman', 'blackmanharris', ...
    'chebwin', 'flattopwin', 'gausswin', 'hamming', 'hann', ...
    'kaiser', 'nuttallwin', 'parzenwin', 'rectwin', ...
    'tukeywin', 'taylorwin', 'triang'};

for i = 1 : length(windows)
    
    window = str2func(windows{i});
    window = window(col_dur_sam);
    overlapping_samples = 0;
    n_dft_points = col_dur_sam;
    [~, F, T, P] = spectrogram(y, window, overlapping_samples, ...
                               n_dft_points, fs);
    figure(2); subplot(4, 4, i);
    surf(T, F, 10*log10(abs(P)), 'edgecolor', 'none');
    title(windows{i});
    colormap(gray); axis tight;
    set(gca,'YTickLabel',[]), set(gca,'XTickLabel',[]);
    view(0, -270)

end

%% Draw spectrogram with single window function (MATLAB implementation)

window_name = 'kaiser';
window = str2func(window_name);
window = window(col_dur_sam);
overlapping_samples = 0;
n_dft_points = col_dur_sam;
[~, F, T, P] = spectrogram(y, window, overlapping_samples, ...
                           n_dft_points, fs);
figure(3);
surf(T, F, 10*log10(abs(P)), 'edgecolor', 'none');
title('Spectrogram');
colormap(gray); axis tight;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
ylim([0 max_freq]);
xlim([0 col_dur_sec*cols]);
view(0, -270);
    
%% Draw spectrogram with single window function (custom implementation)

chunk = col_dur_sam;
overlap = col_dur_sam;
nfft = col_dur_sam;
num = chunk : overlap : length(y);
spectrog = zeros(length(num), nfft/2);

for i = length(num):-1:1,
    framed = y( (num(i)-chunk) + 1 : num(i) );
    windowed = framed .* kaiser(col_dur_sam)';
    fourrier = fft(windowed, nfft);
    spectrog(i, :) = log(abs(fourrier(1:nfft/2)));
end

figure(4); imshow(flipud(spectrog'), [])
title('Spectrogram');
