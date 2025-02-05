% Load the audio file
[audio_in, fs] = audioread('song2.wav');

% Convert to mono if stereo
if size(audio_in, 2) > 1
    audio_in = mean(audio_in, 2);
end

% Remove DC offset
audio_in = audio_in - mean(audio_in);

% Length of the audio signal
N = length(audio_in);

%% 1. Discrete-Time Signal Generation and Operations
% Generate a discrete-time sinusoidal signal in chunks
chunk_size = fs * 10; % Process 10-second chunks
sin_signal = zeros(size(audio_in)); % Preallocate

f_sin = 440; % Frequency in Hz
for start_idx = 1:chunk_size:N
    end_idx = min(start_idx + chunk_size - 1, N);
    n = start_idx:end_idx; % Indices for this chunk
    sin_signal(start_idx:end_idx) = 0.1 * sin(2 * pi * f_sin * n / fs);
end

% Add the sinusoidal signal to the audio
audio_with_sin = audio_in + sin_signal;

%% 2. Discrete-Time System
% Simple low-pass filter (moving average)
M = 5; % Number of points
b = (1/M) * ones(1, M); % Filter coefficients
audio_filtered = filter(b, 1, audio_with_sin);

%% 3. Frequency Domain Analysis
% Apply windowing to reduce artifacts in FFT
window = hann(N); % Hanning window
FFT_audio = fft(audio_filtered .* window);
freq_axis = (0:N-1) * (fs / N);
freq_axis = freq_axis - fs * (freq_axis > fs/2);

% Define frequency thresholds
low_thresh = 300;  % Hz for instrumentals
high_thresh = 5000; % Hz for vocals

% Create masks for filtering
instrumental_mask = (abs(freq_axis) < low_thresh);
vocal_mask = (abs(freq_axis) >= low_thresh) & (abs(freq_axis) <= high_thresh);

% Apply masks directly to the FFT components
FFT_instrumental = FFT_audio;
FFT_instrumental(~instrumental_mask) = 0;

FFT_vocal = FFT_audio;
FFT_vocal(~vocal_mask) = 0;

% Inverse FFT to get time-domain signals
instrumental = real(ifft(FFT_instrumental));
vocal = real(ifft(FFT_vocal));

% Normalize the signals to avoid clipping
instrumental = instrumental / max(abs(instrumental));
vocal = vocal / max(abs(vocal));

%% 4. FIR and IIR Filter Combination
% Design a high-pass FIR filter for vocal extraction
hp_cutoff = 400; % Hz
normalized_cutoff = hp_cutoff / (fs / 2);
fir_filter = fir1(50, normalized_cutoff, 'high');
fir_vocal = filter(fir_filter, 1, vocal);

% Design a band-pass IIR filter for vocal enhancement
bp_cutoff = [400, 3500]; % Hz
normalized_bp = bp_cutoff / (fs / 2);
[iir_b, iir_a] = butter(4, normalized_bp, 'bandpass');

% Apply FIR followed by IIR
fir_iir_vocal = filter(iir_b, iir_a, fir_vocal);

%% 5. Spectrogram Generation
figure;
subplot(2,1,1);
spectrogram(instrumental, 1024, [], [], fs, 'yaxis');
title('Spectrogram of Instrumental');
colorbar;

subplot(2,1,2);
spectrogram(fir_iir_vocal, 1024, [], [], fs, 'yaxis');
title('Spectrogram of FIR + IIR Filtered Vocal');
colorbar;

%% Save Outputs
audiowrite('instrumental.wav', instrumental, fs);
audiowrite('vocal_fir_iir.wav', fir_iir_vocal, fs);

%% Play Processed Audio
%disp('Playing the instrumental...');
%sound(instrumental, fs);
%pause(5); % Wait for instrumental playback to finish

%disp('Playing the combined FIR and IIR filtered vocal%...');
%sound(fir_iir_vocal, fs);

%disp('Separation completed and saved as instrumental.wav and vocal_fir_iir.wav.');%
