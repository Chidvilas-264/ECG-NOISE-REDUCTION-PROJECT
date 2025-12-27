clc;
clear;
close all;
% AN OPTIMIZED FIR FILTER DESIGN FOR NOISE REMOVAL IN ECG SIGNAL
fs = 125; % Sampling frequency (in Hz)
Ts = 1/fs; % Sampling period (in seconds)
duration = 60; % Signal duration (in seconds)
% Loading the ECG signal from Datasets
x = load("data_3.mat"); 
y = x.val/2000; % Scaling the signal
y = y'; % Transpose the signal to ensure the correct format
% Number of samples in the signal
N = length(y); % Length of the ECG signal from the dataset
% Construct the continuous time vector for 1 minute
t_continuous = linspace(0,duration,N); % Continuous time vector corresponding to the signal
% Ensure that both t_continuous and y have the same length
assert(length(t_continuous) == length(y),"Time vector and ECG signal must have same length");
% Plot the ECG signal with its continuous time vector for entire 1 minute
figure;
plot(t_continuous,y(:,1));
title("ECG signal in continuous over 1 minute");
xlabel("Time(s)");
ylabel("Amplitude");
grid on;
% Applying sampling theorem to convert continuous time signal to discrete time signal
t_discrete = 0:Ts:max(t_continuous); % Time vector for discrete time signal
% Sample the discrete time ECG signal using interpolation
ecg_signal_discrete = interp1(t_continuous,y,t_discrete,"linear");
% Plot the original continuous and discrete sampled ECG signal for entire 1 minute
figure;
plot(t_continuous,y(:,1),"b","LineWidth",2,"DisplayName","Continuous ECG");
hold on;
stem(t_discrete,ecg_signal_discrete,"r","MarkerSize",2,"DisplayName","Discrete ECG");
title("Continuous vs Discrete sampled ECG signal for 1 minute");
xlabel("Time (s)");
ylabel("Amplitude");
legend;
grid on;
xlim([0,60]);
% Applying Z-transform to convert from time domain to frequency domain
syms z n;
ecg_z_transform = sum(ecg_signal_discrete.*z.^(-n));
% Designing a bandpass FIR filter for noise removal
f_pass = [0.67 50]; % Passband frequencies in Hz
f_nyquist = fs/2;
f_pass_norm = f_pass/f_nyquist; % Normalized passband
order = 150; % Filter order
b = fir1(order,f_pass_norm,"bandpass"); % Design of bandpass filter
% Filtering the ECG signal
filtered_ecg_signal = filter(b,1,ecg_signal_discrete);
% Plotting the filtered ECG signal
figure;
plot(t_discrete,filtered_ecg_signal);
title("Filtered ECG signal (Noise Removed)");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
% Overlaying original and filtered signals for comparison
figure;
plot(t_discrete,ecg_signal_discrete,"b","DisplayName","Original ECG");
hold on;
plot(t_discrete,filtered_ecg_signal,"r","DisplayName","Filtered ECG");
title("Original vs Filtered ECG signal");
xlabel("Time (s)");
ylabel("Amplitude");
legend;
grid on;
% Detecting P-wave, QRS-complex, T-wave Mathematically
% P-wave detection
[p_wave_peaks,p_wave_locs] = findpeaks(filtered_ecg_signal(:,1),"MinPeakHeight",0.15,"MinPeakDistance",0.2*fs);
% QRS complex detection
[qrs_peaks,qrs_locs] = findpeaks(filtered_ecg_signal(:,1),"MinPeakHeight",0.6,"MinPeakDistance",0.2*fs);
% T-wave detection
[t_wave_peaks,t_wave_locs] = findpeaks(filtered_ecg_signal(:,1),"MinPeakHeight",0.3,"MinPeakDistance",0.4*fs);
% Plotting detected features on the filtered ECG signal
figure;
plot(t_discrete,filtered_ecg_signal);
hold on;
plot(p_wave_locs/fs,p_wave_peaks,"rv","MarkerFaceColor","r","DisplayName","P-wave");
plot(qrs_locs/fs,qrs_peaks,"gs","MarkerFaceColor","g","DisplayName","QRS-complex");
plot(t_wave_locs/fs,t_wave_peaks,"mo","MarkerFaceColor","m","DisplayName","T-wave");
title("Detected P-wave, QRS-complex, and T-wave on Filtered ECG Signal");
xlabel("Time (s)");
ylabel("Amplitude");
legend;
grid on;
% Adjusting plot parameters for clarity
xlim([0,60]);
ylim([-1,1]);
rr_intervals = diff(qrs_locs)/fs; % R-R intervals (in seconds)
heart_rate = 60./rr_intervals; % Heart rate in BPM
if mean(heart_rate) < 60
    disp("Bradycardia : Heart rate is too slow");
elseif mean(heart_rate) > 100
    disp("Tachycardia : Heart rate is too fast");
else
    disp("Normal Heart rate");
end
% Calcualting the frequenies of p-wave, qrs-complex and t-wave
p_wave_freq = 1 ./ mean(diff(p_wave_locs) / fs);
disp("frequency of p-wave is :");disp(p_wave_freq);
qrs_freq = 1 ./ mean(diff(qrs_locs) / fs);
disp("frequency of qrs-complex is :");disp(qrs_freq);
t_wave_freq = 1 ./ mean(diff(t_wave_locs) / fs);
disp("frequency of t-wave is :");disp(t_wave_freq);
%FOR HEALTHY PERSON THE FREQUENCY VALUES ARE GIVEN AS : 
% P-WAVE : 0.5 TO 10 Hz
% QRS-COMPLEX : 10 TO 100 Hz
% T-WAVE : 0.5 TO 10 Hz
if (p_wave_freq > 0.5 && p_wave_freq < 10) && (qrs_freq > 10 && qrs_freq < 100) && (t_wave_freq > 0.5 && t_wave_freq < 10)
    disp("The patient is healthy.");
else
    disp("The patient is not healthy.");
end
