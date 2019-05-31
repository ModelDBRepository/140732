function [impedance,f,freq,pp] = Zprofile (ZAPin,Fs,Izap)
zap = ZAPin;

% time vector and frequency range
T = 1/Fs;   % time step (seconds)                                       
len = length(zap(:,1)); % column length (number of data points in a trace)
N = 2^(nextpow2(len));  % Number of points used for calculation of FFT
range = N/2;            % Range for Z(f) plot
f = Fs*(0:range-1)/N;   % Frequency range - up to ~ half of Fs (i.e., Fs/2)
d = Fs/N;               % increment in frequency vector (unit-Hz)
fi = 14*d; ff = 393*d;  % start and end frequencies for peak analysis (0.5 to 15 Hz)
freq = fi/d:ff/d;       % frequency range for plots

%% compute FFT of voltage trace and return the FFT 
% get average zap response about 0.0 mV
T_ms = T*1000; % time step in ms
col = length(zap(1,:)); % number of columns (row length)
x = zap(:,2:col);    % vector for computing avg trace
xbar = mean(x,2); % trace of averaged data
y = 100/T_ms; % number of samples in 100 ms
one = ones(1,y);
len1 = length(one);
subtr = one*xbar(1:y);
subtract = subtr/len1;  % amount that needs to be subtracted for Vm to oscillate about zero
                        % based on the average Vm taken for the first 100
                        % ms of the sample, which is before the zap current injection
zap0 = xbar - subtract; % Vm response to zap, plotted about 0.0 mV
FFTv = fft(zap0,N);     % compute FFT(V)

%% compute FFT of current input and return the FFT
FFTi =  fft(Izap',N);

%% calculate impedance(f), resonance freqency, and Q-factor
impedance = (1/1000)*(FFTv./FFTi); % complex impedance = Zreal + jZimaginary (MOhm)
Zf = abs(impedance); % magnitude of impedance vector = |Z|(f) in MOhm
[pp,p] = csaps(f(freq),Zf(freq),0.7); % cubic spline fit
fit = fnplt(pp); % row1: freq, row2: Z
index = find(fit(2,:) == max(fit(2,:))); % freq index for Z = Zmax
Fres = fit(1,index); % resonant frequency
if length(Fres) > 1
    Fres = mean(Fres);
end
Fres; % resonant frequency
Q = fit(2,index(1))/fit(2,1); % Q factor

