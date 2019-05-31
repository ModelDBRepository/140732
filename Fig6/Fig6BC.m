%% figure 6...the first model fig
% This code produces panels B and C in seven separate graphs

% This m-file calls Zsim1.m to get the voltage trace for the standard model and the 
% trace when Ih is absent (red and black traces in Fig 6C). Zsim2.m and Zsim3.m get 
% the voltage traces for the blue and green plots from Fig 6C, respectively

% Zprofile.m is called to generate Z-curves, raw and smoothed

gh_max = 0.027; % maximal h-conductance (mS/cm2)
Fs = 10000; % sampling rate (Hz)

%% Panels B: ZAP Voltage traces
%%%%% times and ZAP stimulation parameters
L = 15.0; % length of ZAP stimulus (sec)
ts = 500/1000; % start pulse time (sec)
tf = ts + L; % terminate pulse time (sec)
tmax = tf + 800/1000; % total time (s) [second term is the extra time]
Ipulse = 0.00002; % pulse magnitude (uA), 0.00002 uA = 20 pA
Fi = 0; % start freq
Ff = 15; % end freq
a = 2*pi*(Ff - Fi)/(2*L);
b = 2*pi*Fi;
% swept sine wave: Ipulse*sin(a*t.^2 + b*t) 
tzap = 0:(1/Fs):tmax; % time vector for the zap stimulus in seconds
Izap = zeros(1,length(tzap)); % space of the Izap current input
ii = find(tzap == ts); % index for start time of Izap
ff = find(tzap == tf)-1; % index for termination of Izap
Izap(ii:ff) = Ipulse*sin(a*(tzap(ii:ff)-tzap(ii)).^2 + b*(tzap(ii:ff)-tzap(ii))); % ZAP stimulus

%%%%% get voltage responses for each model variant
[time,Vm1,Vzd] = Zsim1 (gh_max,Izap,tzap,Fs); % top two traces: 1st - no Ih, 2nd - double exp Ih
Vm2 = Zsim2 (gh_max,Izap,tzap,Fs); % third trace, single exp, independent act/deact kinetics
Vm3 = Zsim3 (gh_max,Izap,tzap,Fs); % fourth trace, single exp, 1 fn for act/deact kinetics

%%%%% plot voltage resonses
figure; plot(time,Vzd,'r','Linewidth',2); axis([0.4 15.7 -80 -59]); title('No Ih'); xlabel('Time (sec)'); ylabel('Vm (mV)');
figure; plot(time,Vm1,'k','Linewidth',2); axis([0.4 15.7 -80 -59]); title('Standard Model'); xlabel('Time (sec)'); ylabel('Vm (mV)');
figure; plot(time,Vm2,'b','Linewidth',2); axis([0.4 15.7 -80 -59]); title('Single exp, 2 functions');xlabel('Time (sec)'); ylabel('Vm (mV)');
figure; plot(time,Vm3,'g','Linewidth',2); axis([0.4 15.7 -80 -59]); title('Single exp, 1 function');xlabel('Time (sec)'); ylabel('Vm (mV)');

%% Panels C: ZAP analysis (impedance and phase)

ZAPin_ZD = [time' Vzd'];
ZAPin_Vm1 = [time' Vm1'];
ZAPin_Vm2 = [time' Vm2'];
ZAPin_Vm3 = [time' Vm3'];

%%% compute Z = FFT(V)/FFT(I)
[ZZzd,f,freq,Zzd] = Zprofile(ZAPin_ZD,Fs,Izap);
[ZZ1,f,freq,Z1] = Zprofile(ZAPin_Vm1,Fs,Izap);
[ZZ2,f,freq,Z2] = Zprofile(ZAPin_Vm2,Fs,Izap);
[ZZ3,f,freq,Z3] = Zprofile(ZAPin_Vm3,Fs,Izap);

%%%%%%%%%%%%%%%%%%%% impedance plots
figure; hold on
fnplt(Zzd,'r');
fnplt(Z1,'k');
fnplt(Z2,'b');
fnplt(Z3,'--g');
hold off
title('Impedance profiles'); ylabel('|Impedance| (MOhm)'); xlabel('Frequency (Hz)');
axis([0 15.5 180 510]);

%%%%%%%%%%%%%%%%%%% phase plots
R = real(ZZzd); X = imag(ZZzd);
phizd = atand(X./R); % No Ih

R = real(ZZ1); X = imag(ZZ1);
phi1 = atand(X./R); % standard model 

R = real(ZZ2); X = imag(ZZ2);
phi2 = atand(X./R); % single exp, 2 fn

R = real(ZZ3); X = imag(ZZ3);
phi3 = atand(X./R); % single exp, 1 fn

figure; hold on
plot(f(freq),phizd(freq),'r','linewidth',2);
plot(f(freq),phi1(freq),'k','Linewidth',2);
plot(f(freq),phi2(freq),'b','Linewidth',2);
plot(f(freq),phi3(freq),'--g','Linewidth',2);
hold off
title('Phase Plots'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)');
axis([0 15.5 -70 15]);

