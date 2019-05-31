%% figure 7
% This code produces panels A, C, and D in three separate graphs

% 1st parameter set:
Vh = -90.71;  % V1/2 for Ih activation
gEsyn = -0.0022; % synaptic conductance
k = 12.51;       % slope factor;

%%%%%%%%%%%%%%% parameters and initial conditions for the model
Vrest = -70; % starting potential (mV)
d = 40; % diameter is microns
Cm = 1.0; % specific capacitance (uF/cm2)
gL = 0.04; % leak conductance (mS/cm2)
gh_max = 0.1; % maximal h-conductance (mS/cm2) at -120 mV
d = d/10000; % diameter in centimeters
area = pi*(d^2); % Sphere A = 4pi*r^2 = pi*d^2 (cm^2)
Cm = Cm*area; % total capacitance (uF)
gL = gL*area; % total leak conduectance (mS)
EL = -75; % leak reversal potential (mV)
I1 = -0.0504; % DC for Vm = -70
I2 = 0.02; % DC for Vm = -70

%%%%%% Ih
gh_max = gh_max*area; % maximal h-conductance (mS)
Eh = -33.7; % reversal potential for Ih (mV)
A = 0.96; % voltage dependent component of activation curve
% fraction of fast activation:
Vh1 = -101.1;
k1 = 9.701;
Vh2 = -65;
k2 = 0.3813;
B = 0.4499; % vertical translation term

% sampling parameters, data storage vectors
tmax = 1000; % total time (ms)
Fs = 10000; % sampling rate (Hz)
dt = (1/Fs)*1000; % time step (ms)
time = 0:dt:tmax; % time vector (ms)
V_Iht = zeros(1,length(time)); % voltage vector
V_ZD = zeros(1,length(time)); % second voltage vector
V_Iht0 = zeros(1,length(time)); % second voltage vector
Xinf = zeros(1,length(time)); % steady state activation level
Xf = zeros(1,length(time)); % fast activation gate for Ih
Xs = zeros(1,length(time)); % slow activation gate for Ih
X = zeros(1,length(time)); % gate variable
Ih = zeros(1,length(time)); % h-current
gh = zeros(1,length(time)); % h-conductance
Ih_t0 = zeros(1,length(time)); % h-current

%% synaptic input parameters
spk = zeros(2,length(time)); % matrix for spike train 
spk(2,:) = time;
Npulse = 5; % number of pulses
freq_pulse = 50; % frequency of input train (Hz)
ISI = (1/freq_pulse)*1000; % interstimulus interval (ms)
Ts = 200; % time (in ms) for the first pulse
tspk = Ts:ISI:Ts+(Npulse-1)*ISI; % vector of spike times
% synaptic conductance waveform
trise = 1; % rise time constant
tfall = 7; % decay time constant
Esyn = 0; % synaptic reversal potential
g = exp(-time./tfall) - exp(-time./trise); % waveform
tp = (trise*tfall/(trise - tfall))*log(trise/tfall);
gwave = g./(exp(-tp/tfall) - exp(-tp/trise)); % scaled waveform
%%% fill the first row of spk with ones occuring at spike times
for i = 1:length(tspk)
    ind = find(spk(2,:) == tspk(i));
    spk(1,ind) = 1;
end
gsyn = conv(gwave,spk(1,:)); % convolve spike train with gsyn
gin = gEsyn.*gsyn(1:length(time)); % clip extra values at end, this is the sanaptic conductance input

%% initial conditions
V_Iht(1) = Vrest; % first membrane potential trace
V_ZD(1) = Vrest; % sencond membrane potential trace
V_Iht0(1) = Vrest; % sencond membrane potential trace
Xinf(1) = A/(1+exp((V_Iht(1)-Vh)/k)) + (1-A); % initialize steady state conductance
Xf(1) = Xinf(1); % set fast gating variable to steady state value
Xs(1) = Xinf(1); % set slow gating variable to steady state value
F = ((0.6376-B)*(1+exp((V_Iht(1)-Vh1)/k1))^(-1))+((0.6233-B)*(1+exp(-(V_Iht(1)-Vh2)/k2))^(-1))+B; % fractional contribution of fast component
X(1) = Xf(1)*F + Xs(1)*(1-F); % initialize gate variable

%% run the simulation (forward Euler)
for i = 2:length(time)
    
    %%%% Ih(t) computation
    Taf = 129.5/(12.93*exp(V_Iht(i-1)/22.09) + 0.2166*exp(-V_Iht(i-1)/40.07));   % fast activation time constant    
    Tas = 122.1/(1.955*exp(V_Iht(i-1)/22.45) + 0.01528*exp(-V_Iht(i-1)/34.69));  % slow activation time constant    
    Tdf = 0.3843*V_Iht(i-1) + 47.34;                                     % fast deactivation time constant    
    Tds = 30/(320.2*exp(V_Iht(i-1)/7.243) + 0.05197*exp(-V_Iht(i-1)/63.85));     % slow deactivation time constant    
    F = ((0.6376-B)*(1+exp((V_Iht(i-1)-Vh1)/k1))^(-1))+((0.6233-B)*(1+exp(-(V_Iht(i-1)-Vh2)/k2))^(-1))+B; % fractional contribution of fast component
    Xinf(i) = A/(1+exp((V_Iht(i-1)-Vh)/k)) + (1-A);                         % steady state conductance        
    if X(i-1) <= Xinf(i) %%% if X(t-dt) < Xinf(t)
        Xf(i) = Xf(i-1) + (dt/Taf)*(Xinf(i) - Xf(i-1));
        Xs(i) = Xs(i-1) + (dt/Tas)*(Xinf(i) - Xs(i-1));
        X(i) = Xf(i)*F + Xs(i)*(1 - F); %%% X(t)
        gh(i) = gh_max*X(i);
        Ih(i) = gh(i)*(V_Iht(i-1) - Eh);            
    else
        Xf(i) = Xf(i-1) + (dt/Tdf)*(Xinf(i) - Xf(i-1));
        Xs(i) = Xs(i-1) + (dt/Tds)*(Xinf(i) - Xs(i-1));
        X(i) = Xf(i)*F + Xs(i)*(1 - F); %%% X(t)
        gh(i) = gh_max*X(i);
        Ih(i) = gh(i)*(V_Iht(i-1) - Eh);            
    end 
    
    %%%% Ih(t0) computation
    Ih_t0(i) = gh(2)*(V_Iht0(i-1) - Eh);    
    
    %%%% compute voltage: Cm(dVm/dt) = -(Ih + Ileak) + Iappllied %%%%
    V_Iht(i) = V_Iht(i-1) - (dt/Cm)*(Ih(i) + gL*(V_Iht(i-1) - EL)) + gin(i)*(V_Iht(i-1) - Esyn) + I1;
    V_ZD(i) = V_ZD(i-1) - (dt/Cm)*(gL*(V_ZD(i-1) - EL)) + gin(i)*(V_ZD(i-1) - Esyn) + I2;
    V_Iht0(i) = V_Iht0(i-1) - (dt/Cm)*(Ih_t0(i) + gL*(V_Iht0(i-1) - EL)) + gin(i)*(V_Iht0(i-1) - Esyn) + I1;  
end


%% Summation Indices, SI = 100*(EPSP5 - EPSP1)/EPSP1 
ti = Ts; % first pulse time
tf = ti + (Npulse-1)*ISI; % last pulse time
i_i = find(time == ti); % first pulse time index
i_f = find(time == tf); % last pulse time index
in = Fs*ISI/1000; % number of samples between each pulse

%%% zd SI
EPSP1 = max(V_ZD(i_i:i_i+in)) - Vrest;
EPSP5 = max(V_ZD(i_f:i_f+in)) - Vrest;
SI_zd = (EPSP5 - EPSP1)/EPSP1;

%%% Ih(t) SI
EPSP1 = max(V_Iht(i_i:i_i+in)) - Vrest;
EPSP5 = max(V_Iht(i_f:i_f+in)) - Vrest;
SI_Ih = (EPSP5 - EPSP1)/EPSP1;

%%% Ih(t0) SI
EPSP1 = max(V_Iht0(i_i:i_i+in)) - Vrest;
EPSP5 = max(V_Iht0(i_f:i_f+in)) - Vrest;
SI_Ih2 = (EPSP5 - EPSP1)/EPSP1;
     
%%%% percent effect of gh(t0) = SI_D:
d1 = SI_zd - SI_Ih;
d2 = SI_zd - SI_Ih2;
SI_D = (1 - d2/d1)*100 

%% plots
%%%%% Fig 7A
figure; hold on
plot(time,V_ZD,'r','Linewidth',2);
plot(time,V_Iht0,'Color',[0.6 0.6 0.6],'Linewidth',2);
plot(time,V_Iht,'k','Linewidth',2);
hold off
axis([190 450 -72 -54]);title('Fig 7A');xlabel('Time (ms)');ylabel('V_m (mV)');

%%%%% Fig 7aC
figure; hold on
plot(time,1000000*Ih,'k','Linewidth',1.5);
plot(time,1000000*Ih_t0,'Color',[0.6 0.6 0.6],'Linewidth',1.5);
hold off
axis([190 450 -40 -18]);title('Fig 7Ca');xlabel('Time (ms)');ylabel('I_h (pA)'); 

%%%%% Fig 8C
figure; hold on
plot(time,(10^9)*gh,'k','Linewidth',2);
plot(time,(10^9)*gh(2),'Color',[0.6 0.6 0.6],'Linewidth',2);
hold off
axis([190 450 650 1050 ]);title('Fig 7Cb');xlabel('Time (ms)');ylabel('g_h (pS)'); 
legend('g_h(t)','g_h(t0)');



