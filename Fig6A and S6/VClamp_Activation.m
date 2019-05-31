
%%%%%% passive parameters:
d = 40; % diameter is microns
gh_max = 0.027; % maximal h-conductance (mS/cm2) at -120 mV
d = d/10000; % diameter in centimeters
area = pi*(d^2); % Sphere A = 4pi*r^2 = pi*d^2 (cm^2)

%%%%%% Ih
gh_max = gh_max*area; % maximal h-conductance (mS)
Eh = -33.7; % reversal potential for Ih (mV)
A = 0.96; % voltage dependent component of activation curve
Vh = -90.7; % voltage at half activation (mV), i.e., V1/2
k = 12.5; % slope factor;
% fraction of fast activation
Vh1 = -101.1;
k1 = 9.701;
Vh2 = -65;
k2 = 0.3813;
B = 0.4499; % vertical translation term

%%%%% times and stimulation parameters
ts = 200; % start pulse time (ms)
tf = 1400; % terminate pulse time (ms)
tmax = tf + 200; % total time (ms)
Vpulse = [-60 -70 -80 -90 -100 -110 -120]; % pulse potentials (mV)
Npulse = length(Vpulse); % number of voltage steps
Vhold = -50; % holding potential

%%%%% sampling parameters, data storage vectors
Fs = 10000; % sampling rate (Hz)
dt = (1/Fs)*1000; % time step (ms)
time = 0:dt:tmax; % time vector (ms)
Ih = zeros(Npulse,length(time)); % h-current


%% initial conditions
Vm = Vhold; % starting voltage clamp value
Xinf = A/(1+exp((Vm-Vh)/k)) + (1-A); % initialize steady state conductance
Xf = Xinf; % set fast gating variable to steady state value
Xs = Xinf; % set slow gating variable to steady state value
F = ((0.6376-B)*(1+exp((Vm-Vh1)/k1))^(-1))+((0.6233-B)*(1+exp(-(Vm-Vh2)/k2))^(-1))+B; % fractional contribution of fast component
X = Xf*F + Xs*(1-F); % initialize gate variable

%% run the simulation (forward Euler)
for j = 1:Npulse
    for i = 2:length(time)
        
        t = time(i); % time (ms)
        
        % square step
        if (t >= ts) && (t <= tf)
            Vm = Vpulse(j);
        else
            Vm = Vhold;
        end
        
        %%%% Ih parameters and computation
        Taf = 129.5/(12.93*exp(Vm/22.09) + 0.2166*exp(-Vm/40.07));   % fast activation time constant    
        Tas = 122.1/(1.955*exp(Vm/22.45) + 0.01528*exp(-Vm/34.69));  % slow activation time constant    
        Tdf = 0.3843*Vm + 47.34;                                     % fast deactivation time constant    
        Tds = 30/(320.2*exp(Vm/7.243) + 0.05197*exp(-Vm/63.85));     % slow deactivation time constant    
        F = ((0.6376-B)*(1+exp((Vm-Vh1)/k1))^(-1))+((0.6233-B)*(1+exp(-(Vm-Vh2)/k2))^(-1))+B; % fractional contribution of fast component
        Xinf = A/(1+exp((Vm-Vh)/k)) + (1-A);                         % steady state conductance              
        
        if X <= Xinf %%% if X(t-dt) < Xinf(t)
            Xf = Xf + (dt/Taf)*(Xinf - Xf);
            Xs = Xs + (dt/Tas)*(Xinf - Xs);
            X = Xf*F + Xs*(1 - F); %%% X(t)
            Ih(j,i) = gh_max*X*(Vm - Eh);            
        else
            Xf = Xf + (dt/Tdf)*(Xinf - Xf);
            Xs = Xs + (dt/Tds)*(Xinf - Xs);
            X = Xf*F + Xs*(1 - F);
            Ih(j,i) = gh_max*X*(Vm - Eh); % Ih in uA       
        end                       
    end
end

figure;
plot(time,1000000*Ih(1:Npulse,:),'k','Linewidth',3);
axis([180 1420 -120 1]);title('Simulated Voltage-Clamp (Fig 6A and Fig S6)');xlabel('Time (ms)');ylabel('Ih (pA)');

