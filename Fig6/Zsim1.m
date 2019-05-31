function [time,Vm,V2] = Zsim1 (gh_max,Izap,tzap,Fs)

% basic modeling with ZAP input:
%%%%% this is an updated version of the model where one function is used
%%%%% for the fraction of the fast current component (F)
%%%%% Based on VCAct2.m and ZAP3.m

%%%%%% passive parameters:
d = 40; % diameter is microns
d = d/10000; % diameter in centimeters
area = pi*(d^2); % Sphere A = 4pi*r^2 = pi*d^2 (cm^2)
Cm = 1; % specific capacitance (uF/cm2)
Cm = Cm*area; % total capacitance (uF)
gL = 0.04; % leak conductance (mS/cm2)
gL = gL*area; % total leak conduectance (mS)
EL = -75; % leak reversal potential (mV)

%%%%%% Ih
gh_max = gh_max*area; % maximal h-conductance (mS)
Eh = -33.7; % reversal potential for Ih (mV)
A = 0.96; % voltage-dependent component of activation curve
Vh = -90.7; % voltage at half activation (mV), i.e., V1/2
k = 12.5; % slope factor;
% fraction of fast activation
Vh1 = -101.1;
k1 = 9.701;
Vh2 = -65;
k2 = 0.3813;
B = 0.4499; % vertical translation term

%%%%% sampling parameters, data storage vectors
dt = (1/Fs)*1000; % time step (ms)
time = tzap; % time vector (sec)
Ih = zeros(1,length(time)); % h-current
Vm = zeros(1,length(time)); % voltage vector for Vm1 - with Ih
V2 = zeros(1,length(time)); % second voltage vector for Vzd - no Ih
Xf = zeros(1,length(time)); % fast activation gate for Ih
Xs = zeros(1,length(time)); % slow activation gate for Ih

%% initial conditions
Vstart = -70; % initial voltage
I1 = 0.0005; % DC for cell with Ih
I2 = 0.0198; % DC for cell without Ih
Vm(1) = Vstart; % starting voltage value
V2(1) = Vstart; % starting voltage value
Xinf = A/(1+exp((Vm(1)-Vh)/k)) + (1-A); % initialize steady state conductance
Xf(1) = Xinf; % set fast gating variable to steady state value
Xs(1) = Xinf; % set slow gating variable to steady state value
F = ((0.6376-B)*(1+exp((Vm(1)-Vh1)/k1))^(-1))+((0.6233-B)*(1+exp(-(Vm(1)-Vh2)/k2))^(-1))+B; % fractional contribution of fast component
X = Xf*F + Xs*(1-F); % initialize gate variable

%% run the simulation (forward Euler)
 for i = 2:length(time)
     
     %%%% Ih parameters and computation
     Taf = 129.5/(12.93*exp(Vm(i-1)/22.09) + 0.2166*exp(-Vm(i-1)/40.07));   % fast activation time constant    
     Tas = 122.1/(1.955*exp(Vm(i-1)/22.45) + 0.01528*exp(-Vm(i-1)/34.69));  % slow activation time constant    
     Tdf = 0.3843*Vm(i-1) + 47.34;                                     % fast deactivation time constant    
     Tds = 30/(320.2*exp(Vm(i-1)/7.243) + 0.05197*exp(-Vm(i-1)/63.85));     % slow deactivation time constant    
     F = ((0.6376-B)*(1+exp((Vm(i-1)-Vh1)/k1))^(-1))+((0.6233-B)*(1+exp(-(Vm(i-1)-Vh2)/k2))^(-1))+B; % fractional contribution of fast component
     Xinf = A/(1+exp((Vm(i-1)-Vh)/k)) + (1-A);                         % steady state conductance              
        
     if X <= Xinf %%% if X(t-dt) < Xinf(t)
         Xf(i) = Xf(i-1) + (dt/Taf)*(Xinf - Xf(i-1));
         Xs(i) = Xs(i-1) + (dt/Tas)*(Xinf - Xs(i-1));
         X = Xf(i)*F + Xs(i)*(1 - F); %%% X(t)
         Ih(i) = gh_max*X*(Vm(i-1) - Eh);            
     else
         Xf(i) = Xf(i-1) + (dt/Tdf)*(Xinf - Xf(i-1));
         Xs(i) = Xs(i-1) + (dt/Tds)*(Xinf - Xs(i-1));
         X = Xf(i)*F + Xs(i)*(1 - F); %%% X(t)
         Ih(i) = gh_max*X*(Vm(i-1) - Eh); % Ih in uA            
     end 
        
     %%%% compute voltage: Cm(dVm/dt) = -(Ih + Ileak) + Iappllied %%%%
     Vm(i) = Vm(i-1) - (dt/Cm)*(Ih(i) + gL*(Vm(i-1) - EL)) + (dt/Cm)*Izap(i) + I1;  % for Vm1 - with Ih
     V2(i) = V2(i-1) - (dt/Cm)*(gL*(V2(i-1) -EL)) + (dt/Cm)*Izap(i) + I2; % for Vzd - no Ih
     
 end
 
end
