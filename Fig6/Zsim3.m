function [Vm] = Zsim3 (gh_max,Izap,tzap,Fs)

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
A = 0.96; % voltage dependent component of activation curve
Vh = -90.7; % voltage at half activation (mV), i.e., V1/2
k = 12.5; % slope factor;

%%%%% sampling parameters, data storage vectors
dt = (1/Fs)*1000; % time step (ms)
time = tzap; % time vector (sec)
Vm = zeros(1,length(time)); % voltage vector
Ih = zeros(1,length(time)); % h-current

%% initial conditions
Vstart = -70; % initial voltage
I1 = 0.0005; % DC for cell with Ih
Vm(1) = Vstart; % starting voltage value
Xinf = A/(1+exp((Vm(1)-Vh)/k)) + (1-A); % initialize steady state conductance
X = Xinf;

%% run the simulation (forward Euler)
 for i = 2:length(time)         
               
     %%%% Ih parameters and computation
     Th = 79.99/(15*exp(Vm(i-1)/20) + 0.002454*exp(-Vm(i-1)/20));   % fast activation time constant                    
     Xinf = A/(1+exp((Vm(i-1)-Vh)/k)) + (1-A);                        % steady state conductance          
     X = X + (dt/Th)*(Xinf - X);
     Ih(i) = gh_max*X*(Vm(i-1) - Eh);
        
        
     %%%% compute voltage: Cm(dVm/dt) = -(Ih + Ileak) + Iappllied %%%%
     Vm(i) = Vm(i-1) - (dt/Cm)*(Ih(i) + gL*(Vm(i-1) - EL)) + (dt/Cm)*Izap(i) + I1;
     
 end

end