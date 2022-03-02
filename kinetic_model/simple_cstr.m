% 15wt% MEA Cp ~ 3.4145 J/(g*K) @ mol loading 0.1, @ 25C. https://pubs.acs.org/doi/pdf/10.1021/je960314v
% CO2 1.9795E4+7.3437E1 *T + -5.6019E-2 *T^2 + 1.7153E-5 * T^3 from 300K-1088.6K @ Aspen Databank
% H_abs = -32000J/mol (CO2) `

function [sys,x0,str,tss]=simple_cstr(t,x,u,flag,Param,X_ss)

global d cp V CAin DH1 DH2 UA k01 k02 EA1 EA2 Tin F Tu CDin PA PB PC PD PHEAT y0 r X0;

%
% An s-function to model a simple CSTR with first order kinetics and the following
% reaction ( A + B --> C) with a variable height well mixed vessel.  The equations that describe
% ( CO2G + MEA -> Absoebed)
% the system are:
% dV/dt  = F_in - F
% R1     = k1*Cd*Ca
% R2     = k2*Ca
% dCa/dt = F_in/V*(Ca_in-Ca)-R1-R2
% dCb/dt = F_in/V*(Cb_in-Cb)+R1
% dCc/dt = F_in/V*(Cc_in-Cc)+R2
% dCd/dt = F_in/V*(Cd_in-Cd)-R1
% dT /dt = - T/V*(F_in - F) + 1/V *(F_in*T_in - F*T) + 1/d/Cp*(-DH1)*R1+ 1/d/Cp*(-DH2)*R2 + 1/V/d/Cp*U*A*(Tc - T)
%================================================================================

switch flag,

case 0,	% Initialize the states and sizes
   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);
   
    % ****************
  	%  Outputs
  	% ****************   
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param);
   
   % ****************
   % Update
   % ****************
case 1,	% Obtain derivatives of states
   
   sys = mdlDerivatives(t,x,u,Param);

otherwise,
   sys = [];
end

% ******************************************
% Sub-routines or Functions
% ******************************************

% ******************************************
% Initialization
% ******************************************
function [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);

global d cp V CAin DH1 DH2 UA k01 k02 EA1 EA2 Tin F Tu CDin PA PB PC PD PHEAT y0 r X0;

% This handles initialization of the function.
% Call simsize of a sizes structure.
sizes = simsizes;
sizes.NumContStates  = 6;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 7;     % outputs of model (manipulated variables)
sizes.NumInputs      = 8;     % inputs of model [r;Y;Ud]
sizes.DirFeedthrough = 1;     % System is not causal, there is explicit dependence of outputs on inputs.
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % Initialize the discrete states (inputs).
%x0 = [0.3, 0.2, 0, 0.15, 0, 298];
str = [];	                  % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

% ******************************************
%  Outputs
% ******************************************
function sys = mdlOutputs(t,x,u,Param);

global DH1 DH2 k01 k02 EA1 EA2 alpha1 alpha2 beta1 beta2;
global P_bot P_top Vliq Vvap rou Cpliq;
P_bot = u(7); %Pa
P_top = u(8); %Pa

%SWITCH = Param(1);
Vm = (8.314*x(6))/((P_bot+P_top)/2);% molar specific volume [m3/mol]
% outputs:
sys(1) = x(1); %Fliq_in
sys(2) = x(2); %Fvap_in
sys(3) = x(3); %Ca(CO2g) or P_co2
sys(4) = x(4); %Cb(MEA)
sys(5) = x(5); %Cc(MEA-CO2 binded)
sys(6) = x(6); %T
sys(7) = x(3)*Vm; %Output in mol%
% sys(8) = (PB*F*x(3)+PC*F*x(5)-PA*F*(u(2)+1)-PD*F*(u(8)+1)-PHEAT*(UA*((u(4)+200)-x(4))));

% ******************************************
% Derivatives
% ******************************************
function sys = mdlDerivatives(t,x,u,Param)

global DH1 DH2 k01 k02 EA1 EA2 alpha1 alpha2 beta1 beta2;
global P_bot P_top Vliq Vvap rou Cpliq;

Vsection = 1.872*0.25*pi*(25.4/100)^2;
Vliq = Vsection*0.035;
Vvap = Vsection*0.985*(1-0.035);
rou = 1004; %kg/m3
Cpliq = 3414.5; % J/(kg*K)
DH1 = -32000; % J/ mol(CO2)
DH2 = 32000;

EA1 = 69050; %J/mol
EA2 = EA1+abs(DH1); %J/mol

%SWITCH = Param(1);
k01    = Param(1);
alpha1 = Param(2);
beta1  = Param(3);
k02    = Param(4);
alpha2 = Param(5);
beta2  = Param(6);
P_bot = u(7); %Pa
P_top = u(8); %Pa
% entering liq flow
Fliq_in = u(1);
% entering vap flow
Fvap_in = u(2);
% Entering concentration of A:
Ca   = u(3);
% Entering concentration of B:
Cb   = u(4);
% Entering concentration of C:
Cc   = u(5);
% Feed Temperature(liq)
T_in = u(6);

%Derivatives:
k1 = k01 * exp (-EA1/8.314/(alpha1*x(6)+beta1));
k2 = k02 * exp (-EA2/8.314/(alpha2*x(6)+beta2));
R1 = k1*x(3)*x(4); % mol/(m3 s)
R2 = k2*x(5);
Vm = (8.314*x(6))/((P_bot+P_top)/2); % molar specific volume [m3/mol]
sys(1) = Fliq_in-x(1); % Fliq out
sys(2) = Fvap_in-x(2)+(R2-R1)*Vliq*Vm; % Fvap out
sys(3) = Fvap_in/Vvap*(Ca-x(3))+(-R1+R2)*(Vliq/Vvap); %Ca
sys(4) = Fliq_in/Vliq*(Cb-x(4))-R1+R2; % Cb
sys(5) = Fliq_in/Vliq*(Cc-x(5))+R1-R2; %Cc
sys(6) = (Fliq_in/Vliq)*(T_in-298)-(x(1)/Vliq)*(x(6)-298)+(1/(rou*Cpliq))*(R1*(-DH1)+R2*(-DH2));