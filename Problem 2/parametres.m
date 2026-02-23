clear all;
close all;


% Paràmetres de la xarxa d'alterna

Ap=90e3*sqrt(2)/sqrt(3)
f=50;
sc=0; 

%Paràmetres de connexió a xarxa
Rx=4.8*0.669
Lx=0.052*0.669

%Temps de mostreig
T=1/12e3;


%Impedància d'acoblament del convertidor
Lac=30e-3
Rac=3.68/10


%Paràmetres del controlador del llaç de corrent
Ll=Lac;
rl=Rac;
tau=1/3100;
Kp=Ll/tau;
Ki=rl/tau;

t0=0;
t1=1.5;
% max_step=T/(20);
% min_step=T/(20)/100;
% tol=1e-3;


%Consignes de potència reactiva
ida3=0
ida4=0


% Paràmetres de la xarxa de contínua
R1=0.5;
R2=0.25;
R3=0.4;
L1=0.005;
L2=0.0025;
L3=0.004
K3=1/20
K4=K3
P1=100e6
P2=100e6
E0=145e3
id=0


%Entrada de potencia
it=[0 0.2 0.200001 0.35 0.3500001 t1];
iv=[0 0 100e6 100e6 0 0];


