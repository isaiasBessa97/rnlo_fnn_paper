close all, clear all, clc
%% Loading data
data = load("dataset\BID001_HPPC_31052024.xlsx");
z = data(:,2)';
u = data(:,3)';
%% Loading parameters
load("DS_002_RCpar.mat")
R1 = R1_c;
R2 = R2_c;
C1 = C1_c;
C2 = C2_c;
Qn = 3.08;
Ts = 1; %Sample time
%% Discrete model parameters 
a11 = 1-Ts/(R1*C1);
a22 = 1-Ts/(R2*C2);
b1 = Ts/C1;
b2 = Ts/C2;
b3 = -Ts/(3600*Qn);
A = [a11 0   0;
     0   a22 0;
     0   0   1];
B = [b1;
     b2;
     b3];
