close all, clear all, clc
%% System declaration
load("DS_002_RCpar.mat")
Qn = 3.08;
R0 = R0_c;
R1 = R1_c;
R2 = R2_c;
C1 = C1_c;
C2 = C2_c;
Ts = 1;
a11 = 1-Ts/(R1*C1);
a22 = 1-Ts/(R2*C2);
bu1 = Ts/C1;
bu2 = Ts/C2;
bu3 = -Ts/(3600*Qn);
A = [a11 0 0;0 a22 0;0 0 1];
Bu = [bu1;bu2;bu3];
Du = [-R0];
ev = 3*sqrt(2);
Bv = 0.0154*Bu;
Dv = [0.036];

epsi = 0;
Bphi = [0;0;0];
Dphi = [1];

ord = length(pVoc)-1;
dphi = [ord:-1:1]*(pVoc(1:ord).*0.5.^[ord-1:-1:0])';
C = [-1 -1 dphi];
%% LMI variables declaration
nx = size(A,1);
nu = size(Bu,2);
nv = size(Bv,2);
ny = size(C,1);
nphi = size(Bphi,2);

P = sdpvar(nx,nx,'sym');
G = sdpvar(nx,nx,'full');
% W = sdpvar(nx,ny);
mu = sdpvar(1,1);
kap = sdpvar(1,1);

X3 = sdpvar(1,[5 1;1 0;2 0]);

options = sdpsettings('verbose',0,'solver','mosek');