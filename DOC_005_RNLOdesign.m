close all, clear all, clc
%% Description
% Script developed to design an observer gain L ISS with respect to
% disturbance signal
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
W = sdpvar(nx,ny);
mu = sdpvar(1,1);
kap = sdpvar(1,1);
options = sdpsettings('verbose',0,'solver','mosek');
%% Constraints
x1_lim = [-1e-1 1e-1];
x2_lim = [-1e-1 1e-1];
x3_lim = [-0.9 0.9];

ak = [1/x1_lim(1) 1/x1_lim(2) 0           0           0           0;
      0           0           1/x2_lim(1) 1/x2_lim(2) 0           0;
      0           0           0           0           1/x3_lim(1) 1/x3_lim(2)];
%% Parameters settings
la1 = 2;
tau1 = 1e-2;
lc1 = 6e-1;%6e-1
P01 = [(4/x1_lim(2))^2 0 0;
       0 (4/x2_lim(2))^2 0;
       0 0 (4/x3_lim(2))^2];
%% LMI formulation
LMIs = [];
LMIs = [kap >= 0] + [mu >= 0] + LMIs;
LMI1 = [P - la1*eye(nx) >= 0];
LMIs = LMIs + LMI1;
for ii = 1:size(ak,2)
    LMI2 = [[1+kap            (1+kap)*ak(:,ii)';
             (1+kap)*ak(:,ii) P] >= 0];
    LMIs = LMIs + LMI2;
end
zxw = zeros(nx,nv);
zxp = zeros(nx,nphi);
zwp = zeros(nv,nphi);
LMI3 = [[P-G-G'          G*A+W*C                  G*Bv+W*Dv      G*Bphi+W*Dphi;
        (G*A+W*C)'       (tau1-1)*P+mu*lc1*eye(nx)  zxw            zxp;
        (G*Bv+W*Dv)'     zxw'                    -tau1*eye(nv)    zwp;
        (G*Bphi+W*Dphi)' zxp'                     zwp'          -mu*eye(nphi)] <= 0];
LMI4 = [(1+kap)*P01 - P >= 0];
LMIs = LMIs + LMI3 +  LMI4;
%% Solving
result = optimize(LMIs,-kap,options)
P1 = double(P);
G1 = double(G);
W1 = double(W);
mu1 = double(mu);
kap1 = double(kap)
L1 = -inv(G1)*W1;
%% Save result
save("DS_003_rnloGain","P1","G1","W1","mu1","la1","L1","P01","kap1",...
     "lc1","tau1","ev","epsi")