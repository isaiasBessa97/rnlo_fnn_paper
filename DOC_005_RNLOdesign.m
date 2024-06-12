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
a11 = 1-Ts/(R1_c*C1_c);
a22 = 1-Ts/(R2_c*C2_c);
bu1 = Ts/C1_c;
bu2 = Ts/C2_c;
bu3 = -Ts/(3600*Qn);
A = [a11 0 0;0 a22 0;0 0 1];
Bu = [bu1;bu2;bu3];
Du = [-R0];
nw = sqrt(2);
Bw = 0.308*Bu;
Dw = [0.042];
Bphi = [0;0;0];
Dphi = [1];

ord = length(pVoc)-1;
dphi = [ord:-1:1]*(pVoc(1:ord).*0.5.^[ord-1:-1:0])';
C = [-1 -1 dphi];
%% LMI variables declaration
nx = size(A,1);
nu = size(Bu,2);
nw = size(Bw,2);
ny = size(C,1);
nphi = size(Bphi,2);

P = sdpvar(nx,nx,'sym');
G = sdpvar(nx,nx,'full');
W = sdpvar(nx,ny);
mu = sdpvar(1,1);
la = sdpvar(1,1);
options = sdpsettings('verbose',0,'solver','mosek');
%% Constraints
x1_lim = [-6e-2 6e-2];
x2_lim = [-6e-2 6e-2];
x3_lim = [-0.9 0.9];

ak = [1/x1_lim(1) 1/x1_lim(2) 0           0           0           0;
      0           0           1/x2_lim(1) 1/x2_lim(2) 0           0;
      0           0           0           0           1/x3_lim(1) 1/x3_lim(2)];
%% Parameters settings
kap = 0.9;
tau = 1e-10;
lc = 1e-20;
P0 = [(2/6e-2)^2 0 0;
      0 (2/6e-2)^2 0;
      0 0 (2/0.9)^2];
%% LMI formulation
LMIs = [];
LMIs = [la >= 0] + [mu >= 0] + LMIs;
LMI1 = [P - la*eye(nx) >= 0];
LMIs = LMIs + LMI1;
for ii = 1:size(ak,2)
    LMI2 = [[1+kap            (1+kap)*ak(:,ii)';
             (1+kap)*ak(:,ii) P] >= 0];
    LMIs = LMIs + LMI2;
end
zxw = zeros(nx,nw);
zxp = zeros(nx,nphi);
zwp = zeros(nw,nphi);
LMI3 = [[P-G-G'          G*A+W*C                  G*Bw+W*Dw      G*Bphi+W*Dphi;
        (G*A+W*C)'       (tau-1)*P+mu*lc*eye(nx)  zxw            zxp;
        (G*Bw+W*Dw)'     zxw'                    -tau*eye(nw)    zwp;
        (G*Bphi+W*Dphi)' zxp'                     zwp'          -mu*eye(nphi)] <= 0];
LMI4 = [(1+kap)*P0 - P >= 0];
LMIs = LMIs + LMI4;
%% Solving
result = optimize(LMIs,-la,options)
P = double(P);
G = double(G);
W = double(W);
mu = double(mu);
la = double(la);
L = -inv(G)*W;
W'