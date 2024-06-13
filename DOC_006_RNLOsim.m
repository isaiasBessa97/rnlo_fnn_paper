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
nw = sqrt(2);
Bw = 0.308*Bu;
Dw = [0.042];
Bphi = [0;0;0];
Dphi = [1];
C = [-1 -1 0];
%% Loading datas
data = load("dataset\BID002_HPPC_01062024.xlsx");
nf = length(data);
vt = data(2:nf-100,2);
it = data(2:nf-100,3);
soc(1) = 1;
for ii = 2:length(vt)
    soc(ii) = soc(ii-1) - (1/(3600*Qn))*it(ii);
end
load("DS_003_rnloGain.mat")
%% Configuring input signals
t = 0:length(vt)-1;
u = it;
w = nw*rand(length(t),1);
%% Initial conditions
x1_hat(:,1) = [0;0;0.5];
ordR0 = length(pR0)-1;
cond = "c"; %"c" to R0 constant, "p" to R0 polinomial 
if cond == "c"
    R0(1) = R0_c;
elseif cond == "p"
    R0(1) = pR0*(x1_hat(3).^[ordR0:-1:0])';
end
Du = [-R0(1)];
voc1_hat(1) = pVoc*(x1_hat(3).^[ordR0:-1:0])';
y1_hat(1) = C*x1_hat(:,1) + Du*u(1) + voc1_hat;
%% Simulation
for ii = 2:length(vt)
    x1_hat(:,ii) = A*x1_hat(:,ii-1)+Bu*u(ii-1)+L1*(vt(ii-1)-y1_hat(ii-1));
    if cond == "c"
        R0(ii) = R0_c;
    elseif cond == "p"
        R0(ii) = pR0*(x1_hat(3,ii).^[ordR0:-1:0])';
    end    
    Du = [-R0(ii)];
    voc1_hat(ii) = pVoc*(x1_hat(3,ii).^[ordR0:-1:0])';
    y1_hat(ii) = C*x1_hat(:,ii)+Du*u(ii)+voc1_hat(ii);
end
%% Plot results
figure()
plot(t,vt,"k-","linewidth",2)
hold on
plot(t,y1_hat,"r-.","linewidth",2)
hold off
set(gca,"TickLabelInterpreter","latex","FontSize",16)
xlabel("Time (s)","FontSize",16,"Interpreter","latex")
ylabel("Voltage (V)","FontSize",16,"Interpreter","latex")
legend({"Measured","RNLO"},"Fontsize",14,"interpreter","latex")
ylim([2.5 4.3])

figure()
plot(t,soc,"k-","linewidth",2)
hold on
plot(t,x1_hat(3,:),"r-.","linewidth",2)
hold off
set(gca,"TickLabelInterpreter","latex","FontSize",16)
xlabel("Time (s)","FontSize",16,"Interpreter","latex")
ylabel("SOC","FontSize",16,"Interpreter","latex")
legend({"Measured","RNLO"},"Fontsize",14,"interpreter","latex")
%%
save("DS_007_RNLOresult","y1_hat","x1_hat")