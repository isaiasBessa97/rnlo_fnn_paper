close all, clear all, clc
%% Description
% Script created to simulate the battery model identified using discrete
% equations.
% Sample time: 1s
%% Load data
data = load("dataset\BID003_HPPC_01062024.xlsx");
vt = data(:,2);
it = data(:,3);
soc(1) = 1;
Qn = 3.08;
for ii = 2:length(vt)
    soc(ii) = soc(ii-1) - it(ii)/(3600*Qn);
end
t = 0:1:length(vt)-1;
%% Discrete system
load("DS_002_RCpar.mat")
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
C = [-1 -1 0];
%% Initial conditions
u = it;
x(:,1) = [0;0;1];
ordR0 = length(pR0)-1;
cond = "c"; %"c" to R0 constant, "p" to R0 polinomial 
if cond == "c"
    R0(1) = R0_c;
elseif cond == "p"
    R0(1) = pR0*(x(3).^[ordR0:-1:0])';
end
Du = [-R0(1)];
voc(1) = pVoc*(x(3).^[ordR0:-1:0])';
y(1) = C*x(:,1) + Du*u(1) + voc;
for ii = 2:length(vt)
    x(:,ii) = A*x(:,ii-1)+Bu*u(ii-1);
    if cond == "c"
        R0(ii) = R0_c;
    elseif cond == "p"
        R0(ii) = pR0*(x(3,ii).^[ordR0:-1:0])';
    end    
    Du = [-R0(ii)];
    voc(ii) = pVoc*(x(3,ii).^[ordR0:-1:0])';
    y(ii) = C*x(:,ii)+Du*u(ii)+voc(ii);
end
%% Plot results
figure()
plot(t,vt,"k-","linewidth",2)
hold on
plot(t,y,"r-.","linewidth",2)
hold off
set(gca,"TickLabelInterpreter","latex","FontSize",16)
xlabel("Time (s)","FontSize",16,"Interpreter","latex")
ylabel("Voltage (V)","FontSize",16,"Interpreter","latex")
legend({"Measured","Model"},"Fontsize",14,"interpreter","latex")
ylim([2.5 4.2])