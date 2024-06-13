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

ev = sqrt(2);
Bv = 0.308*Bu;
Dv = [0.042];

epsi = 0.02*sqrt(2);
Bpsi = [0;0;1];
Dpsi = [0];

Bw = [ev*Bv epsi*Bpsi];
Dw = [ev*Dv epsi*Dpsi];

Bphi = [0;0;0];
Dphi = [1];

ord = length(pVoc)-1;
dphi = [ord:-1:1]*(pVoc(1:ord).*0.5.^[ord-1:-1:0])';
C = [-1 -1 0];
%% Open models
load("DS_006_modelFNN_v0.mat")
load("DS_004_rnlonnGain.mat")
load("DS_007_RNLOresult.mat")
data = load("dataset\BID004_HPPC_02062024.xlsx");
nf = length(data);
vt = data(2:nf-100,2);
it = data(2:nf-100,3);
soc(1) = 1;
for ii = 2:length(vt)
    soc(ii) = soc(ii-1) - (1/(3600*Qn))*it(ii);
end
%% Configuring input signals
t = 0:length(vt)-1;
u = it;
%% Initial conditions
x3_hat(:,1) = [0;0;0.5];
ordR0 = length(pR0)-1;
cond = "c"; %"c" to R0 constant, "p" to R0 polinomial 
if cond == "c"
    R0(1) = R0_c;
elseif cond == "p"
    R0(1) = pR0*(x3_hat(3).^[ordR0:-1:0])';
end
Du = [-R0(1)];
voc3_hat(1) = pVoc*(x3_hat(3).^[ordR0:-1:0])';
y3_hat(1) = C*x3_hat(:,1) + Du*u(1) + voc3_hat;
%% Simulation
for ii = 2:length(vt)
    if ii-m > 0
        ub(1:m+d) = u(ii:-1:ii-m);
        yb(1:m+d) = vt(ii:-1:ii-m);
        x1b_c(1:m) = x3_hat(1,ii-d:-1:ii-m);
        x2b_c(1:m) = x3_hat(2,ii-d:-1:ii-m);
        x3b_c(1:m) = x3_hat(3,ii-d:-1:ii-m);
        yb_c(1:m) = y3_hat(ii-d:-1:ii-m);
        pcb(1:m) = psi(ii-d:-1:ii-m);

        pc_cell = net({[ub';yb';x1b_c';x2b_c';x3b_c';yb_c';pcb']});
        psi(ii) = pc_cell{1}(1);
        if (psi(ii) >= 0.04) || (psi(ii) <= -0.04) 
            if psi(ii) >= 0.04
                psi(ii) = 0.04;
            else
                psi(ii) = -0.04;
            end
        else
            psi(ii) = psi(ii);
        end
    else
        psi(ii) = 0;           
    end

    x3_hat(:,ii) = A*x3_hat(:,ii-1)+Bu*u(ii-1)+...
                   L2*(vt(ii-1)-y3_hat(ii-1))+Bpsi*psi(ii-1);
    if cond == "c"
        R0(ii) = R0_c;
    elseif cond == "p"
        R0(ii) = pR0*(x3_hat(3,ii).^[ordR0:-1:0])';
    end    
    Du = [-R0(ii)];
    voc3_hat(ii) = pVoc*(x3_hat(3,ii).^[ordR0:-1:0])';
    y3_hat(ii) = C*x3_hat(:,ii)+Du*u(ii)+voc3_hat(ii);
end
%% Plot results
figure()
plot(t,vt,"k-","linewidth",2)
hold on
plot(t,y1_hat,"b-.","linewidth",2)
plot(t,y3_hat,"r:","linewidth",2)
hold off
set(gca,"TickLabelInterpreter","latex","FontSize",16)
xlabel("Time (s)","FontSize",16,"Interpreter","latex")
ylabel("Voltage (V)","FontSize",16,"Interpreter","latex")
legend({"Measured","RNLO","RNLO+NN"},"Fontsize",14,"interpreter","latex")
ylim([2.5 4.3])

figure()
plot(t,soc,"k-","linewidth",2)
hold on
plot(t,x1_hat(3,:),"b-.","linewidth",2)
plot(t,x3_hat(3,:),"r:","linewidth",2)
hold off
set(gca,"TickLabelInterpreter","latex","FontSize",16)
xlabel("Time (s)","FontSize",16,"Interpreter","latex")
ylabel("SOC","FontSize",16,"Interpreter","latex")
legend({"Measured","RNLO","RNLO+NN"},"Fontsize",14,"interpreter","latex")

figure()
plot(t,psi,"k-","linewidth",2)
set(gca,"TickLabelInterpreter","latex","FontSize",16)
xlabel("Time (s)","FontSize",16,"Interpreter","latex")
ylabel("$\psi$","FontSize",16,"Interpreter","latex")