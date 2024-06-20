close all, clear all, clc
%% Setting colors
c1 = [001 156 225]/255;
c2 = [000 086 224]/255;
c3 = [000 224 144]/255;
c4 = [000 015 224]/255;
%% Loading data
data = load("dataset\BID004_RANDCh_30052024.xlsx");
N = length(data);
vt = data(2:N-60,2)';
it = data(2:N-60,3)';
soc(1) = 1;
Qn = 3.08;
Ts = 1; %Sample time
for ii = 2:length(it)
    soc(ii) = soc(ii-1)-Ts*it(ii)/(3600*Qn);
end
t = 0:Ts:length(it)-1;
%% Loading parameters
load("DS_002_RCpar.mat")
R1 = R1_c;
R2 = R2_c;
C1 = C1_c;
C2 = C2_c;
Np = length(pVoc)-1; %Polynomial degree
pdVoc = [Np:-1:1].*pVoc(1:9);
pdR0 = [Np:-1:1].*pR0(1:9);
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
%% Initialization
xsoc_0 = [0.50 0.75 1.25 1.50];
for jj = 1:length(xsoc_0)
u(1) = it(1);
x_hat(:,1) = [0;0;xsoc_0(jj)];
soc_hat{jj}(1) = x_hat(3,1);
e_hat{jj}(1) = soc(1)-x_hat(3,1);
R0_hat(1) = pR0*((x_hat(3,1)).^[Np:-1:0])';
voc_hat(1) = pVoc*((x_hat(3,1)).^[Np:-1:0])';
y_hat(1) = -x_hat(1,1)-x_hat(2,1)-R0_hat(1)*u(1)+voc_hat(1);
P{1} = [0.001 0 0;0 0.001 0;0 0 0.1];
R = 0.2;
Q = 1e-6*eye(3);
In = eye(size(A,1));
%% Extrapolation/Prediction
x_ext(:,1) = A*x_hat(:,1) + B*u(1);
P_ext{1} = A*P{1}*A' + Q;
    for kk = 2:length(it)
        %% Measure
        z(kk) = vt(kk);
        u(kk) = it(kk);
        %% Update/Correction
        dR0(kk-1) = pdR0*((x_ext(3,kk-1)).^[Np-1:-1:0])';
        dVoc(kk-1) = pdVoc*((x_ext(3,kk-1)).^[Np-1:-1:0])';
        R0(kk) = pR0*((x_ext(3,kk-1)).^[Np:-1:0])';
        voc(kk) = pVoc*((x_ext(3,kk-1)).^[Np:-1:0])';
        C = [-1 -1 (-dR0(kk-1)*u(kk-1)+dVoc(kk-1))];
        g(kk) = -x_ext(1,kk-1)-x_ext(2,kk-1)-R0(kk)*u(kk-1)+voc(kk);
    
        K{kk} = P_ext{kk-1}*C'*inv(C*P_ext{kk-1}*C'+R);
        x_hat(:,kk) = x_ext(:,kk-1) + K{kk}*(z(kk)-g(kk));
        P{kk} = (In-K{kk}*C)*P_ext{kk-1}*(In-K{kk}*C)'+K{kk}*R*K{kk}';
    
    
        R0_hat(kk) = pR0*((x_hat(3,kk)).^[Np:-1:0])';
        voc_hat(kk) = pVoc*((x_hat(3,kk)).^[Np:-1:0])';
        y_hat(kk) = -x_hat(1,kk)-x_hat(2,kk)-R0_hat(kk)*u(kk)+voc_hat(kk);
    
        e_hat{jj}(kk) = soc(kk) - x_hat(3,kk) ;
        soc_hat{jj}(kk) = x_hat(3,kk);
        %% Extrapolation/Prediction
        x_ext(:,kk) = A*x_hat(:,kk) + B*u(kk);
        P_ext{kk} = A*P{kk}*A' + Q;
    end
end
%% Plot results

figure()
plot(t,soc,"k-","linewidth",2)
hold on
plot(t,soc_hat{1},":","linewidth",2,"Color",c1)
plot(t,soc_hat{2},":","linewidth",2,"Color",c2)
plot(t,soc_hat{3},":","linewidth",2,"Color",c3)
plot(t,soc_hat{4},":","linewidth",2,"Color",c4)
hold off
ylim([0 1.5])
xlim([0 length(it)])
set(gca,"TickLabelInterpreter","latex","FontSize",20)
ylabel("SOC","Interpreter","latex","FontSize",20)
xlabel("Time (s)","Interpreter","latex","FontSize",20)
legend({"Measured","$x_0 = 1.50$","$x_0 = 1.25$","$x_0 = 0.50$",...
        "$x_0 = 0.25$"},"interpreter","latex","fontsize",16)
axes('Position',[.20 .24 .22 .22])
indexOfInterest = (t >= 0) & (t <= 1000);
plot(t(indexOfInterest),soc(indexOfInterest),"k-","linewidth",2)
hold on
plot(t(indexOfInterest),soc_hat{1}(indexOfInterest),":","linewidth",2,"Color",c1)
plot(t(indexOfInterest),soc_hat{2}(indexOfInterest),":","linewidth",2,"Color",c2)
plot(t(indexOfInterest),soc_hat{3}(indexOfInterest),":","linewidth",2,"Color",c3)
plot(t(indexOfInterest),soc_hat{4}(indexOfInterest),":","linewidth",2,"Color",c4)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',14)

figure()
plot(t,e_hat{1},"-.","linewidth",2,"Color",c1)
hold on
plot(t,e_hat{2},"-.","linewidth",2,"Color",c2)
plot(t,e_hat{3},"-.","linewidth",2,"Color",c3)
plot(t,e_hat{4},"-.","linewidth",2,"Color",c4)
hold off
ylim([-0.50 0.50])
xlim([0 length(it)])
set(gca,"TickLabelInterpreter","latex","FontSize",20)
ylabel("$e_{\mathrm{SOC}}$","Interpreter","latex","FontSize",20)
xlabel("Time (s)","Interpreter","latex","FontSize",20)
legend({"$x_0 = 1.50$","$x_0 = 1.25$","$x_0 = 0.50$",...
        "$x_0 = 0.25$"},"interpreter","latex","fontsize",16)
axes('Position',[.24 .24 .24 .24])
indexOfInterest = (t >= 0) & (t <= 1000);
plot(t(indexOfInterest),e_hat{1}(indexOfInterest),":","linewidth",2,"Color",c1)
hold on
plot(t(indexOfInterest),e_hat{2}(indexOfInterest),":","linewidth",2,"Color",c2)
plot(t(indexOfInterest),e_hat{3}(indexOfInterest),":","linewidth",2,"Color",c3)
plot(t(indexOfInterest),e_hat{4}(indexOfInterest),":","linewidth",2,"Color",c4)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',14)