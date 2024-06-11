close all, clear all, clc
%% Discription
% Script developed to computer parameters of the 2nd order ECM for
% lithium-ion batteries:
% R0: internal resistance
% R1,C1: electrochemical  polarization effect
% R2,C2: concentration polarization effect
% voc: Open-circuit voltage
%% Load data
load("DS_001_treatedData.mat")
%% Separating datas for R0 and voc
for ii = 1:length(Pp)
    for jj = 1:length(Itp{ii})
        R0_meas_a{ii}(jj) = abs(Vtp{ii}(1,jj) - Vtp{ii}(2,jj))/...
                    abs(Itp{ii}(1,jj) - Itp{ii}(2,jj));
        R0_meas_c{ii}(jj) = abs(Vtp{ii}(3,jj) - Vtp{ii}(4,jj))/...
                    abs(Itp{ii}(3,jj) - Itp{ii}(4,jj));
    end 
end
R0_meas = [R0_meas_a{1} R0_meas_a{2} R0_meas_a{3} R0_meas_a{4}...
           R0_meas_c{1} R0_meas_c{2} R0_meas_c{3} R0_meas_c{4}];
soc_R0 = [SOCp{1}(1,:) SOCp{2}(1,:) SOCp{3}(1,:) SOCp{4}(1,:)...
          SOCp{1}(3,:) SOCp{2}(3,:) SOCp{3}(3,:) SOCp{4}(3,:)];

voc_meas = [Vtp{1}(5,:) Vtp{2}(5,:) Vtp{3}(5,:) Vtp{4}(5,:)];
soc_voc = [SOCp{1}(5,:) SOCp{2}(5,:) SOCp{3}(5,:) SOCp{4}(5,:)];
%% R0 identification
% Two types of R0: constant (R0_c) and polynomial (pR0)
ordR0 = 9;
R0_c = mean(R0_meas);
pR0 = polyfit(soc_R0,R0_meas,ordR0);
%% voc identification
ordVoc = 9;
pVoc = polyfit(soc_voc,voc_meas,ordVoc);
%% SOC simulated
soc_sim = 0.1:0.001:1;

X_R0 = soc_sim'.^[ordR0:-1:0];
R0_p = pR0*X_R0';

X_voc = soc_sim'.^[ordVoc:-1:0];
voc_p = pVoc*X_voc';
%% R1,R2,C1,C2 identifaction setting
fOpt = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[2.2 0 0 0 0],...
                'Upper',[4.5 1e-3 1e-5 5e1 5e-2],...
                'StartPoint',[0.7988 0.0027 0.5212 0.7114 0.2222]); %[0.7988 0.0027 0.5212 0.7114 0.2222]
fTyp = fittype('a0-a1*exp(-b1*x)-a2*exp(-b2*x)','options',fOpt);
%% R1,R2,C1,C2 identifaction procedure
for ii = 1:length(Pp)
    for jj = 1:length(Itp{ii})-1
        Vce = Vt{ii}(Pp{ii}(3,jj):Pp{ii}(5,jj+1));
        tce = linspace(0,length(Vce),length(Vce));
        pVce = fit(tce',Vce,fTyp);
        V1 = max([pVce.a1 pVce.a2]);
        V2 = min([pVce.a1 pVce.a2]);
        tau1 = inv(max([pVce.b1 pVce.b2]));
        tau2 = inv(min([pVce.b1 pVce.b2]));
        R1_meas_ce{ii}(jj) = V1/(Itp{ii}(3,jj)*...
            (1-exp(-((Tp{1}(3,1)-Tp{1}(1,1))/tau1))));
        C1_meas_ce{ii}(jj) = tau1/R1_meas_ce{ii}(jj);
        R2_meas_ce{ii}(jj) = V2/(Itp{ii}(3,jj)*...
            (1-exp(-((Tp{1}(3,1)-Tp{1}(1,1))/tau2))));
        C2_meas_ce{ii}(jj) = tau2/R2_meas_ce{ii}(jj);
    end
end
R1_meas = [R1_meas_ce{1} R1_meas_ce{2} R1_meas_ce{3} R1_meas_ce{4}];
R2_meas = [R2_meas_ce{1} R2_meas_ce{2} R2_meas_ce{3} R2_meas_ce{4}];
C1_meas = [C1_meas_ce{1} C1_meas_ce{2} C1_meas_ce{3} C1_meas_ce{4}];
C2_meas = [C2_meas_ce{1} C2_meas_ce{2} C2_meas_ce{3} C2_meas_ce{4}];

R1_c = mean(R1_meas);
R2_c = mean(R2_meas);
C1_c = mean(C1_meas);
C2_c = mean(C2_meas);
%% Plot figures
figure()
plot(soc_R0,R0_meas,"rx",'LineWidth',2)
hold on
plot(soc_sim,R0_p,'k-','linewidth',2)
plot(soc_sim,R0_c*ones(length(soc_sim),1),'b-','linewidth',2)
hold off
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel("SOC",'Interpreter','latex','FontSize',16)
ylabel("$R_0$ ($\Omega$)",'Interpreter','latex','FontSize',16)
legend({"Measured","Poly","Const"},"interpreter","latex","fontsize",14)

figure()
plot(soc_voc,voc_meas,"rx",'LineWidth',2)
hold on
plot(soc_sim,voc_p,'k-','linewidth',2)
hold off
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel("SOC",'Interpreter','latex','FontSize',16)
ylabel("OCV (V)",'Interpreter','latex','FontSize',16)
legend({"Measured","Model"},"interpreter","latex","fontsize",14,...
       "location","northwest")
%% Saving data
save("DS_002_RCpar","pR0","pVoc","R0_c","R1_c","R2_c","C1_c","C2_c")