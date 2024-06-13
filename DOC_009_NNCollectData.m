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
%% Loading datas
directory = "dataset";
filesAndFolders = dir(directory);
filesInDir = filesAndFolders(~([filesAndFolders.isdir]))
batTest = "RANDCh";
numOfFiles = length(filesInDir);
ii = 1;
kk = 0;
while (ii <= numOfFiles-1)
    fileName = filesInDir(ii).name;
    if strfind(fileName,batTest) == 8
        data = load(directory + "\" + fileName);
        nf = length(data);
        kk = kk+1;
        vt{kk} = data(3:nf-600,2); % Terminal voltage
        it{kk} = data(3:nf-600,3); % Input current
        soc{kk}(1) = 1;
        for jj = 2:length(it{kk})
            soc{kk}(jj) = soc{kk}(jj-1) - (1/(3600*Qn))*it{kk}(jj);
        end
    end
    ii = ii + 1;
end
vtTrain = {vt{1};vt{2};vt{3};vt{5}};
itTrain = {it{1};it{2};it{3};it{5}};
socTrain = {soc{1};soc{2};soc{3};soc{5}};
vtTest = {vt{4};vt{6}};
itTest = {it{4};it{6}};
socTest = {soc{4};soc{6}};
load("DS_004_rnlonnGain.mat")
%% Collecting training data
for ii = 1:length(vtTrain)
    u{ii} = itTrain{ii};
    t{ii} = 0:length(u{ii})-1;
    psi{ii} = 1*(epsi*2/sqrt(2))*(rand(length(itTrain{ii}),1)-0.5);
    x2_hat{ii}(:,1) = [0;0;0.5];
    targ{ii}(1) = soc{ii}(1)-x2_hat{ii}(3,1);
    ordR0 = length(pR0)-1;
    cond = "c"; %"c" to R0 constant, "p" to R0 polinomial 
    if cond == "c"
        R0(1) = R0_c;
    elseif cond == "p"
        R0(1) = pR0*(x2_hat(3).^[ordR0:-1:0])';
    end
    Du = [-R0(1)];
    voc2_hat{ii}(1) = pVoc*(x2_hat{ii}(3).^[ordR0:-1:0])';
    y2_hat{ii}(1) = C*x2_hat{ii}(:,1) + Du*u{ii}(1) + voc2_hat{ii}(1);
    for kk = 2:length(vtTrain{ii})
        x2_hat{ii}(:,kk) = A*x2_hat{ii}(:,kk-1)+Bu*u{ii}(kk-1)+ ...
                       L2*(vtTrain{ii}(kk-1)-y2_hat{ii}(kk-1))+ ...
                       Bpsi*psi{ii}(kk-1);
        targ{ii}(kk) = socTrain{ii}(kk)-x2_hat{ii}(3,kk);
        if cond == "c"
            R0(kk) = R0_c;
        elseif cond == "p"
            R0(kk) = pR0*(x2_hat{ii}(3,kk).^[ordR0:-1:0])';
        end    
        Du = [-R0(kk)];
        voc2_hat{ii}(kk) = pVoc*(x2_hat{ii}(3,kk).^[ordR0:-1:0])';
        y2_hat{ii}(kk) = C*x2_hat{ii}(:,kk)+Du*u{ii}(kk)+voc2_hat{ii}(kk);
    end 
end
uTrain = [u{1}' u{2}' u{3}' u{4}'];
yTrain = [vtTrain{1}' vtTrain{2}' vtTrain{3}' vtTrain{4}'];
yTrain_hat = [y2_hat{1} y2_hat{2} y2_hat{3} y2_hat{4}];
xTrain_hat = [x2_hat{1} x2_hat{2} x2_hat{3} x2_hat{4}];
psiTrain = [psi{1}' psi{2}' psi{3}' psi{4}'];
targTrain = [targ{1} targ{2} targ{3} targ{4}];
%% Saving data
save("DS_005_RNLONN_dataTrain","uTrain","yTrain","yTrain_hat",...
     "xTrain_hat","psiTrain","targTrain")