close all, clear all, clc
%% Description
% Code cretated to open data collected in the testbanch
%% Choosing the directory
directory = "dataset";
filesAndFolders = dir(directory);
filesInDir = filesAndFolders(~([filesAndFolders.isdir]))
%% Tes code
batTest = "HPPC";
numOfFiles = length(filesInDir);
%% Initialization
ii = 1;
kk = 0;
Qn = 3.08; %Nominal capacity
%% Open datasets
while (ii <= numOfFiles-1)
    fileName = filesInDir(ii).name;
    if strfind(fileName,batTest) == 8
        data = load(directory + "\" + fileName);
        kk = kk+1;
        t{kk} = data(:,1); % Time
        vt{kk} = data(:,2); % Terminal voltage
        it{kk} = data(:,3); % Input current
        soc{kk}(1) = 1;
        for jj = 2:length(it{kk})
            soc{kk}(jj) = soc{kk}(jj-1) - (1/(3600*Qn))*it{kk}(jj);
        end
    end
    ii = ii + 1;
end
%% Plot figures
% Plot current
figure()
plot(t{1},it{1},'k-','linewidth',2)
hold on
plot(t{2},it{2},'r--','linewidth',2)
plot(t{3},it{3},'b-.','linewidth',2)
plot(t{4},it{4},'m:','linewidth',2)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',18)
xlabel("Time (s)","FontSize",20,"Interpreter","latex")
ylabel("Current (A)","Fontsize",20,"Interpreter","latex")
legend({"BID001","BID002","BID003","BID004"},"Fontsize",14,...
    "Interpreter","latex","location","southeast")
grid on, grid minor

% Plot voltage
figure()
plot(t{1},vt{1},'k-','linewidth',2)
hold on
plot(t{2},vt{2},'r--','linewidth',2)
plot(t{3},vt{3},'b-.','linewidth',2)
plot(t{4},vt{4},'m:','linewidth',2)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',18)
xlabel("Time (s)","FontSize",20,"Interpreter","latex")
ylabel("Voltage (V)","Fontsize",20,"Interpreter","latex")
legend({"BID001","BID002","BID003","BID004"},"Fontsize",14,...
    "Interpreter","latex","location","southeast")
grid on, grid minor

% Plot soc
figure()
plot(t{1},soc{1},'k-','linewidth',2)
hold on
plot(t{2},soc{2},'r--','linewidth',2)
plot(t{3},soc{3},'b-.','linewidth',2)
plot(t{4},soc{4},'m:','linewidth',2)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',18)
xlabel("Time (s)","FontSize",20,"Interpreter","latex")
ylabel("SOC","Fontsize",20,"Interpreter","latex")
legend({"BID001","BID002","BID003","BID004"},"Fontsize",14,...
    "Interpreter","latex","location","southwest")
grid on, grid minor