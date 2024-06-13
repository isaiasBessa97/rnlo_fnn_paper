close all, clear all, clc
%% Open data
load("DS_005_RNLONN_dataTrain.mat")
%% Dataset configuration
d = 1; % delay size
m = 8; % buffer size
n = length(uTrain); %original size input
for ii = 1+m:n
%% Buffer design
    %Measured values
    ub(1:m+d) =  uTrain(ii:-1:ii-m);
    yb(1:m+d) =  yTrain(ii:-1:ii-m);
    %Estimated values
    x1b_hat(1:m) = xTrain_hat(1,ii-d:-1:ii-m);
    x2b_hat(1:m) = xTrain_hat(2,ii-d:-1:ii-m);
    x3b_hat(1:m) = xTrain_hat(3,ii-d:-1:ii-m);
    yb_hat(1:m) =  yTrain_hat(ii-d:-1:ii-m);
    %Randomly generate auxiliary signal
    pcb(1:m) = psiTrain(ii-d:-1:ii-m);
%% Target design
    tb1 = targTrain(ii);
%% Input (nx) and output (ny) length
    nx = length(ub)+length(yb)+length(x1b_hat)+length(x2b_hat)+...
        +length(x3b_hat)+length(yb_hat)+length(pcb);
    ny = length(tb1);
%% Input (X) and Output (Y) to training
    X{ii-m} = [ub';yb';x1b_hat';x2b_hat';x3b_hat';yb_hat';pcb'];
    T{ii-m} = [tb1];
end
%% Train procedure
%Network structure
net = feedforwardnet([10],'trainlm');
%Output activation function
net.layers{2}.transferFcn = 'tansig'
net.layers{1}.transferFcn = 'purelin'
%Setting parameters
net.trainParam.min_grad = 1e-12;
net.trainParam.goal = 1e-6;
net.trainParam.epochs = 200;
%Train network
net = train(net,X,T);
%% Save model
save("DS_006_modelFNN_v0","net","m","d")