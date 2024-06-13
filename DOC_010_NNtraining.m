close all, clear all, clc
%% Open data
load("DS_005_RNLONN_dataTrain.mat")
%% Dataset configuration
d = 1; % delay size
m = 7; % buffer size
for ii = 1:length(uTrain)
    n = length(uTrain{ii}); %original size input
    for jj = 1+m:n
    %% Buffer design
        %Measured values
        ub(1:m+d) =  uTrain{ii}(jj:-1:jj-m);
        yb(1:m+d) =  yTrain{ii}(jj:-1:jj-m);
        %Estimated values
        x3b_hat(1:m) = xTrain_hat{ii}(3,jj-d:-1:jj-m);
        yb_hat(1:m) =  yTrain_hat{ii}(jj-d:-1:jj-m);
        %Randomly generate auxiliary signal
        pcb(1:m) = psiTrain{ii}(jj-d:-1:jj-m);
    %% Target design
        tb1 = targTrain{ii}(jj);
    %% Input (nx) and output (ny) length
        nx = length(ub)+length(yb)+length(x3b_hat)+length(yb_hat)+length(pcb);
        ny = length(tb1);
    %% Input (X) and Output (Y) to training
        X{jj-m} = [ub';yb';x3b_hat';yb_hat';pcb'];
        T{jj-m} = [tb1];
    end
end
%% Train procedure
%Network structure
net = feedforwardnet([5],'trainlm');
%Output activation function
net.layers{2}.transferFcn = 'tansig'
net.layers{1}.transferFcn = 'purelin'
%Setting parameters
net.trainParam.min_grad = 1e-15;
net.trainParam.goal = 1e-8;
net.trainParam.epochs = 500;
%Train network
net = train(net,X,T);
%% Save model
save("DS_006_modelFNN_v0","net","m","d")