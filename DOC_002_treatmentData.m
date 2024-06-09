close all, clear all, clc
%% Description
% This script is responsabel for conditioner the datas collected and
% separate each pulse in five points:
% (a): Last point before input current be connected
% (b): First point after input current be connected
% (c): Last point before input current be desconnected
% (d): First point after input current be desconnected
% (e): Point where the battery is rested
%% Constant parameters
Qn = 3.08; %Nominal capacity
%% Open data for BID001
data1 = load("dataset\BID001_HPPC_31052024.xlsx");
t1 = data1(:,1);
vt1 = data1(:,2);
it1 = data1(:,3);
soc1(1) = 1;
for ii = 2:length(t1)
    soc1(ii) = soc1(ii-1) - (1/(3600*Qn))*it1(ii);
end
pa1 = [5 656:652:46948];
pb1 = pa1+1;
pc1 = [56 pb1(2:end) + 53];
pd1 = pc1 + 1;
pe1 = pa1 - 1;
Pp1 = [pa1;pb1;pc1;pd1;pe1];
%% Open data for BID002
data2 = load("dataset\BID002_HPPC_01062024.xlsx");
t2 = data2(:,1);
vt2 = data2(:,2);
it2 = data2(:,3);
soc2(1) = 1;
for ii = 2:length(t2)
    soc2(ii) = soc2(ii-1) - (1/(3600*Qn))*it2(ii);
end
pa2 = [5 655:652:35211 35862 36512 37163 37813 38464 39115 39767 40416 ...
       41066 41716 42366 43016 43666 44315 44965 45615];
pb2 = pa2+1;
pc2 = [57 [pb2(2:55) + 53] 35916 36566 37217 37867 38518 39169 39821 ...
    40470 41120 41770 42420 43070 43720 44369 45019 45669];
pd2 = pc2 + 1;
pe2 = pa2 - 1;
Pp2 = [pa2;pb2;pc2;pd2;pe2];