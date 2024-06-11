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
pc1 = [58 pb1(2:end) + 53];
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
%% Open data for BID003
data3 = load("dataset\BID003_HPPC_01062024.xlsx");
t3 = data3(:,1);
vt3 = data3(:,2);
it3 = data3(:,3);
soc3(1) = 1;
for ii = 2:length(t3)
    soc3(ii) = soc3(ii-1) - (1/(3600*Qn))*it3(ii);
end
pa3 = [5 656:652:46296 46947:651:48249];
pb3 = pa3 + 1;
pc3 = [58 pa3(2:end) + 54];
pd3 = pc3 + 1;
pe3 = pa3 - 1;
Pp3 = [pa3;pb3;pc3;pd3;pe3];
%% Open data for BID004
data4 = load("dataset\BID004_HPPC_02062024.xlsx");
t4 = data4(:,1);
vt4 = data4(:,2);
it4 = data4(:,3);
soc4(1) = 1;
for ii = 2:length(t4)
    soc4(ii) = soc4(ii-1) - (1/(3600*Qn))*it4(ii);
end
pa4 = [5 656:652:1960 2613:653:48323];
pb4 = pa4 + 1;
pc4 = [58 pb4(2:end) + 53];
pd4 = pc4 + 1;
pe4 = pa4 - 1;
Pp4 = [pa4;pb4;pc4;pd4;pe4];
%% Agrouping data
Pp = {Pp1;Pp2;Pp3;Pp4};
Tp = {t1(Pp1);t2(Pp2);t3(Pp3);t4(Pp4)};
Vtp = {vt1(Pp1);vt2(Pp2);vt3(Pp3);vt4(Pp4)};
Itp = {it1(Pp1);it2(Pp2);it3(Pp3);it4(Pp4)};
SOCp = {soc1(Pp1);soc2(Pp2);soc3(Pp3);soc4(Pp4)};
Tt = {t1;t2;t3;t4};
Vt = {vt1;vt2;vt3;vt4};
It = {it1;it2;it3;it4};
SOCt = {soc1;soc2;soc3;soc4};
%% Saving data
save("DS_001_treatedData","Pp","Tp","Vtp","Itp","SOCp","Tt","Vt","It",...
     "SOCt")