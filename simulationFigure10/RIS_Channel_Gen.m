clear all
close all
clc
N = unique(floor(sqrt(logspace(0,7,80))));
for i = 1:length(N)
    n = N(i);
    try
        load(['Channels/','z_N=',num2str(n.^2),'.mat'])
    catch
        RIS = RIS_Channel(n,n);
        RIS.RIS_CH();
        z = RIS.z;
        save(['Channels/','z','_N=',num2str(RIS.N_RIS_row*RIS.N_RIS_col)],"z")
    end
end