%%This Matlab script generates Figure 10 in the paper:
%
%Morteza Tavana, Meysam Masoudi, Emil Björnson, “Energy Harvesting
%Maximization for Reconfigurable Intelligent Surfaces Using Amplitude
%Measurements,” IEEE Transactions on Communications, vol. 72, no. 4, pp.
%2201-2215, April 2024.
%
%Note: You need to run RIS_Channel_Gen.m to create and store the channels
%before running this script.
%
%Download article: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10356096
%
%This is version 1.0 (Last edited: 2024-04-15)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.



clear all
close all
clc
L = 3;
P_TX = 1; % Watt
Phi = (0:L-1)*2*pi/L;
A = [ones(1,L);cos(Phi);sin(Phi)]';
b = exp(1j*Phi);
iter_Alg = 3; % Number of the iterations (M) of the proposed algorithm
iter_Alg_RND = 10; % Number of the iterations (M) of the random algorithm
SNR_dB = [-20,-10,0,10];
SNR = 10.^(SNR_dB/10);
N_set = unique(floor(sqrt(logspace(0,7,80)))).^2;
RX_power_converged_proposed = zeros(length(SNR),length(N_set));
RX_power_converged_random = zeros(length(SNR),length(N_set));
Harvested_power_converged_proposed = zeros(length(SNR),length(N_set));
Harvested_power_converged_random = zeros(length(SNR),length(N_set));
RIS = RIS_Channel(1,1);
for l = 1:length(N_set)
    l
    N = N_set(l);
    load(['Channels/','z_N=',num2str(N),'.mat'])
    z = z.';
    for ll = 1:length(SNR)
        if N<=100
            iter = 5000;
        elseif N<2000000
            if ll<3
                iter = max(10,floor(1000000/N));
            else
                iter =  max(1,floor(1000000/N));
            end
        else
            iter = 1;
        end
        %iter = 1;%floor(iter/100);
        Z=repmat(z(:),1,iter);
        y = zeros(N,iter_Alg+1,iter);
        y_rnd = zeros(N,iter_Alg_RND*L+1,iter);
        theta = 2*pi*rand(N,iter);
        theta_rnd = theta;
        sigma = sqrt(P_TX*sum(abs(Z).^2,1)/N/SNR(ll));
        myfun= @(x) P_TX*abs(z*exp(1j*x)+sigma./sqrt(2).*(randn(1,iter)+1j*randn(1,iter))).^2;
        myfun_ZN= @(x) P_TX*abs(sum(Z.*exp(1j*x),1)).^2;
        y(:,1,:) = repmat(myfun(theta),N,1);
        y_rnd(:,1,:) = y(:,1,:);
        %%    Proposed Algorithm
        tmp1 = z*exp(1j*theta);
        invA = inv(A);
        ePhi = (exp(1j*(Phi(:)))-1);
        for j = 1:iter_Alg
            noise = sigma./sqrt(2).*(randn(N,iter)+1j*randn(N,iter));
            etheta = exp(1j*(theta));
            for i = 1:N
                y_tmp = P_TX*abs(tmp1+z(i)*(ePhi*etheta(i,:))+sigma./sqrt(2).*(randn(L,iter)+1j*randn(L,iter))).^2;
                x = invA*y_tmp;
                tmp1 = tmp1+z(i)*(etheta(i,:).*(exp(1j*(atan2(x(3,:),x(2,:))))-1));
                theta(i,:) = theta(i,:)+atan2(x(3,:),x(2,:));
                y(i,j+1,:) = P_TX*abs(tmp1+noise(i,:)).^2;
            end
        end
        %%    Random Algorithm
        tmp1 = z*exp(1j*theta_rnd);
        for j = 1:(iter_Alg_RND*L)
            tmp_rnd = 2*pi*rand(N,iter);
            noise_1 = sigma./sqrt(2).*(randn(N,iter)+1j*randn(N,iter));
            noise_2 = sigma./sqrt(2).*(randn(N,iter)+1j*randn(N,iter));
            for i = 1:N
                tmp2 = tmp1+z(i)*(exp(1j*(tmp_rnd(i,:)))-exp(1j*(theta_rnd(i,:))));
                y_tmp_1 = P_TX*abs(tmp1+noise_1(i,:)).^2;
                y_tmp_rnd = P_TX*abs(tmp2+noise_2(i,:)).^2;
                for kk = 1:iter
                    if y_tmp_rnd(kk) > y_tmp_1(kk)
                        theta_rnd(i,kk) = tmp_rnd(i,kk);
                        y_rnd(i,j+1,kk) = y_tmp_rnd(kk);
                        tmp1(kk)=tmp2(kk);
                    else
                        y_rnd(i,j+1,kk) = y_tmp_1(kk);
                    end
                end
            end
        end
        Harvested_power_converged_proposed(ll,l)= mean(RIS.Harvesting_efficiency(myfun_ZN(theta),30,0.07,0.1));
        Harvested_power_converged_random(ll,l)= mean(RIS.Harvesting_efficiency(myfun_ZN(theta_rnd),30,0.07,0.1));
        RX_power_converged_proposed(ll,l)= mean(myfun_ZN(theta));
        RX_power_converged_random(ll,l)= mean(myfun_ZN(theta_rnd));
    end
end


%% Original Code

NN=[1:1:10, 12:2:100,120:20:1000, 1200:200:10000, 12000:2000:100000];
RX_Power_Approx = zeros(length(NN),1);
Inc_Power_Approx = zeros(length(NN),1);
Harvested_Power = zeros(length(NN),1);
for i = 1:length(NN)
    n = NN(i);
    d1 = RIS.d_RIS*RIS.lambda*n/2;
    [RX_Power_Approx(i), Inc_Power_Approx(i)] = RIS.Infinit_surface(d1);
    Harvested_Power(i) = RIS.Harvesting_efficiency(RX_Power_Approx(i),30,0.07,0.1);
end

semilogx(NN.^2, 30+10*log10(Harvested_Power),'-k')
xlabel('Number of RIS elements')
ylabel('Mean achieved power (dBm)')
hold on
grid on
box on
for ll = 1:length(SNR)
    semilogx(N_set,30+10*log10(Harvested_power_converged_proposed(ll,:)),'--b',N_set,30+10*log10(Harvested_power_converged_random(ll,:)),'-.r')
end
legend('Genie-aided','Proposed method','Random method')
