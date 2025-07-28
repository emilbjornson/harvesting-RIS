%%This Matlab script generates Figure 8 in the paper:
%
%Morteza Tavana, Meysam Masoudi, Emil Björnson, “Energy Harvesting
%Maximization for Reconfigurable Intelligent Surfaces Using Amplitude
%Measurements,” IEEE Transactions on Communications, vol. 72, no. 4, pp.
%2201-2215, April 2024.
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
N = 100; % Number of RIS elements
L_set = [3,10,30,100]; % Considered number of measurement phases per element
iter_Alg = 10; % Number of iterations (M) of the proposed algorithm
iter = 1e3; % Number of channel realizations to compute the mean normalized achieved power (MNAP)
SNR_dB = -30:2:30;
SNR = 10.^(SNR_dB/10);
sigma = sqrt(1./SNR);
y_m = zeros(length(L_set),length(SNR));
y_m_rnd = zeros(length(L_set),length(SNR));
figure;
hold on
for ll = 1:length(L_set)
    L = L_set(ll);
    Phi = (0:L-1)*2*pi/L;
    A = [ones(1,L);cos(Phi);sin(Phi)]';
    b = exp(1j*Phi);
    for l = 1:length(SNR)
        [ll,l]
        y = zeros(iter_Alg+1,N,iter);
        y_rnd = zeros(iter_Alg*L+1,N,iter);
        theta = 2*pi*rand(N,iter);
        theta_rnd = theta;
        Z = 1/sqrt(2)*(randn(N,iter)+1j*randn(N,iter));
        myfun= @(x) abs(sum(Z.*exp(1j*x),1)+sigma(l)/sqrt(2)*(randn(1,iter)+1j*randn(1,iter))).^2;
        y(1,:,:) = repmat(myfun(theta),N,1);
        y_rnd(1,:,:) = y(1,:,:);
        m = sum(abs(Z),1).^2*(1+sigma(l)^2/N);
        %%    Proposed Algorithm
        for j = 1:iter_Alg
            theta_old = theta;
            for i = 1:N
                y_tmp = zeros(L,iter);
                for kkk = 1:L
                    tmp = theta;
                    tmp(i,:) = theta(i,:) + Phi(kkk);
                    y_tmp(kkk,:)=myfun(tmp);
                end
                x = A\y_tmp;
                theta(i,:) = theta(i,:)+atan2(x(3,:),x(2,:));
                y(j+1,i,:) = myfun(theta);
            end
        end
        y_norm = y./ reshape(m,1,1,[]);
        y_norm = permute(y_norm,[2 1 3]);
        y_norm = [reshape(y_norm(1,1,:),1,iter); reshape(y_norm(:,2:end,:),[],iter)];
        y_norm_mat = y_norm;
        %%    Random Algorithm
        for j = 1:(iter_Alg*L)
            theta_old_rnd = theta_rnd;
            for i = 1:N
                theta_tmp_rnd = theta_rnd;
                theta_tmp_rnd(i,:) = 2*pi*rand(iter,1);
                y_tmp_1 = myfun(theta_rnd);
                y_tmp_rnd = myfun(theta_tmp_rnd);
                for kk = 1:iter
                    if y_tmp_rnd(kk) > y_tmp_1(kk)
                        theta_rnd(i,kk) = theta_tmp_rnd(i,kk);
                        y_rnd(j+1,i,kk) = y_tmp_rnd(kk);
                    else
                        y_rnd(j+1,i,kk) = y_tmp_1(kk);
                    end
                end
            end
        end
        y_norm_rnd = y_rnd./reshape(m,1,1,[]);
        y_norm_rnd = permute(y_norm_rnd,[2 1 3]);
        y_norm_rnd = [reshape(y_norm_rnd(1,1,:),1,iter); reshape(y_norm_rnd(:,2:end,:),[],iter)];
        y_norm_mat_rnd =  y_norm_rnd;
        y_proposed_mean = mean(y_norm_mat,2);
        y_proposed_mean_rnd = mean(y_norm_mat_rnd,2);
        y_m(ll,l)= y_proposed_mean(end);
        y_m_rnd(ll,l)= y_proposed_mean_rnd(end);
    end
    if L==1
        plot(SNR_dB,y_m(ll,:),'-b','linewidth',1)
    elseif L==10
        plot(SNR_dB,y_m(ll,:),'--b','linewidth',1)
    elseif L==30
        plot(SNR_dB,y_m(ll,:),':b','linewidth',1)
    else
        plot(SNR_dB,y_m(ll,:),'-.b','linewidth',1)
    end
end
plot(SNR_dB,mean(y_m_rnd,1),'r','linewidth',1)
%% Plot
grid on
box on
xlabel('SNR (dB)')
ylabel('MNAP at \vartheta')
legend('Proposed L=3','Proposed L=10','Proposed L=30','Proposed L=100','Random')