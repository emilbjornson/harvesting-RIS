%%This Matlab script generates Figure 9 in the paper:
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
%%
N = 10; % Number of RIS elements
M = 4; % Number of possible discrete phases
Omega = (0:M-1)*2*pi/M; % set of possible discrete phases
Phi = Omega(1:3); % Set of the measurement phases
A = [1,1,1;cos(Phi(1)),cos(Phi(2)),cos(Phi(3));...
    sin(Phi(1)),sin(Phi(2)),sin(Phi(3))]';
Possible_Theta = (fullfact(M*ones(1,N)));
iter_Alg = 1+100; % Number of the iterations (M) of the proposed algorithm
iter = 1e3; % Number of channel realizations to compute the mean normalized achieved power (MNAP)
m_BF = zeros(iter,1);
y_norm_mat = zeros(iter_Alg*N,iter);
y_norm_mat_rnd = zeros(iter_Alg*N,iter);
for kk = 1:iter
    kk
    y = zeros(iter_Alg+1,N);
    y_rnd = zeros(iter_Alg+1,N);
    theta = ones(N,1)*Omega(1);
    theta_rnd = ones(N,1)*Omega(1);
    Z = (randn(N,1)+1j*randn(N,1));
    myfun= @(x) abs(sum(Z.*exp(1j*x(:)))).^2;
    y(1,:) = myfun(theta);
    y_rnd(1,:) = myfun(theta_rnd);
    m = sum(abs(Z))^2;
    %% Proposed Algorithm
    for j = 1:iter_Alg
        theta_old = theta;
        for i = 1:N
            y_tmp = zeros(3,1);
            theta_1 = theta;
            theta_2 = theta;
            theta_3 = theta;
            theta_1(i) = theta(i) + Phi(1);
            theta_2(i) = theta(i) + Phi(2);
            theta_3(i) = theta(i) + Phi(3);
            y_tmp(1) = myfun(theta_1);
            y_tmp(2) = myfun(theta_2);
            y_tmp(3) = myfun(theta_3);
            x = A\y_tmp;
            alpha = atan2(x(3),x(2))+theta(i);
            tmp = mod(Omega-alpha,2*pi);
            tmp = min(tmp,2*pi-tmp);
            [~,Index_tmp] = min(tmp);
            theta(i) = Omega(Index_tmp);
            y(j+1,i) = myfun(theta);
        end
        if theta==theta_old
            y = y(1:j+1,:);
            break;
        end
    end
    y_norm = y(:,:)./m;
    y_norm = [y_norm(1,1); reshape(y_norm(2:end,:)',[],1)];
    y_norm_mat(1:length(y_norm),kk) =  y_norm;
    if length(y_norm)<iter_Alg
        y_norm_mat(length(y_norm)+1:end,kk) =  y_norm(end);
    end
    %% Random Algorithm
    for j = 1:iter_Alg
        theta_old_rnd = theta_rnd;
        for i = 1:N
            %y_tmp_rnd = myfun(theta_rnd);
            theta_tmp_rnd = theta_rnd;
            while true
                tmp_phase = Omega(randi([1,M],1));
                if tmp_phase ~=theta_rnd(i)
                    break;
                end
            end
            theta_tmp_rnd(i) = tmp_phase;%Omega(randi([1,M],1));
            y_tmp_1 = myfun(theta_rnd);
            y_tmp_rnd = myfun(theta_tmp_rnd);
            if y_tmp_rnd > y_tmp_1
                theta_rnd(i) = theta_tmp_rnd(i);
                y_rnd(j+1,i) = y_tmp_rnd;
            else
                y_rnd(j+1,i) = y_tmp_1;
            end
        end
        if theta_rnd==theta_old_rnd
            y_rnd = y_rnd(1:j+1,:);
            break;
        end
    end
    y_norm_rnd = y_rnd(:,:)./m;
    y_norm_rnd = [y_norm_rnd(1,1); reshape(y_norm_rnd(2:end,:)',[],1)];
    y_norm_mat_rnd(1:length(y_norm_rnd),kk) =  y_norm_rnd;
    if length(y_norm_rnd)<iter_Alg*N
        y_norm_mat_rnd(length(y_norm_rnd)+1:end,kk) =  y_norm_rnd(end);
    end
    %% Brute Force Method
    y_BF = zeros(size(Possible_Theta,1),1);
    fun_offline = repmat(Z,1,M).*repmat(exp(1j*Omega),N,1);
    fun_offline = fun_offline(:);
    y_BF= abs(sum(fun_offline((Possible_Theta(1:end,:)'-1)*N+(1:N)'))).^2;
    [m_BF(kk),Index_BF] = max(y_BF);
    m_BF(kk) = m_BF(kk)/m;
end
y_proposed_mean = mean(y_norm_mat,2);
y_proposed_mean_rnd = mean(y_norm_mat_rnd,2);
y_optimal = repmat(mean(m_BF),1,length(y_proposed_mean));
%% Plot
step = 1;
MAX = 200;
% Confidence interval
c1 = qfuncinv(2.5e-2)*std(m_BF)/sqrt(iter)*ones(iter_Alg*N,1);
c2 = qfuncinv(2.5e-2)/sqrt(iter)*ones(iter_Alg*N,1);
c3 = qfuncinv(2.5e-2)/sqrt(iter)*ones(iter_Alg*N,1);
for i = 1:(iter_Alg*N)
    c2(i) = c2(i)*std(y_norm_mat(i,:));
    c3(i) = c3(i)*std(y_norm_mat_rnd(i,:));
end
F1 = figure;
hold on
errorbar((0:step:(MAX-1))*3,y_optimal(1:step:MAX),c1(1:step:MAX),'k','linewidth',1)
errorbar((0:step:(MAX-1))*3,y_proposed_mean(1:step:MAX),c2(1:step:MAX),'b','linewidth',1)
errorbar((0:step:(MAX-1)),y_proposed_mean_rnd(1:step:MAX),c3(1:step:MAX),'r','linewidth',1)
grid on
box on
xlabel('Number of measurements')
ylabel('MNAP')
legend('Maximum feasible','Proposed method','Random method')