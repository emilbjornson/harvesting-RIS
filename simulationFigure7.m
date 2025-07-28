%%This Matlab script generates Figure 7 in the paper:
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
N = 100; % Number of the RIS elements
iter = 51; % Number of iterations (M) of the algorithm
seed = 1e3; % Number of channel realizations

y_opt = zeros(iter+1,seed);
y_opt_norm = zeros(iter+1,seed);
y_opt_rnd = zeros(iter+1,seed);
y_opt_rnd_norm = zeros(iter+1,seed);

m = zeros(seed,1);
for kk = 1:seed
    kk
    x_opt_rnd = zeros(N,1);
    x_opt = zeros(N,1);
    x_tmp_pi_2 = zeros(N,1);
    x_tmp_pi = zeros(N,1);
    Z = (randn(N,1)+1j*randn(N,1));
    [y_opt(1,kk),m(kk)] = myfun(x_opt,Z);
    [y_opt_rnd(1,kk),m(kk)] = myfun(x_opt_rnd,Z);
    for j = 1:iter
        for i = 1:N
            %% Proposed Algorithm
            x_tmp_pi_2 = x_opt;
            x_tmp_pi = x_opt;
            x_tmp_pi_2(i) = x_opt(i) + pi/2;
            x_tmp_pi(i) = x_opt(i) + pi;
            f_tmp = myfun(x_opt,Z);
            f_tmp_pi_2 = myfun(x_tmp_pi_2,Z);
            f_tmp_pi = myfun(x_tmp_pi,Z);
            phi = atan2(2*f_tmp_pi_2-f_tmp-f_tmp_pi,f_tmp-f_tmp_pi);
            x_opt(i) = x_opt(i)+phi;
            %% Random Algorithm
            f_opt_rnd = myfun(x_opt_rnd,Z);
            x_tmp_rnd = x_opt_rnd;
            x_tmp_rnd(i) = rand(1)*2*pi;
            f_tmp_rnd = myfun(x_tmp_rnd,Z);
            if f_tmp_rnd>f_opt_rnd
                x_opt_rnd(i) = x_tmp_rnd(i);
            end
        end
        y_opt(j+1,kk) = myfun(x_opt,Z);
        y_opt_rnd(j+1,kk) = myfun(x_opt_rnd,Z);
    end
    y_opt_norm(:,kk) = y_opt(:,kk)./m(kk);
    y_opt_rnd_norm(:,kk) = y_opt_rnd(:,kk)./m(kk);
end

%%
figure;
cdfplot(y_opt_norm(2,:));
hold on
cdfplot(y_opt_rnd_norm(11,:));
cdfplot(y_opt_rnd_norm(51,:));
legend('Proposed method: Iter. 1', 'Random method: Iter. 10', 'Random method: Iter. 50')
xlabel('NAP')
ylabel('CDF')
grid on
box on
%%
function  [y,m] = myfun(x,Z)
m = sum(abs(Z)).^2;
y = abs(sum(Z.*exp(1j*x(:)))).^2;
end