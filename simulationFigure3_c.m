%%This Matlab script generates Figure 3c in the paper:
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
z0 = 1;
absZ = 10; % The result holds of absZ = 1/10
phiZ = 0:2*pi/60:2*pi;
phiZ = phiZ(1:end-1);
iter = 1e4;
L=3;
SNR_dB = [-10,0,10,20];
SNR = 10.^(SNR_dB/10);
Phi = ((0:L-1)*2*pi/L)';
A = [ones(L,1),cos(Phi),sin(Phi)];
A_inv = inv(A'*A)*A';
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    "OptimalityTolerance",1e-10,"FunctionTolerance",1e-10,"StepTolerance",1e-14,"MaxIterations",1e5,"MaxFunctionEvaluations",1e5);
options2 = optimoptions('fmincon','Display', 'off');
%%
MSE_Exact_ML_mean = zeros(length(SNR_dB),length(absZ),length(phiZ));
MSE_Linear_mean = zeros(length(SNR_dB),length(absZ),length(phiZ));
MSE_theory_Linear = zeros(length(SNR_dB),length(absZ),length(phiZ));
MSE_theta_Linear = zeros(length(SNR_dB),length(absZ),length(phiZ));
MSE_theta_Exact_ML = zeros(length(SNR_dB),length(absZ),length(phiZ));
Bias_theta_Linear = zeros(length(SNR_dB),length(absZ),length(phiZ));
Bias_theta_Exact_ML = zeros(length(SNR_dB),length(absZ),length(phiZ));

for kk = 1:length(SNR_dB)
    for jj = 1:length(absZ)
        parfor mm = 1:length(phiZ)
            [kk, jj, mm]
            MSE_Linear = zeros(iter,1);
            MSE_Exact_ML = zeros(iter,1);
            theta_hat_Linear = zeros(iter,1);
            theta_hat_Exact_ML = zeros(iter,1);
            exitflag = zeros(iter,1);
            z = absZ(jj)*exp(1j*phiZ(mm));
            sigma = sqrt((1+abs(z).^2)./2./SNR(kk));
            x = [abs(z0).^2+abs(z).^2, 2*real(z0*conj(z)), 2*imag(z0*conj(z))]';
            theta= atan2(x(3),x(2))*ones(iter,1);
            for i = 1:iter
                  %[kk, jj, mm, i]
                phi_iter = Phi;
                n = sigma/sqrt(2)*(randn(L,1)+1j*randn(L,1));
                y = abs(z0+z*exp(1j*phi_iter)+n).^2;
                %% Linear Estimator
                x_hat_Linear = A_inv*y;
                MSE_Linear(i) = sum((x_hat_Linear-x).^2,1);
                theta_hat_Linear(i) = atan2(x_hat_Linear(3),x_hat_Linear(2));
                %% Exact ML Estimator
                if SNR_dB(kk)<21
                    try
                        f0 = @(w) sum(A*w./sigma^2-log(besseli(0,2/sigma^2*sqrt((A*w).*y),1))-2/sigma^2*sqrt((A*w).*y));
                        w0 = x;
                        [x_hat_exact_ML,~,exitflag(i)] = fmincon(f0,w0,[],[],[],[],[],[],@f1,options2);
                    catch
                        f0 = @(w) double(vpa(sum(sym(A)*sym(w)./sym(sigma)^sym(2)-log(besseli(0,sym(2)/sym(sigma)^sym(2)*sqrt(sym(A)*sym(w)).*sym(y))))));
                        w0 = x;
                        [x_hat_exact_ML,~,exitflag(i)] = fmincon(f0,w0,[],[],[],[],[],[],@f1,options);
                    end
                else
                    f0 = @(w) sum((sqrt(y)-sqrt(A*w)).^2+sigma^2/4*log(A*w));
                    w0 = x;
                    c = 0.25*sigma^4./y;
                    [x_hat_exact_ML,~,exitflag(i)] = fmincon(f0,w0,-A,-c,[],[],[],[],@f1,options2);
                end
                MSE_Exact_ML(i) = sum((x_hat_exact_ML-x).^2,1);
                theta_hat_Exact_ML(i) = atan2(x_hat_exact_ML(3),x_hat_exact_ML(2));
            end
            %% Bias of theta
            Bias_theta_Linear(kk,jj,mm) = mean(atan2(sin(theta_hat_Linear-theta),cos(theta_hat_Linear-theta)));
            Bias_theta_Exact_ML(kk,jj,mm) = mean(atan2(sin(theta_hat_Exact_ML(exitflag>0)-theta(1)),cos(theta_hat_Exact_ML(exitflag>0)-theta(1))));
            %%
            %% MSE of x
            MSE_Linear_mean(kk,jj,mm) = sqrt(mean(MSE_Linear));
            MSE_Exact_ML_mean(kk,jj,mm) = sqrt(mean(MSE_Exact_ML(exitflag>0)));
            MSE_theory_Linear(kk,jj,mm) = 2/3*2*sigma^2*trace(A_inv'*A_inv)+sigma^4*trace(A_inv'*A_inv)+sigma^4;
            %% MSE of theta
            tmp = mod(theta_hat_Linear-theta,2*pi);
            tmp2 = min(tmp,2*pi-tmp);
            MSE_theta_Linear(kk,jj,mm) = sqrt(mean(tmp2.^2));
            tmp = mod(theta_hat_Exact_ML-theta,2*pi);
            tmp2 = min(tmp,2*pi-tmp);
            MSE_theta_Exact_ML(kk,jj,mm) = sqrt(mean(tmp2.^2));
        end
    end
end
%%
c = {'--b', '--r', '--m','--k'};
c2 = {'-b', '-r', '-m','-k'};
%% Plot MSE theta
F1 = figure;
for kk = 1:length(SNR_dB)
    tmpLinear = reshape(MSE_theta_Linear(kk,1,:),[],1);
    tmpML = reshape(MSE_theta_Exact_ML(kk,1,:),[],1);
    tmpLinear = [tmpLinear;tmpLinear(1)];
    tmpLinear = tmpLinear(end:-1:1);
    tmpML = [tmpML;tmpML(1)];
    tmpML = tmpML(end:-1:1);
    plot([phiZ,2*pi],tmpLinear,c{kk},[phiZ,2*pi],tmpML,c2{kk})
    hold on
end
grid on
box on
xlabel('\vartheta (rad)')
ylabel('RMSE of \vartheta')
%% Functions
function [c,ceq] = f1(x)
c = sqrt(x(2)^2+x(3)^2)-x(1);
ceq = [];
end



