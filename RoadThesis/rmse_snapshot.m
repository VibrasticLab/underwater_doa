%Calculate the RMSE v.s. number of snapshots
clc
clear all
close all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[data,fs] =audioread('H1_45derajat_1.wav');
cita=resample(data,8000,16000);
% cita = [-54.8 -28.6 -9.2 10.5 31.4 56.7];
% K = length(cita);%number of sources
K=length(cita);
AngRangleft = -90;
AngRangright = 90;
min_deg = -70;
max_deg = 70;
%Parameters setting
deploy = [0 1 4 6];%Deployment of the sensors
M = length(deploy);%number of the sensors
Kmax = M*(M-1)/2;
f = 200;%carrier frequency
c = 340;%propagation velocity
d = 0.5*c/f;%sensor spacing unit
j = sqrt(-1);
SNR = 0;%signal-to-noise ratio
sigman = 0.1;
angleinterv = 1;%interval of angle searching for proposed method
MonteCarlo = 500;
MNit = 100;%maximum number of iteration
rou = 0.0001;%terminate the iteration when norm(gamma_new-gamma_old,2)/norm(gamma_old,2) < rou
%%
snapshot = 50:50:500;
snaplength = length(snapshot);
rmse1 = zeros(1,snaplength);
k0 = 0;
for snap = 50:50:500
    k0 = k0+1;
    mc1 = 0;
    for mc = 1:MonteCarlo
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Array Output Data Generation
        X = MRA_output(cita,deploy,f,c,d,snap,sigman,SNR);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        tic 
        R = X*X'/snap;
        try
           
        [gamma, lamda, alpha0, beta, mu, Sigma, theta_deg] = NNSBL(R, snap, deploy, f, c, d, MNit, rou, angleinterv, AngRangleft, AngRangright);
        N = length(mu);
        sspec1 = zeros(N,1);
        for k = 1:N
            if (-1*mu(k)/sqrt(2*Sigma(k,k))) < 10
                sspec1(k) = Sigma(k,k) + mu(k)^2 + mu(k)*sqrt(2*abs(Sigma(k,k))/pi)*exp(-1*mu(k)^2/(2*abs(Sigma(k,k))))/erfc(-1*mu(k)/sqrt(2*abs(Sigma(k,k))));
            else
                sspec1(k) = abs(Sigma(k,k));
            end
        end
        doa1 = Peaksearch(sspec1, theta_deg, min_deg, max_deg, K);
        rmse1(k0) = rmse1(k0)+(doa1-cita)*(doa1-cita)';
        fprintf('snapshots: %s; MonteCarlo Time: %s; RMSE: %s. ',num2str(snap),num2str(mc),num2str(sqrt((doa1-cita)*(doa1-cita)'/K)));
        
        catch
            fprintf('snapshots: %s; MonteCarlo Time: %s; ERROR! ',num2str(snap),num2str(mc));
        end
        toc
    end
    rmse1(k0) = sqrt(rmse1(k0)/(K*MonteCarlo));
end
%%
%Show results by figures
figure(2)
plot(snapshot,rmse1,'ro-');
xlabel('Number of Snapshots');
ylabel('RMSE (Degree)');
legend('Proposed');
