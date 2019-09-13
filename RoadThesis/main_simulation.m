%%%Codes for "Source Localization for Sparse Array using Nonnegative Sparse Bayesian Learning": Spatial Spectrum
clear
close all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[data,fs] =audioread('H1_45derajat_1.wav');
cita=resample(data,8000,16000);
%cita = [-54.8 -28.6 -9.2 10.5 31.4 56.7];
K = length(cita);%number of sources
K = length(cita);%number of sources
AngRangleft = -90;
AngRangright = 90;
min_deg = -70;
max_deg = 70;
%Parameters setting
deploy = [0 1 4 6];%Deployment of the sensors
M = length(deploy);%number of the sensors
f = 200;%carrier frequency
c = 340;%propagation velocity
d = 0.5*c/f;%sensor spacing unit
j = sqrt(-1);
snap = 200;%number of snapshots
SNR = -5;%signal-to-noise ratio
sigman = 0.1;
angleinterv = 1;%interval of angle searching
MNit = 100;%maximum number of iteration
rou = 0.0001;%terminate the iteration when norm(gamma_new-gamma_old,2)/norm(gamma_old,2) < rou
%Signal generation
X = MRA_output(cita,deploy,f,c,d,snap,sigman,SNR);
%%
Kmax = M*(M-1)/2;
R = X*X'/snap;
%Proposed (NNSBL)
tic
[gamma, lamda, alpha0, beta, mu, Sigma, theta_deg1] = NNSBL(R, snap, deploy, f, c, d, MNit, rou, angleinterv, AngRangleft, AngRangright);
N = length(mu);
sspec1 = zeros(N,1);
for k = 1:N
    if (-1*mu(k)/sqrt(2*Sigma(k,k))) < 10
        sspec1(k) = Sigma(k,k) + mu(k)^2 + mu(k)*sqrt(2*abs(Sigma(k,k))/pi)*exp(-1*mu(k)^2/(2*abs(Sigma(k,k))))/erfc(-1*mu(k)/sqrt(2*abs(Sigma(k,k))));
    else
        sspec1(k) = abs(Sigma(k,k));
    end
end
sspec1 = sspec1/max(sspec1);
doa1 = Peaksearch(sspec1, theta_deg1, min_deg, max_deg, K);

rmse1 = sqrt((doa1-cita)*(doa1-cita).'/K);
fprintf('RMSE (NNSBL): %s. ',num2str(rmse1));
toc
%Conventional SBL
tic
[gamma, lamda, alpha0, beta, mu, Sigma, theta_deg2] = Conven_SBL(R, snap, deploy, f, c, d, MNit, rou, angleinterv, AngRangleft, AngRangright);
N = length(mu);
sspec2 = zeros(1,N);
for k = 1:N
    if (-1*mu(k)/sqrt(2*Sigma(k,k))) < 10
        sspec2(k) = Sigma(k,k) + mu(k)^2 + mu(k)*sqrt(2*abs(Sigma(k,k))/pi)*exp(-1*mu(k)^2/(2*abs(Sigma(k,k))))/erfc(-1*mu(k)/sqrt(2*abs(Sigma(k,k))));
    else
        sspec2(k) = abs(Sigma(k,k));
    end
end
sspec2 = sspec2/max(sspec2);
doa2 = Peaksearch(sspec2, theta_deg2, min_deg, max_deg, K);
rmse2 = sqrt((doa2-cita)*(doa2-cita).'/K);
fprintf('RMSE (Conventional SBL): %s. ',num2str(rmse2));
toc
plot(theta_deg1,sspec1,'k-','LineWidth',1);
hold on
plot(theta_deg2,sspec2,'k--','LineWidth',1);
for k = 1:K
    hold on
    xv = cita(k)*ones(1,1001);
    yv = 0:0.001:1;
    plot(xv,yv,'k:');
end
hold off
figure(3)
axis([-90 90 0 1]);
set(gca,'XTick',-90:10:90)  
set(gca,'XTickLabel',{'-90','-80','-70','-60','-50','-40','-30','-20','-10','0','10','20','30','40','50','60','70','80','90'})
xlabel('Bearing (Deg)');
ylabel('Normalized Spatial Spectra');
legend('Proposed (NNSBL)','Conventional SBL')
