function [gamma, lamda, alpha0, beta, mu, Sigma, theta_deg] = NNSBL(R, snap, deploy, f, c, d, MNit, rou, r, AngRangleft, AngRangright)
% gamma: variance for the signal in each assumed location
% lamda: parameter of the hyprior for gamma
% alpha0: vairiance of the noise
% beta: off-grid distance for each grid point
% N(mu, Sigma): posterior distribution for the signal energy
% R: sample covariance matrix for the array data
% snap: number of snapshots
% deploy: the deployment of the sensors in units of d
% f: carrier frequency
% c: propagation velocity
% d: the unit of inter-sensor distance
% MNit: maximum number of iteration
% rou: terminate the iteration when norm(gamma_new-gamma_old,2)/norm(gamma_old,2) < rou
% r: interval of grid points
%%
%%%%%%%%%%% Initialization %%%%%%%%%%%%
pa = 10e-9;
pb = 10e-9;
j = sqrt(-1);
M0 = size(R,1);
M = M0^2;
I = eye(M0);
y = zeros(M,1);%data vector
Im = zeros(M,1);%Vectorize an identity matrix
for k = 1:M0
    y(((k-1)*M0+1):k*M0) = R(:,k);%vectorization
    Im(((k-1)*M0+1):k*M0) = I(:,k);
end
Rw = kron(R.',R)/snap;%covariance of the residual
%--------- A & B---------%
theta = ((AngRangleft+r):r:(AngRangright-r))*pi/180;
theta_deg = (AngRangleft+r):r:(AngRangright-r);
N = length(theta);
A0 = zeros(M0,N);%the basic array manifold matrix
Ab = zeros(M,N);
Bb = zeros(M,N);
dp = zeros(M,1);
for k = 1:M0
    dp(((k-1)*M0+1):k*M0) = -1*deploy(k)*ones(M0,1) + deploy';
end
for k = 1:N
    a = exp(-j*2*pi*(f/c)*d*sin(theta(k))*deploy');
    A0(:,k) = a(:);
    a0 = kron(conj(a),a);
    b0 = -j*2*pi*(f/c)*d*(pi/180)*cos(theta(k))*(dp.*a0);
    Ab(:,k) = a0(:);
    Bb(:,k) = b0(:);
end
Aw = [real(Ab);imag(Ab)];
Bw = [real(Bb);imag(Bb)];
%%%%% y = (A + B * diag(beta)) * w + alpha0 * Im + nb, where nb ~ CN(0,Rw)
%---------- gamma ----------%
gamma = ones(N,1) * 10e-12;
%---------- lamda ----------%
lamda = (N-1+pa)/(0.5*ones(1,N)*gamma+pb);
%---------- alpha0 ---------%
alpha0 = 10e-6;
%---------- beta -----------%
beta = zeros(N,1);
%-------------------------------------------------------------------------
Rb = (1/2)*[real(Rw) -1*imag(Rw);imag(Rw) real(Rw)];
Rwv = pinv(Rw);
Rv = 2*[real(Rwv) -1*imag(Rwv);imag(Rwv) real(Rwv)];% Rv = inv(Rb)
FlagIter = 1;
t = 0;
Fi = [real(Ab);imag(Ab)];
Wmu = Fi.'*Rv;
Wa1 = Im.'*real(Rwv)*real(y) - Im.'*imag(Rwv)*imag(y);
Wa2 = Im.'*real(Rwv)*real(Ab) - Im.'*imag(Rwv)*imag(Ab);
Wa3 = Im.'*real(Rwv)*Im;
%%
%Begin iteration
while (FlagIter)
    gamma_old = gamma;
    lamda_old = lamda;
    alpha0_old = alpha0;
    
    yb = [real(y)-alpha0_old*Im;imag(y)];
    %Calculate the posterior mean and covariance
    Gamma = diag(gamma_old);
    if (2*M)<N
        Sigma = Gamma-Gamma*Fi.'*pinv(Rb+Fi*Gamma*Fi.')*Fi*Gamma;%posterior covariance matrix
    else
        Ginv = eye(N);
        for k = 1:N
            Ginv(k,k) = 1/gamma_old(k);
        end
        Sigma = pinv(Fi.'*Rv*Fi+Ginv);
    end
    mu = Sigma*Wmu*yb;%posterior mean vector
    nnmean = zeros(N,1);%mean
    nnmeansq = zeros(N,1);%square mean
    
    for k = 1:N
        if (-1*mu(k)/sqrt(2*Sigma(k,k))) < 10
            nnmean(k) = mu(k) + sqrt(2*Sigma(k,k)/pi)*exp(-1*mu(k)^2/(2*abs(Sigma(k,k))))/erfc(-1*mu(k)/sqrt(2*abs(Sigma(k,k))));
            nnmeansq(k) = abs(Sigma(k,k)) + mu(k)^2 + mu(k)*sqrt(2*abs(Sigma(k,k))/pi)*exp(-1*mu(k)^2/(2*abs(Sigma(k,k))))/erfc(-1*mu(k)/sqrt(2*abs(Sigma(k,k))));
        else
            nnmean(k) = 0;
            nnmeansq(k) = abs(Sigma(k,k));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %gamma
    for k = 1:N
        gamma(k) = -1/(2*lamda_old) + sqrt(1/(4*lamda_old^2)+nnmeansq(k)/lamda_old);
    end
    
    %lamda
    lamda = (N-1+pa)/(0.5*ones(1,N)*gamma+pb);
    
    %alpha0
    alpha0 = (Wa1-Wa2*nnmean)/Wa3;
    
    if alpha0 <= 0
        alpha0 = alpha0_old;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = t+1;
    verr = norm(gamma-gamma_old,2)/norm(gamma_old,2);
    if (verr < rou)||(t > MNit)
        FlagIter = 0;
    end
end
%%
%%%% Calculate the angle offsets: beta
ni = peak_find(nnmeansq);
lenpk = length(ni);
theta_est = theta(ni(:));
nnmean_est = zeros(lenpk,1);
nnmean_est(:) = nnmean(ni(:));
nnmeansq_est = zeros(lenpk,1);
nnmeansq_est(:) = nnmeansq(ni(:));
r_cur = r;
t_recur = 5;
for t = 1:t_recur
    Abr = zeros(M,lenpk);
    Bbr = zeros(M,lenpk);
    for k = 1:lenpk
        a = exp(-j*2*pi*(f/c)*d*sin(theta_est(k))*deploy');
        a0 = kron(conj(a),a);
        b0 = -j*2*pi*(f/c)*d*(pi/180)*cos(theta_est(k))*(dp.*a0);
        Abr(:,k) = a0(:);
        Bbr(:,k) = b0(:);
    end
    Awr = [real(Abr);imag(Abr)];
    Bwr = [real(Bbr);imag(Bbr)];
    
    W = diag(nnmeansq_est);
    yb = [real(y)-alpha0*Im;imag(y)];
    P = (Bwr.'*Rv*Bwr).*W;
    V1 = (Bwr.'*Rv.').*(nnmean_est*yb.');
    V2 = (Bwr.'*Rv*Awr).*W;
    v = V1*ones(2*M,1)-V2*ones(lenpk,1);
    for k = 1:lenpk
        bk = v(k)/P(k,k);
        if bk < -r_cur/2
            beta(ni(k)) = beta(ni(k))-r_cur/2;
            theta_est(k) = theta_est(k)-r_cur/2;
        elseif bk > r_cur/2
            beta(ni(k)) = beta(ni(k))+r_cur/2;
            theta_est(k) = theta_est(k)+r_cur/2;
        else
            beta(ni(k)) = beta(ni(k))+bk;
            theta_est(k) = theta_est(k)+bk;
        end
    end
    r_cur = r_cur/2;
end
%Calculate the posterior mean and covariance
Gamma = diag(gamma);
if (2*M)<N
    Sigma = Gamma-Gamma*Fi.'*pinv(Rb+Fi*Gamma*Fi.')*Fi*Gamma;%posterior covariance matrix
else
    Ginv = eye(N);
    for k = 1:N
        Ginv(k,k) = 1/gamma(k);
    end
    Sigma = pinv(Fi.'*Rv*Fi+Ginv);
end
mu = Sigma*Fi.'*Rv*yb;%posterior mean vector
theta_deg = theta_deg + beta.';
