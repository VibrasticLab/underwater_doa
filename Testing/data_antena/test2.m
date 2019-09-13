%Appendix C: m-file that simulates reconstruction of K targets. 
%m-file to research the possibilities of using compressive sensing (with 
%the MFOCUSS algorithm on a line array of 64 elements with 16 random  
%working elements. 
%clear the workspace and close the figures 
clc; close all; clear all;  
nr=0;  %simulation number to show difference when noise is different.  
%       (this number sets the random generator to the specified start value)   

%PARAMETER DEFINITION: CHANGE THESE TO YOUR WHISHES %The array parameters
[data,fc]=audioread('H1_45derajat_1.wav');
N=data %Total number of elements the dense ULA 
M=16; %Number of random working elements that remain in the sparse array 
K=1;  %Number of targets received by the array 
Nsamp = 1; %number of time snapshots/MMV used 
max_scan_offbroadside=59; %define the maximum off boresight scan angle 

%Define environement parameters 
c = 1500;% physconst('LightSpeed'); 
%fc = 1000;              % Operating frequency 3 GHz 
lambda = c/fc;         % Wavelength  
d=lambda/2;            % distance between elements   

%Define the discretization grid 
scan_angle_min=45; %the left edge of the discretization grid 
scan_angle_max=135; %the right edge of the discretization grid 
coverage=2*max_scan_offbroadside; 
Ns=181; %no of discretization grid points in the search area see table below for examples:
%assumed a discretization grid from -90 to 90 degrees setting Ns to:
% 181 results in a grid spacing of 1 degree 
% 361 results in a grid spacing of 0.5 degrees 
% 721 results in a grid spacing of 0.25 degrees 
% 1801 results in a grid spacing of 0.1 degrees 
%etc   

%received signal parameters 
Noise=-20; %received noise level in dB 
Target_power=1;% power level of the targets, choose 1 to normalize to 0dB,  
% if random strength targets, this defines the maximum power level, rember 
% to (un)comment line 69 for the use of random strength targets 
rng(nr); 
%Define the DOA of the target, give up K values between +max_scan_offbroadside! 
DOAangles = [45 90]; 
%or (un)comment the line below to choose K random ON GRID DOAs 
%DOAangles = randperm(coverage,K)-max_scan_offbroadside 
%use randperm to avoid double DOAs 
% or (un)comment the 2 lines below to choose K random OFF GRID DOAs 
% off_grid_distance=0.5 
%specify the off grid distance for the random DOAs 
% DOAangles = randperm(coverage,K)-max_scan_offbroadside+off_grid_distance 
%use randperm to avoid double DOAs   

%Define the MFOCUSS input parameters 
lambda_m=db2pow(Noise); %regularization parameter, according to the  
%description in the noisy case generally close to the noise variance. But 
%also user defined values can be chosen
p=0.8; %Parameter that trades off between the speed of convergence and the  
%       sparsity in the solution. p=0.8 generally leads to good results. 
r=1000; %maximum number of iterations MFOCUSS can use to find a solution. 

% %create a conventional ULA with 64 elements using the phased array toobox 
ULA64 = phased.ULA('NumElements',N,'ElementSpacing',d); 
ULA64.Element.BackBaffled = true; % this makes that the array only radiates to one side   

% For K random DOA angels: 
angs = DOAangles;   

%% create a signal received by the N element ULA 
pos64 = getElementPosition(ULA64); 
nPower = db2pow(Noise);  %noise power at the element in Watts, Power = 1; % power of the incoming signals 
sPower = Target_power; % power of the incoming signals 
%sPower=Target_power*0.1*randi([1,10],1,K); %(un)comment to use random strength signals. 
% generate a multichannel signal received by the ula 
rng('default'); 
rs = rng(0); 
x64 = sensorsig(pos64/lambda,Nsamp,angs,nPower,sPower);   

%% rewrite in order to use CS 
%rewrite the model to X=AS+W=PSI*Z+W and solve Y=PHI*X=PHI*PSI*Z+PHI*W 
X=x64.'; %this creates X; a N x T matrix. 

%% Create Phi %Create an array with only 16 working elements, but let the first and last 
%one always be working setting goodrs to 11 gives the array used throughout the %thesis 
goodrs=11; 
rng(goodrs); 
% first define Phi as a NxN identity matrix from which random rows can be % removed:
Phi = eye(N); %NxN identity matrix 

%% (un)comment this section this when element 1 and N may also be removed: 
% rng(0); 
% elemToRemove = randperm(N);  
%create a row vector with random values from 1 to the total number of elements 
% elemToRemove = elemToRemove(1:N-M); 
%Take the fist 28 values to be the random failing elements. 
% Phi(elemToRemove,:)=[]; 
%"randomly" make te rows of the removed elements zero % 
%take the measurements using the sensing matrix and the signal, this % 
%resluts in a MxN matrix PHI which contains only values of the non broken %
%elements. using elemToRemove results in a matrix PHI which removes the % 
%removed elements. %
% (Un)comment this section when element 1 and N must remain in the array 
remove_range = setdiff(1:N, [1 N]); %create an array without the first and last element number. 
randomorder = randperm(length(remove_range)); % create a random number order for elements to be removed 
elemToRemove = remove_range(randomorder(1:N-M)); %remove elements at randomorder number from the range (excl. 1&N) 
elemRemain=remove_range(randomorder(N-M+1:end)); %store remaining elements (for display purposes) sort(elemRemain);    
%sort remaining elements (exl. 1 and N) 
Phi(elemToRemove,:)=[]; %remove the random rows from Phi 

%% create the scanangle matrix 
angles=linspace(scan_angle_min,scan_angle_max,Ns); %define the angles representing the gridpoints 
alpha(1,:)=2*pi*(d/lambda)*sind(angles); %create a row vector with the 
% alpha per angle, total Ns alpha values (Ns scanangles) 
%now create Psi by reproducing it for every element with the delay, Psi now
%becomes a NxNs matrix with the steering vectors for every element in a row 
for n=1:N     
    a_row=exp(i*alpha*(n-1));     
    Psi(n,:)=a_row;  %Psi is the angle scanning matrix 
end

Theta=Phi*Psi; %create the observation matrix 
%create the signals received by the sparse array 
Y=Phi*X; %create the MxT measurement matrix, containing time samples of M elements 

%% USE MFOCUSS to find the solution. (the file mfocuss.m needs to be in the 
% same directory as this mfile, or in a "specified path" directory 
[Z_MF,gamma_ind,gamma_est,iters]=mfocuss(Theta,Y,lambda_m,'p',p,'MAX_ITERS' ,r); 

%% calculate the MSE with Frobenius norm (equal to the 2-norm for Nsamp=1) 
max_error_MF_MSE=norm((Psi*Z_MF)-X,'fro')^2;  
MSE=max_error_MF_MSE/norm(X,'fro')^2; 

%% Now that CS is performed is calculated the angle spectrum Py 
for theta = 1:Ns     
    Py(:,theta)= 1/Nsamp * sum(abs(Z_MF(theta,:).^2)); 
end

%% calculate the scanned beam pattern obtained with the true X 
X_meas=abs(X'*Psi); 
mx=max(X_meas); 
X_meas_norm=X_meas/mx; %% Plotting: plot the scanned beam pattern and MFOCUSS solution 
figure(1) %the MFOCUSS solution: 
plot(angles, pow2db(Py),'-*','MarkerSize',10,'LineWidth',1.5); 
grid on;
hold on; 
%the (normalized) scanned beam pattern of the fully dense ULA 
plot(angles,mag2db(X_meas_norm),'-r','LineWidth',1.5); 
hold on; 
%notice to adjust the legend when one of the above lines are (un)commented 
legend(['\bf\fontsize{18}MFOCUSS'],['\bf\fontsize{16}Scanned beam pattern (X*Psi)'],'location','SE')
ylim([-40 5]); %xlim([DOAangles1(1,nnn)-10 10]); 
title(['\bf\fontsize{22}',num2str(K),' targets at ', num2str(sort(DOAangles)),' degrees. Noise level = ',num2str(Noise),'dB.']); 
xlabel('\bf\fontsize{20}Broadside angle (degrees)'); 
ylabel('\bf\fontsize{20}Power (dB)'); 
axis tight; 
grid on; 
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. set(gca,'FontSize',18); %set fontsize of axis numbers 
%% display the location of the targets found with the very simple  
% "detection algorithm" findpeaks. Note that the a priori information of 
% the number of targets K is used! 
[~, locs_MFOCUSS]=findpeaks(abs(Py),'NPeaks',K,'SortStr','descend'); 
DOA=angles(1,locs_MFOCUSS) %display the MSE and number of iterations used by MFOCUSS 
MSE 
iters