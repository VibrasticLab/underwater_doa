%Appendix C: m-file that simulates reconstruction of K targets. 
%m-file to research the possibilities of using compressive sensing (with 
%the MFOCUSS algorithm on a line array of 64 elements with 16 random  
%working elements. 
%clear the workspace and close the figures 
clc; 
%close all; 
clear all;  
nr=0;  %simulation number to show difference when noise is different.  
%       (this number sets the random generator to the specified start value)   

%PARAMETER DEFINITION: CHANGE THESE TO YOUR WHISHES %The array parameters 
N=64; %Total number of elements the dense ULA 
% M=16; %Number of random working elements that remain in the sparse array 
% K=3;  %Number of targets received by the array 
Nsamp = 1; %number of time snapshots/MMV used 
% max_scan_offbroadside=59; %define the maximum off boresight scan angle 

%Define environement parameters 
c = 1500;% physconst('LightSpeed'); 
fc = 1000;              % Operating frequency 3 GHz 
lambda = c/fc;         % Wavelength  
d=lambda/2;            % distance between elements   

%Define the discretization grid 
% scan_angle_min=-45; %the left edge of the discretization grid 
% scan_angle_max=135; %the right edge of the discretization grid 
% coverage=2*max_scan_offbroadside; 
% Ns=181; %no of discretization grid points in the search area see table below for examples:
% %assumed a discretization grid from -90 to 90 degrees setting Ns to:
% 181 results in a grid spacing of 1 degree 
% 361 results in a grid spacing of 0.5 degrees 
% 721 results in a grid spacing of 0.25 degrees 
% 1801 results in a grid spacing of 0.1 degrees 
%etc   

% %received signal parameters 
% Noise=-20; %received noise level in dB 
% Target_power=1;% power level of the targets, choose 1 to normalize to 0dB,  
% % if random strength targets, this defines the maximum power level, rember 
% to (un)comment line 69 for the use of random strength targets 
% rng(nr); 
%Define the DOA of the target, give up K values between +max_scan_offbroadside! 
%DOAangles = [-51 -6 43]; 
%or (un)comment the line below to choose K random ON GRID DOAs 
%DOAangles = randperm(coverage,K)-max_scan_offbroadside 
%use randperm to avoid double DOAs 
% or (un)comment the 2 lines below to choose K random OFF GRID DOAs 
% off_grid_distance=0.5 
%specify the off grid distance for the random DOAs 
% DOAangles = randperm(coverage,K)-max_scan_offbroadside+off_grid_distance 
%use randperm to avoid double DOAs   
% 
% %Define the MFOCUSS input parameters 
% lambda_m=db2pow(Noise); %regularization parameter, according to the  
% %description in the noisy case generally close to the noise variance. But 
% %also user defined values can be chosen
% p=0.8; %Parameter that trades off between the speed of convergence and the  
% %       sparsity in the solution. p=0.8 generally leads to good results. 
% r=1000; %maximum number of iterations MFOCUSS can use to find a solution. 
rng(nr); 
% %create a conventional ULA with 64 elements using the phased array toobox 
ULA64 = phased.ULA('NumElements',N,'ElementSpacing',d); 
ULA64.Element.BackBaffled = true; % this makes that the array only radiates to one side 

DOAangles = [-51 -6 43];
angs = DOAangles;  

Noise=-20; %received noise level in dB 
Target_power=1;% power level of the targets, choose 1 to normalize to 0dB, 

pos64 = getElementPosition(ULA64); 
nPower = db2pow(Noise);  %noise power at the element in Watts, Power = 1; % power of the incoming signals 
sPower = Target_power; % power of the incoming signals 
%sPower=Target_power*0.1*randi([1,10],1,K); %(un)comment to use random strength signals. 
% generate a multichannel signal received by the ula 
rng('default'); 
rs = rng(0); 
x64 = sensorsig(pos64/lambda,Nsamp,angs,nPower,sPower); 


% For K random DOA angels: 
plotResponse(ULA64,fc,c,'RespCut','Az','Format','Polar');