% This programme is to estimate the time delay using CC,PHAT and ML
% Using two input sinal from two microphone are distributed spatially
% This programme refer to yushi zhang paper "A Comparative Study of Time
% Delay Techniques Using Microphone Arrays

clc; clear all; close all;

% calculate the crosscorrelation of two noisy signal;
[noisysignal1,Fs]=audioread('Track 3_020.wav');
[noisysignal2,Fs1]=audioread('Track 4_020.wav');
%load doppler_440hz_1.mat;
%load doppler_440hz_2.mat;

%noisysignal1=wav(1:1000);
%noisysignal2=wav2(1:1000);

Fs=16000;

time=(0:length(noisysignal1)-1)/Fs;

cccorrelation=xcorr(noisysignal1,noisysignal2);
signallength=length(noisysignal1);

% CC method:find the peak postion of the crosscorrelation;
[ccmaximum, cctime]=max(cccorrelation);
ccestimation=abs(signallength-cctime);


% GCC method; using the PHAT filter;
gcc=zeros((signallength*2-1),1);
phatfilter=zeros((signallength*2-1),1);
crossspectrum=fft(cccorrelation);
for n=1:(signallength*2-1);
 phatfilter(n)=abs(crossspectrum(n));
 gcc(n)=crossspectrum(n)/phatfilter(n);
end
gcccorrelation=ifft(gcc);

% GCC method:find the peak postion of the filtered crosscorrelation;
for n=1:(signallength*2-1);
 gcccorrelation(n)=abs(gcccorrelation(n));
end
[gccmaximum,gcctime]=max(gcccorrelation);
gccestimation=abs(signallength-gcctime);


% Maximum Likelihood method:using the ML filter;
coherence=cohere(noisysignal1,noisysignal2);
mlcoherence=zeros((signallength*2-1),1);
coherencelength=length(coherence);
coherencenumber=(signallength*2-1)/coherencelength;
p=fix(coherencenumber);
if p<coherencenumber;
 for n=1:coherencelength;
 for m=1:p;
 mlcoherence(m+(n-1)*p)=coherence(n);
 end
 end
end
for n=((coherencelength*p+1):(signallength*2-1));
 mlcoherence(n)=coherence(coherencelength);
end
if p==coherencenumber;
 for n=1:coherencelength;
 for m=1:coherencenumber;
 mlcoherence(m+(n-1)*coherencenumber)=coherence(n);
 end
 end
end
for n=1:(signallength*2-1);
 squaremlcoherence(n)=(abs(mlcoherence(n))).^2;
ml(n)=squaremlcoherence(n)*crossspectrum(n)/(phatfilter(n)*(1-squaremlcoherence(n)));
end
mlcorrelation=ifft(ml);

% ML method:find the peak postion of ML correlation;
for n=1:(signallength*2-1);
 mlcorrelation(n)=abs(mlcorrelation(n));
end
[mlmaximum,mltime]=max(mlcorrelation);
mlestimation=abs(signallength-mltime);


% show the SNR, and the three estimated time delay value;snr,ccestimation,gccestimation,mlestimation
% plot the three crosscorrelation;
lag=zeros((signallength*2-1),1);
for n=1:(signallength*2-1);
    lag(n)=signallength-n;
end

% Plot CC, GCC, and ML 
subplot(3,1,1);
plot(lag,cccorrelation,'b')
legend('CC');
subplot(3,1,2);
plot(lag,gcccorrelation,'r')
ylabel ('cross-correlation');
legend('PHAT');
subplot(3,1,3);
plot(lag,mlcorrelation,'g')
xlabel ('time lag');
legend('ML'); 

% Plot signal input from 2 microphones 
figure(2)
subplot(2,1,1), plot(time,noisysignal1)
subplot(2,1,2), plot(time,noisysignal2)

%Compute value (sec) time delay of 2 signal by using  CC,GCC, and ML methods
ccestimation_val=ccestimation/Fs
gccestimation_val=gccestimation/Fs
mlestimation_val=mlestimation/Fs