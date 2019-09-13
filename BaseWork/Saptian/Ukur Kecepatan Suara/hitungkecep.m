% Program ini digunakan untuk menghitung kec. suara menggunakan CC,PHAT and ML
% Menggunakan 2 microphone yg terdistribusi spasial sehingga terdapat waktu
% tunda. Input 2 sinyal yang terdelay akibat distribusi spasial mikrofon
% Panjang sinyal uji sebesar 0.03 sekon
clc; clear all; close all;

%% Hitung kros korelasi sinyal input
[sinyal1,Fs]=audioread('Track 2_003.wav');  %Sinyal yg ditangkap mik 1 //Track 2//
[sinyal2,Fs1]=audioread('Track 3_003.wav'); %Sinyal yg ditangkap mik 2 //Track 3//

time=(0:length(sinyal1)-1)/Fs;              %Buat waktu dari panjang sinyal

cccorrelation=xcorr(sinyal1,sinyal2);       %Korelasi 2 sinyal
signallength=length(sinyal1);

%% CC method:untuk mencari posisi puncak;
[ccmaximum, cctime]=max(cccorrelation);
ccestimation=abs(signallength-cctime); %Hasil estimasi waktu tunda menggunakan CC


%% GCC method; menggunakan PHAT filter;
gcc=zeros((signallength*2-1),1);
phatfilter=zeros((signallength*2-1),1);
crossspectrum=fft(cccorrelation);
for n=1:(signallength*2-1);
 phatfilter(n)=abs(crossspectrum(n));
 gcc(n)=crossspectrum(n)/phatfilter(n);
end
gcccorrelation=ifft(gcc);

%% GCC method:untuk mencari puncak sinyal yang ter-PHAT;
for n=1:(signallength*2-1);
 gcccorrelation(n)=abs(gcccorrelation(n));
end
[gccmaximum,gcctime]=max(gcccorrelation);
gccestimation=abs(signallength-gcctime);	%Hasil estimasi waktu tunda menggunakan PHAT


%% Maximum Likelihood method:menggunakan ML filter;
coherence=mscohere(sinyal1,sinyal2);
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
mlcorrelation=ifft(ml);	%Hasil estimasi waktu tunda menggunakan ML

%% ML method:untuk mencari puncak sinyal yang ter-ML;
for n=1:(signallength*2-1);
 mlcorrelation(n)=abs(mlcorrelation(n));
end
[mlmaximum,mltime]=max(mlcorrelation);
mlestimation=abs(signallength-mltime);

%%Mencari nilai lags dari sinyal yang terkorelasi
lag=zeros((signallength*2-1),1);
for n=1:(signallength*2-1);
    lag(n)=signallength-n;
end

%% Plot CC, GCC, and ML 
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

%% Plot signal input dari 2 microphones 
figure(2)
subplot(2,1,1), plot(time,sinyal1), title('Sinyal Hidrofon 1'),xlabel('waktu (s)'),ylabel('Amplitudo'),xlim([0 0.2])
subplot(2,1,2), plot(time,sinyal2), title('Sinyal Hidrofon 2'),xlabel('waktu (s)'),ylabel('Amplitudo'),xlim([0 0.2])


%% Hitung delay dari 2 sinyal menggunakan CC,PHAT, dan ML methods
ccestimation_val=ccestimation/Fs;
gccestimation_val=gccestimation/Fs;
mlestimation_val=mlestimation/Fs;

%% Hitung kecepatan suara bawah air jika s bernilai 9.42 m rumus v=s/t
kecepatan1=9.41/ccestimation_val
kecepatan2=9.41/gccestimation_val
kecepatan3=9.41/mlestimation_val
