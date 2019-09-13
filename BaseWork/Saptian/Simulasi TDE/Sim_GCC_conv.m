%%Program ini mensimulasikan TDE dari sinyal uji tone 1000 Hz
%%Yang dikonvolusikan dengan impulse response ruangan
%%Dengan menambahkan noise kedalam sinyal uuji
%%Dengan tujuan dapat mengestimasi waktu tunda berdasarkan fungsi SNR
clc; clear all; close all
%% Sinyal yang diuji
%load mtlb
[a,fs]=audioread('Track4_020.wav'); %Sinyal impuls ruangan
Fs=16000;
t=0:1/Fs:0.5;
f=1000;
signal=sin(2*pi*f*t);
%signal=mtlb;
%signal=resample(signal,16000,Fs);
l1=length(signal);
p=1;
for n=1:l1;
    if 2.^n<(l1*2);
        p=p+1;
    end
end
l2=2.^p;
signal1=zeros(l2,1);
for n=1:l1;
    signal1(n)=signal(n);
end 
delay=50;
%signal1=conv(signal1,a);
%signal1=resample(signal1,16000,Fs);
signal2=lagmatrix(signal1,delay);
signal2(1:delay)=0;
signal1=signal1/max(signal1);
signal2=signal2/max(signal2);

%% noise yang ditambahkan pada sinyal
[a2,fs3]=audioread('whitenoise.wav'); %Sinyal pengganggu
a2=resample(a2,16000,fs3);
%a=resample(a,16000,fs);
%imps=a;%(1:length(signal1));
%sinyal yang ditambahkan ke sinyal matlab
%convsignal1=conv(signal1,imps);
%convsignal2=conv(signal2,imps);
l=length(signal1);
noise1=a2(1:l);
noise2=a2(100001:l+100000);
%noise1=randn(l,1);
%noise2=randn(l,1);
noise1=noise1/max(noise1);
noise2=noise2/max(noise2);
%% Menghitung SNR 
scale=0:0.01:10;
scalenumber=length(scale);
snr=zeros(1,scalenumber);
for t=1:scalenumber;
    originalsignal1=signal1*scale(t);
    originalsignal2=signal2*scale(t);
    noisysignal1=originalsignal1+noise1;
    noisysignal2=originalsignal2+noise2;
    signalpower=0;
    noisepower=0;
    for n=1:l2;
        signalpower=signalpower+(abs(originalsignal1(n)))^2;
        noisepower=noisepower+(abs(noise1(n)))^2+(abs(noise2(n)))^2;
    end
    signalpower=signalpower/l1;
    noisepower=noisepower/(l1*2);
    snr(t)=20*log10(signalpower/noisepower);

%% Menghitung waktu delay menggunakan cross correlation
[cc,lags]=xcorr(noisysignal1,noisysignal2);
[ccmaximum, cctime]=max(cc);
ccestimation(t)=abs(cctime-l);
% Menghitung waktu delay menggunakan generalized cross correlation phase
% transform
gcc=zeros((l*2)-1,1);
phatfilter=zeros((l*2)-1,1);
crossspectrum=fft(cc);
for n=(1:(l*2)-1);
    phatfilter(n)=abs(crossspectrum(n));
    gcc(n)=crossspectrum(n)/phatfilter(n);
end
gcccorrelation=ifft(gcc);
% GCC metode:mencari puncak dari filter crosscorrelation;
for n=1:(l*2-1);
gcccorrelation(n)=abs(gcccorrelation(n));
end
[gccmaximum,gccphattime]=max(gcccorrelation);
%[gccphatmaximum, gccphattime]=max(cc);
gccphatestimation(t)=abs(gccphattime-l);
end

%% plot signal 1 dan signal 2
subplot(1,2,1), myspectrogram(noisysignal1, 16000, [18 1], @hamming, 1024, [-45 -2], [1 -0.97], 'jet', true, 'per');
ylabel('Frekuensi')
xlabel('Sampel')
title('Signal 1')
subplot(1,2,2), myspectrogram(noisysignal2, 16000, [18 1], @hamming, 1024, [-45 -2], [1 -0.97], 'jet', true, 'per'); 
ylabel('Frekuensi')
xlabel('Sampel')
title('Signal 2')
figure(2)
subplot(1,2,1), myspectrogram(signal1, 16000, [18 1], @hamming, 1024, [-45 -2], [1 -0.97], 'jet', true, 'per');
ylabel('Frekuensi')
xlabel('Waktu')
title('Signal 1')
subplot(1,2,2), myspectrogram(signal2, 16000, [18 1], @hamming, 1024, [-45 -2], [1 -0.97], 'jet', true, 'per'); 
ylabel('Frekuensi')
xlabel('Waktu')
title('Signal 2')
figure(3)
subplot(1,2,1);
plot(lags,cc,'b')
legend('CC');
xlabel('Lags')
ylabel('Amplitude');
subplot(1,2,2);
plot(lags,gcccorrelation,'r')
legend('PHAT');
xlabel('Lags')
ylabel('Amplitude');
figure(4)
subplot(1,2,1);
stem(snr,ccestimation,'b')
ylabel('Waktu tunda')
xlabel('SNR')
legend('CC');
subplot(1,2,2);
stem(snr,gccphatestimation,'r')
ylabel('Waktu tunda')
xlabel('SNR')
legend('PHAT');

