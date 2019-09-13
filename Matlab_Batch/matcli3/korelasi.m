%Program ini digunakan untuk mencari hubungan/ korelasi dari 2 buah sinyal
%Pada TA ini digunakan untuk validasi impulse respon
%Dengan cara mengkorelasikan siyal rekaman suara anjing bawah air dan sinyal anjing hasil konvolusi
%Impulse response pada kolam, sinyal dikatakan identik jika nilai c mendekati 1 dan plot xcorr simetris
clc; clear all; close all;
[x,fs1]=audioread('Anjing terekam_01.wav'); 	%Sinyal yang terekam
[y,fs2]=audioread('Anjing konvolusi_01.wav');	%Sinyal hasil konvolusi
[acor,lag] = xcorr(x,y);						%Korelasi silang sinyal terekam dan hasil konv
subplot(1,2,1),plot(lag,acor)					%Plot korelasi
grid on
xlabel('Lags')
ylabel('Amplitude')
title('Korelasi Sinyal')
[acor,lag] = xcorr(y);							%Autokorelasi sinyal terekam
subplot(1,2,2),plot(lag,acor)					%Plot autokorelasi
grid on
xlabel('Lags')
ylabel('Amplitude')
title('Korelasi Sinyal')
title('Autokorelasi Sinyal Terekam')
c=corrcoef(x,y)									%Nilai koefisien korelasi 2 sinyal