%Program ini digunakan untuk mengkonvolusikan sinyal suara dengan impulse response
%terdapat 3 input sinyal 
clc; clear; close all;
[y,Fs] = audioread('Sweep_Channel_1.wav');  		%  Sinyal impulse
[y2,Fs2] = audioread('Track_2_007.wav');  %  Sinyal asli yang akan dikonvolusi
ynew1=y(:,1)'; 	%membaca kolom pertama matriks jika stereo, kalau mono tidak usah
ynew2=y2(:,1)'; %membaca kolom pertama matriks jika stereo, kalau mono tidak usah 
hasil=conv(ynew2,ynew1); %proses konvolusi
out1normalized = hasil/max(hasil); %normalisasi sinyal terkonvolusi
audiowrite('/home/viblab/Documents/margiasih/olahdatasparse/konvolusi/konvolutor.wav',out1normalized,44100) %membuat sinyal terkonvolusi menjadi wav, fs=44100
