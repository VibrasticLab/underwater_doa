%Ekstraksi Impulse Response Menggunakan Metode Exponential Sine Swep
%Metode ini dikembangkan oleh Angelo Farina dari Italia
%Judul paper "Simultaneous Measurement of Impulse Response and Distortion %  with a Swept Sine Technique"
% http://pcfarina.eng.unipr.it/
%Pada script ini terdapat 3 input dengan data dan Fs yang sama panjang, yakni: Sinyal ESS, sinyal ESS terekam, invers ESS
clear all;
close all;
s=audioread('Sweep Channel 1.wav'); 		%Sinyal ESS yang dibangkitkan diruanagan
s=s';
[s_hall,fs]=audioread('Track 2_007.wav'); 	%Sinyal ESS terekam pada ruangan
s_hall=s_hall';
sinv=audioread('Inverse Filter.wav'); 		%Sinyal invers ESS
sinv=sinv';
s_hall=s_hall(1:length(sinv)); 				%menyamakan panjang siny. rekaman dg inverse filter
s=s(1:length(sinv));

s_hall=s_hall/max(s_hall);					%Normalisasi sinyal sehingga amplitud maks = 1
SINV=fft(sinv);
S_HALL = fft(s_hall); 
S=fft(s);

%---------------proses konvolusi sinyal terekam and inverse filter----------------- %
x1=sinv;
x2=s_hall;
x1l=[x1 zeros(1,length(x1)-1)]; %zero pad input vektor untuk menghindari circular convolution
x2l=[x2 zeros(1,length(x2)-1)];

X1L=fft(x1l);
X2L=fft(x2l);
YL=X1L.*X2L; 		%metode konvolusi menggunakan perkalian sinyal pada domain frekuensi 
imp=real(ifft(YL)); %sinyal impulse terekstraksi 
imp=imp/max(imp);	%Normalisasi
time=(0:length(imp)-1)/fs;

figure(1)
NFFT=length(s_hall); 

f = fs/2*linspace(0,1,NFFT/2); 								%create plotting vector
subplot(3,1,1), loglog(f,abs(S(1:NFFT/2)/max(abs(S)))) 		%Plot  amplitude spektrum sinyal sweep
title('Spektrum Sinyal Exponential Sine Sweep')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')
subplot(3,1,2), loglog(f,abs(S_HALL(1:NFFT/2))/max(abs(S_HALL))) % Plot amplitude spektrum sinyal terekam
title('Spektrum Sinyal Rekaman Exponential Sine Sweep')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')

subplot(3,1,3), loglog(f,abs(SINV(1:NFFT/2))/max(abs(SINV))) 	% Plot amplitude spektrum invers filter
title('Spektrum Sinyal Inverse Filter Exponential Sine Sweep')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')


figure(2)
subplot(3,1,1), plot(s_hall) %Plot domain waktu sinyal ESS
title('Sinyal Rekaman Exponential Sine Sweep')
ylabel('Amplitude')
xlabel('Time (ms)')
subplot(3,1,2), plot(s) %Plot domain waktu sinyal terekam
title('Sinyal Exponential Sine Sweep')
ylabel('Amplitude')
xlabel('Time (ms)')
subplot(3,1,3), plot(sinv)% Plot domain waktu sinyal invers
title('Exponential Sine Sweep Inverse Filter')
ylabel('Amplitude')
xlabel('Time (ms)')

figure(3)
plot(imp) %plot impulse response hasil konvolusi antara siny. terekam dg invers filternya
title('Impulse Response Hasil Ekstraksi')
ylabel('Amplitude (dB)')
xlabel('Time (ms)')
IR_abs = abs(imp);
IR_abs = IR_abs/max(IR_abs); %normalisasi

figure(4)
time=(0:length(IR_abs)-1)/fs;
plot(time,20*log10(IR_abs)); %Plot impulse response dalam skala logaritmik
title('Impulse Response Hasil Ekstraksi')
ylabel('Amplitude (dB)')
xlabel('Time (ms)')

