%Fs Hz (samples per second) is the rate at the speech signal is sampled
%Fs=2000;
close all;
clear;
clc;

disp('load data')
[x_ori,Fs_ori]=audioread('/home/viblab/Documents/margiasih/data_rekaman/Track1_031.wav');

%cut sample data
%minValue = 1;
%maxValue = 100000;
%x1 = x_ori(minValue:maxValue,1);

x = resample(x_ori, 7, 10);
Fs=8000;

%sound(x1,Fs_ori);
%figure(1)
%stem(x);title('Recorded input signal Hydrophone');
%xlabel('Length of the input audio signal');
%ylabel('Amplitude of the input audio signal');

%Discrete cosine transform of the recorded signal
a0=dct(x);
panjang=length(a0);
%figure(2)
%stem(a0)
%axis([0 panjang -1 1]);
%title('Discrete cosine transform of the recorded signal');
%xlabel('Length of the DCT spectrum');
%ylabel('Amplitude of the DCT spectrum');

% Thresholding the spectrum to make it sparse
for i=1:1:panjang
if a0(i,1)<=0.04 && a0(i,1)>=-0.06
a0(i,1)=0;
else
a0(i,1)=a0(i,1);
end
end
%a0;
%figure(3)
%stem(a0)
%axis([0 panjang -1 1]);
%title('The Threshold spectrum');
%xlabel('The length of the threshold spectrum');
%ylabel('Amplitude of the threshold spectrum');
 
% Sparsity of the spectrum(K)and Length of the signal (N)
K=500;
N=panjang;
% Random measurement matrix
disp('Creating measurment matrix...');
A = randn(K,N);
A = orth(A')';
%figure(4)
%imagesc(A)
%colorbar;
%colormap('lines');
%title('Random Measurement matrix');
%disp('Done.');

% observations vector
y = A*a0;
%figure(5)
%plot(y)
%title('Observation Vector');

%initial guess = min energy
x0 = A'*y;

%solve the LP
tic
xp = l1eqpd(x0, A, [], y, 1e-2);
toc
%figure(6)
%plot(xp)
%axis([0 panjang -0.6 0.6]);
%title(' Reconstructed Spectrum using l1-minimization');

% Inverse dicrete cosine transform of reconstructed signal (IDCT)
Xrec=idct(xp);
%wavplay(Xrec,Fs)
%figure(7)
%stem(Xrec)
%title('Reconstructed signal at the receiver');
%xlabel('Length of the reconstructed signal using IDCT');
%ylabel('Amplitude of the reconstructed signal using IDCT');

% Calculating Absolute error between the reconstructed and actual signal
err=(max(abs(Xrec-x)));
%stem(err);
%figure(8)
%title(' Absolute Error of Reconstructed spectrum and Threshold spectrum ');
%xlabel('Length of the Maximum Absolute Error');
%ylabel('Maximum Absolute error')

audiowrite('/home/viblab/Documents/margiasih/Resample/resample7per10/Track1_031_re.wav',Xrec,Fs_ori);

