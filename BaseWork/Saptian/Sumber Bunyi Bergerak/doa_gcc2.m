% direction estimation (azimuth phi) for 1 dim. microphone arrays
% using generalized cross correlation and speech activity detection
%
% (for use in oversampled filterbank)
%
% alg     PHAT-GCC
% x1,x2   microphone signals
% dx      microphone distance in meters
% N       signal frame length (512, if omitted)
% Fs      sampling frquency in Hz (16000, if omitted)
% phi     azimuth in degrees
%
%   Copyright 2006 Gerhard Doblinger, Vienna University of Technology
%   g.doblinger@tuwien.ac.at
%   http://www.nt.tuwien.ac.at/about-us/staff/gerhard-doblinger/
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.

% algorithm: C.H.Knapp, G.C.Carter, "the generalized correlation method for 
%            estimation of time delay",
%            IEEE Trans ASSP, vol. ASSP-24, Aug. 1976, pp. 320-327. 
clc; clear all; close all;

[x1,Fs] = audioread('Track 1_030.wav'); %Sinyal yang ditangkap hidrofon 1
[x2,Fs2] = audioread('Track 2_030.wav'); %Sinyal yang ditangkap hidrofon 2
%[X,Fs] = audioread('x_e_16.wav');
%x1=X(:,1);
%x2=X(:,8);

dx=0.3; %Jarak antar hidrofon

N = 2048;

Lov = 4;
M = round(N/Lov);                            % frame hop size

% create time window function

a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1)));

x1 = x1/max(x1);
x2 = x2/max(x2);
Nx = min(length(x1),length(x2));
Ndata = 1+ceil((Nx-N)/M);
Nf = N;
Nfh = Nf/2+1;
Sxy = zeros(Nfh,1);
Sx = zeros(Nfh,1);
Sy = zeros(Nfh,1);
delay = zeros(Ndata,1);

dx = abs(dx);
OV = 4;                                      % oversampling factor for interpolation filter
vs = 1492.1;                                 % acoustic waves propagation speed
Nd = 2+ceil(dx/vs*Fs);                       % max. delay (delay offset to obtain overall positive
                                             % delays in Cxy)
Nfo = OV*Nf;
Ndo = OV*Nd;

L = 2*Nd;

alpha = 0.8;                                 % forgetting factor of spectral power averaging
alpha1 = 1-alpha;
doa_threshold = 0.09;                        % speech activity threshold
                                             % CHANGE, if necessary
delay_old = Ndo;

Cmat = zeros(L*OV,Ndata);

m = 0;                                    % PHAT-GCC algorithm
  for n = 1:M:Nx-N+1;
    m = m+1;
    n1 = n:n+N-1; 
    X1 = fft(x1(n1).*w,Nf);
    X2 = fft(x2(n1).*w,Nf);
    X1 = X1(1:Nfh);
    X2 = X2(1:Nfh);
    Sxy = alpha*Sxy + alpha1*X1 .* conj(X2);
    Cxy = OV*real(ifft(Sxy./(abs(Sxy)+1e-4),Nfo));
    Cxy = [Cxy(Nfo-Ndo+1:Nfo) ; Cxy(1:Ndo)];
    Cmat(:,m) = Cxy; 
    [Cxymax,imax] = max(Cxy);
    if Cxymax > doa_threshold
       delay(m) = imax-1;
       delay_old = delay(m);
    else
       delay(m) = delay_old;
    end
  end
delay = (delay(1:m)-Ndo)/(OV*Fs);            % correct delay offset
phi = 180/pi*real(acos(vs*delay/dx)); 


close all
figurebackcolor = 'white';
pos = [0.01 0.5 0.49 0.42];
fp1 = figure('numbertitle','off','name','DOA estimation',...
	     'Units','normal','Position',pos);
colordef(fp1,figurebackcolor);
t = M/Fs*[0:length(delay)-1];
plot(t,phi)
xlabel('Time t in sec.'), ylabel('Azimuth \phi in deg.');
title(['d = ', num2str(100*dx),' cm,  ', 'F_s = ', num2str(Fs), ' Hz']);
grid on, 
axis tight
%axis([0 50 0 180]);
%set(gca,'PlotBoxAspectRatio',[1 0.5 0.5]);

% plot Cxy

pos = [0.5055 0.5 0.49 0.42];
fp1 = figure('numbertitle','off','name','Generalized cross correlation',...
	     'Units','normal','Position',pos);
colordef(fp1,figurebackcolor);

tau = 1000*linspace(-Nd/Fs,Nd/Fs,L*OV);
imagesc(t,tau,Cmat); colormap('jet'), colorbar
%map = colormap('gray');
%colormap(flipud(map));
%imagesc(t,tau,Cmat); % colorbar;
set(gca,'YDir','normal');
ylabel('Delay \tau in msec.');
xlabel('Time t in sec.');

figure(3)
plot(t,phi,'x')
h = lsline;
set(h(1),'color','r')
y=-4.2424*t+160;
hold on
plot(t,y,'g--')
xlabel('Waktu (s)'), ylabel('Sudut \phi (deg)');
title(['d = ', num2str(100*dx),' cm,  ', 'Kec. Kapal = ', num2str(8.8/max(t)), 'm/s']);
grid on, 
xlim([0 33])
legend('Sudut Estimasi','Regresi Sudut Est.','Lintasan')




