% This programme is to estimate the time delay using CC,PHAT and ML
% Using two input sinal from two microphone are distributed spatially
% This programme refer to yushi zhang paper "A Comparative Study of Time
% Delay Techniques Using Microphone Arrays
function [ccestimation_val, phatestimation_val]=tde(noisysignal1, noisysignal2, fs1)

%signalpower1=0;
%signalpower2=0;
%for n=1:length(noisysignal1);
%    signalpower1=signalpower1+(abs(noisysignal1(n)))^2;
%end
%for i=1:length(noisysignal2);
%    signalpower2=signalpower2+(abs(noisysignal2(i)))^2;
%end

time=(0:length(noisysignal1)-1)/fs1;

cccorrelation=xcorr(noisysignal1,noisysignal2);
signallength=length(noisysignal1);

%% CC method:find the peak postion of the crosscorrelation;
%[ccmaximum, cctime]=max(cccorrelation);
%ccestimation=(signallength-cctime);
for l=1:length(cccorrelation)
    h(l)=cccorrelation(length(cccorrelation)+1-l);
end
[d,w]=max(h);
ccestimation=(signallength-w);



% GCC method; using the PHAT filter;
gcc=zeros((signallength*2-1),1);
phatfilter=zeros((signallength*2-1),1);
crossspectrum=fft(h);
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
gccestimation=(signallength-gcctime);


% show the SNR, and the three estimated time delay value;snr,ccestimation,gccestimation,mlestimation
% plot the three crosscorrelation;
lag=zeros((signallength*2-1),1);
for n=1:(signallength*2-1);
    lag(n)=signallength-n;
end

%if signalpower1>=signalpower2
    %Compute value (sec) time delay of 2 signal by using  CC,GCC, and ML methods
    ccestimation_val=ccestimation/fs1;
    phatestimation_val=gccestimation/fs1;
%else signalpower1<=signalpower2
%    ccestimation_val=-(ccestimation)/fs1;
%    phatestimation_val=-(gccestimation)/fs1;
%end
    
end
