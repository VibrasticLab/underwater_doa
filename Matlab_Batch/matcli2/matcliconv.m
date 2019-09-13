function matcliconv(wavin,wavout)

[y,Fs] = audioread('konvolutor.wav'); 
[y2,Fs2] = audioread(wavin);
ynew1=y(:,1)';
ynew2=y2(:,1)'; 
hasil=conv(ynew2,ynew1);
out1normalized = hasil/max(hasil); 
audiowrite(wavout,out1normalized,44100) 

quit
