 %Program ini digunakan untuk mencari mse dari dua buah sinyal
%Semakin kecil nilai MSE maka 2 buah sinyal dianggap identik
clc; clear all; close all;
[in,Fs]=audioread('Angjing Terekam_01.wav');
[out,Fs1]=audioread('Angjing Konvolusi_01.wav');
if length(in)<length(out);
    out=out(1:length(in));
else
    in=in(1:length(out));
end
in=in./max(abs(in));
out=out/max(abs(out));
in=in(:);
out=out(:);
z=sum((in-out).^2);
error=z/length(in)  %Hasil MSE
