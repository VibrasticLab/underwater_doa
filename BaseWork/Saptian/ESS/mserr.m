%Program ini digunakan untuk mencari mse dari dua buah sinyal
%Semakin kecil nilai MSE maka 2 buah sinyal dianggap identik
clc; clear; close all;
[in,Fs]=audioread('Anjing_terekam_01.wav');
[out,Fs1]=audioread('Anjing_konvolusi_01.wav');
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
error=z/length(in);  %Hasil MSE
disp(error)
