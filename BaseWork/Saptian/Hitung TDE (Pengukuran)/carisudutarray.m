clear all;
close all;
clc;
direktori=pwd;
addpath(direktori);
listfile=dir('*.wav');
hasil_estimasi=zeros(100,12);
hasil_sudut=zeros(100,12);
hasil_rata=zeros(100,2);
fprintf('Proses data TA Saptian\nEstimasi Waktu Tunda...');
for i=1:100;
    fprintf('%02d%%',floor((i)/(100)*100));
    name1=sprintf('Track 1_017.wav',i);%Sinyal yang ditangkap hidrophone 1
    name2=sprintf('Track 2_017.wav',i);%Sinyal yang ditangkap hidrophone 2
    name3=sprintf('Track 3_017.wav',i);%Sinyal yang ditangkap hidrophone 3
    name4=sprintf('Track 4_017.wav',i);%Sinyal yang ditangkap hidrophone 4
    [a, fs1]=audioread(name1);
    [b, fs2]=audioread(name2);
    [c, fs3]=audioread(name3);
    [d, fs4]=audioread(name4);
    %% estimasi waktu delay
    [hasil_estimasi(i,1),hasil_estimasi(i,2)]=tde(a, b, fs1);%Estimasi wktu tunda H1 & H2
    [hasil_estimasi(i,3),hasil_estimasi(i,4)]=tde(a, c, fs1);%Estimasi wktu tunda H1 & H3
    [hasil_estimasi(i,5),hasil_estimasi(i,6)]=tde(a, d, fs1);%Estimasi wktu tunda H1 & H4
    [hasil_estimasi(i,7),hasil_estimasi(i,8)]=tde(b, c, fs1);%Estimasi wktu tunda H2 & H3
    [hasil_estimasi(i,9),hasil_estimasi(i,10)]=tde(b, d, fs1);%Estimasi wktu tunda H2 & H4
    [hasil_estimasi(i,11),hasil_estimasi(i,12)]=tde(c, d, fs1);%Estimasi wktu tunda H3 & H4
    %% estimasi sudut datang terhadap array  mic 1: hasil_estimasi(i,ganjil)=estimasi CC & hasil_estimasi(i,ganjil)=estimasi PHAT
    %a,b ; a,c ;a,d
    [hasil_sudut(i,1)]=angles(hasil_estimasi(i,1),1492.1,0.3,0.45);% 1492.1=kecepatan; 0.3=jarak a-b; 0.45=jarak titik tengah array
    [hasil_sudut(i,2)]=angles(hasil_estimasi(i,2),1492.1,0.3,0.45);% 1492.1=kecepatan; 0.3=jarak a-b; 0.45=jarak titik tengah array
    [hasil_sudut(i,3)]=angles(hasil_estimasi(i,3),1492.1,0.6,0.45);% 1492.1=kecepatan; 0.3=jarak a-c; 0.45=jarak titik tengah array
    [hasil_sudut(i,4)]=angles(hasil_estimasi(i,4),1492.1,0.6,0.45);% 1492.1=kecepatan; 0.3=jarak a-c; 0.45=jarak titik tengah array
    [hasil_sudut(i,5)]=angles(hasil_estimasi(i,5),1492.1,0.9,0.45);% 1492.1=kecepatan; 0.3=jarak a-d; 0.45=jarak titik tengah array
    [hasil_sudut(i,6)]=angles(hasil_estimasi(i,6),1492.1,0.9,0.45);% 1492.1=kecepatan; 0.3=jarak a-d; 0.45=jarak titik tengah array
    %b,c ; b,d
    [hasil_sudut(i,7)]=angles(hasil_estimasi(i,7),1492.1,0.3,0.15);
    [hasil_sudut(i,8)]=angles(hasil_estimasi(i,8),1492.1,0.3,0.15);
    [hasil_sudut(i,9)]=angles(hasil_estimasi(i,9),1492.1,0.6,0.15);
    [hasil_sudut(i,10)]=angles(hasil_estimasi(i,10),1492.1,0.6,0.15);
    %c,d
    [hasil_sudut(i,11)]=angles(hasil_estimasi(i,11),1492.1,0.3,0.45);
    [hasil_sudut(i,12)]=angles(hasil_estimasi(i,12),1492.1,0.3,0.45);
    %%Penentuan sudut utama
    [hasil_rata(i,1)]=rata(hasil_sudut(i,1),hasil_sudut(i,3),hasil_sudut(i,5),hasil_sudut(i,7),hasil_sudut(i,9),hasil_sudut(i,11));
    [hasil_rata(i,2)]=rata(hasil_sudut(i,2),hasil_sudut(i,4),hasil_sudut(i,6),hasil_sudut(i,8),hasil_sudut(i,10),hasil_sudut(i,12));
    fprintf('\b\b\b');
end
fprintf('\bSelesai!\n');
