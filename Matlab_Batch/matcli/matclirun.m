function matclirun(alamatin,alamatout,samplow,samphigh)

fprintf('input %s dan output %s \n',alamatin,alamatout)
fprintf('sample low %s dan high %s \n',samplow,samphigh)

disp('load data')
[x_ori,Fs_ori]=audioread(alamatin);

x = resample(x_ori, str2double(samplow), str2double(samphigh));
Fs=8000;

a0=dct(x);
panjang=length(a0);

for i=1:1:panjang
if a0(i,1)<=0.04 && a0(i,1)>=-0.06
a0(i,1)=0;
else
a0(i,1)=a0(i,1);
end
end

K=500;
N=panjang;
disp('Creating measurment matrix...');
A = randn(K,N);
A = orth(A')';
disp('Done.');

y = A*a0;
x0 = A'*y;

disp('Reconstruction');
tic
xp = l1eqpd(x0, A, [], y, 1e-2);
toc

Xrec=idct(xp);

err=(max(abs(Xrec-x)));

audiowrite(alamatout,Xrec,Fs_ori);

quit
