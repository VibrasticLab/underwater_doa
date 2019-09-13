function X = MRA_output(cita,deploy,f,c,d,snap,sigman,SNR)
%Signal generation
M = length(deploy);
N = length(cita);
A = zeros(M,N);
j = sqrt(-1);
for k = 1:N
    A(:,k) = exp(-j*2*pi*(f/c)*d*sin(cita(k)*pi/180)*deploy');
end
s = (1/sqrt(2))*(randn(N,snap)+j*randn(N,snap));
Xnf0 = A*s;
sigma2 = (sigman^2)*(10^(SNR/10));
X = sqrt(sigma2)*Xnf0+sigman*(1/sqrt(2))*(randn(M,snap)+j*randn(M,snap));