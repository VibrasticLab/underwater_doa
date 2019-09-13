function doa_est = Peaksearch(sspec, theta, minval, maxval, K)
%Assume that theta is arranged in the ascending order, and sspec is the corresponding value
%doa_est is arranged in the ascending order in the end
%Find the searching range
doa_est = zeros(1,K);
N = length(theta);
N1 = 1;
N2 = 1;
for k = 1:N
    if theta(k) < minval
        N1 = N1 + 1;
    end
    if theta(k) < maxval
        N2 = N2 + 1;
    end
end
theta_sel = theta(N1:N2);
s = sspec(N1:N2);
lengths = length(s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss = zeros(2,lengths-2);
t = 0;
for k = 2:lengths-1
    if (s(k)>s(k-1))&&(s(k)>s(k+1))
        t=t+1;
        ss(1,t) = s(k);
        ss(2,t) = k;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doaraw = zeros(1,K);
if t == 0
    doa_est = 90*ones(1,K);
elseif t >= K
    st = zeros(2,t);
    st(:,:) = ss(:,1:t);
    [s0,dix] = dsort(abs(st(1,:)));
    k = 1;
    while k <= K
        doaraw(k) = st(2,dix(k));
        k = k+1;
    end
    doaraw = sort(doaraw);
    doa_est(:) = theta_sel(doaraw(:));
else
    for k = 1:t
        doaraw(k) = ss(2,k);
    end
    doaraw(1:t) = sort(doaraw(1:t));
    for k = t+1:K
        doaraw(k) = doaraw(t);
    end
    doa_est(:) = theta_sel(doaraw(:));
end