function ni = peak_find(s)
lengths = length(s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss = zeros(2,lengths-2);
t = 0;
for k = 2:lengths-1
    if ((s(k)>s(k-1))&&(s(k)>s(k+1)))
        t=t+1;
        ss(1,t) = s(k);
        ss(2,t) = k;
    end
end
ni = zeros(1,t);
for k = 1:t
    ni(k) = ss(2,k);
end