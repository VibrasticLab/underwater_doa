function [sudut]=rata2(a,b,c,d,e,f)
data=[a b c d e f];
meann = mean(data);
stdd = std(data);
I = bsxfun(@gt, abs(bsxfun(@minus, data, meann)), stdd);
out = find(I);
data(I)=[];
sudut=mean(data);
end