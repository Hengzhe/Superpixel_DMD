function U=DecU(X,Nw,Ba,Bp)
% to decode the integer X to complex vector with length Nw
X=X-1;
x=bitget(X,1:Nw*(Ba+Bp));
Ad=zeros(1,Nw);
Pd=zeros(1,Nw);
da=1/(2^Ba-1);
dp=2*pi/2^Bp;
for w=1:Nw
    Ad(w)=sum(x((w-1)*(Ba+Bp)+1:(w-1)*(Ba+Bp)+Ba).*(2.^(0:Ba-1)));
    Pd(w)=sum(x((w-1)*(Ba+Bp)+Ba+1:(w-1)*(Ba+Bp)+Ba+Bp).*(2.^(0:Bp-1)));
end
U=(da*Ad).*exp(1i*dp*Pd);

end