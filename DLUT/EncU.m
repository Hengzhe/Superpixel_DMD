function X=EncU(U,Nw,Ba,Bp)
% Encode the amplitude and phase of U to integer.
%Ba and Bp are the bit number of amplitude and phase. 
if Nw~=length(U)
    disp('EncU Error: Nw and U length not match');
end
if max(abs(U))>1
    disp('EucU WARNING: abs(U) exceed 1.')
end
da=1/(2^Ba-1);
Ad=round(abs(U)./da);
dp=2*pi/2^Bp;
Pd=round(mod(angle(U),2*pi)./dp);

x=zeros(1,Nw*(Ba+Bp));
for w=1:Nw
    x((w-1)*(Ba+Bp)+1:(w-1)*(Ba+Bp)+Ba)=bitget(Ad(w),1:Ba);
    x((w-1)*(Ba+Bp)+Ba+1:(w-1)*(Ba+Bp)+Ba+Bp)=bitget(Pd(w),1:Bp);
end
X=2.^(0:Nw*(Ba+Bp)-1)*x';
X=X+1;
end