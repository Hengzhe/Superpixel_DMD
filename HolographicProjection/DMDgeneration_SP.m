function [DMD,Error]=DMDgeneration_SP(Utar,wavelength,wavelength0,Npx,KDTmodel,spList)
%==========================================================================
% generate DMD pattern based on the multi-color superpixel method. 
%
%[DMD,Error]=DMDgeneration_SP(Utar,wavelength,wavelength0,Npx,KDTmodel,spList)
%
%e.g. DMDgeneration_SP(rand(100,100,2),[0.65,0.445],1.3,7,KDTmodel,spList)
%return the DMD pattern and evaluated Error, corresponding a all-bright picture
%for wavelength 0.65um and 0.445um, with reference wavelength 1.3um,
%superpixel size 7*7, KDTmodel and spList are from intended LUT. 
%
%written by Hengzhe Yan, 2022/12/25. 
%==========================================================================

[M,N,Nw]=size(Utar);
tmp=zeros(M*N,Nw);

%to multiply an offset phase for each superpixel.
[n,m]=meshgrid((0:N-1),(0:M-1));
P=zeros(M,N,Nw);
for w=1:Nw
P(:,:,w)=exp(-1i*2*pi*(wavelength0/wavelength(w))*(1/Npx).*(n+Npx*m));
end
Utar=Utar.*P;

%reshape and using knnsearch to find Idx
for w=1:Nw
   tmp(:,w)=reshape(Utar(:,:,w),M*N,1); 
end
Utar=[real(tmp),imag(tmp)];

[Idx,E]=knnsearch(KDTmodel,Utar);
Error=E;

%From Idx get superpixel and reshape back
Idx=reshape(Idx,M,N);
tmp=zeros(M*Npx,N*Npx);
for i=1:M
   for j=1:N
       sp=spList(Idx(i,j),:);
       tmp((i-1)*Npx+1:i*Npx,(j-1)*Npx+1:j*Npx)=reshape(sp,Npx,Npx);
   end
end
DMD=tmp;
end