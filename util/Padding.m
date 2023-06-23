function B=Padding(A,n)
%Padding(A,n), padding n zeros for two sides of matrix A
[M,N,Nw]=size(A);
if Nw==1
tmp=zeros(M+2*n,N+2*n);
tmp(n+1:n+M,n+1:n+N)=A;
B=tmp;
else
   tmp2=zeros(M+2*n,N+2*n,Nw);
   for w=1:Nw
      tmp2(:,:,w)=Padding(A(:,:,w),n); 
   end
   B=tmp2;
end
end