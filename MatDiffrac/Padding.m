function B=Padding(A,n)
%Padding(A,n), padding n zeros for two sides of matrix A
[M,N,L]=size(A);
tmp=zeros(M+2*n,N+2*n,L);
tmp(n+1:n+M,n+1:n+N,:)=A;
B=tmp;
end