function A=InvPadding(B,n)
%the inverse of Padding, InvPadding(Padding(A,n),n)=A
[M,N,~]=size(B);
A=B(n+1:M-n,n+1:N-n,:);
end