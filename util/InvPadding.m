function A=InvPadding(B,n)
%the inverse of Padding, InvPadding(Padding(A,n),n)=A
[M,N,Nw]=size(B);
for w=1:Nw
A(:,:,w)=B(n+1:M-n,n+1:N-n,w);
end
end