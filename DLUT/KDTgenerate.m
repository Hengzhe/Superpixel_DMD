mf=matfile('LUT-3.mat');

FieldList=mf.FieldList;
spList=mf.spList;

X=[real(FieldList),imag(FieldList)];

tic
model=KDTreeSearcher(X,'distance','Euclidean');
toc

