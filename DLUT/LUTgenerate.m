%to generate the DigLUT, Digitize in polar coordinate, 
%The generated field (not digital field) is stored. 

% Parameters that need modification:
%       Npx, wavelength, wavelength0, MaxA, Ba, Bp

addpath gatbx-master

%parameters. 
Npx=10;

[j,i]=meshgrid(0:Npx-1,0:Npx-1);
Phase0= 2*pi/(Npx^2)*(Npx*i+j);
Phase0=reshape(Phase0,1,Npx*Npx);

wavelength=[0.650,0.52,0.445];%wavelength
Nw=length(wavelength);
wavelength0=2.05; %reference wavelength/um

Upx=zeros(Nw,Npx*Npx);
for w=1:length(wavelength)
    Upx(w,:)=exp(1i*(wavelength0/wavelength(w))*Phase0);
end

MaxA=8; %maximum amplitude. 
Ba=1; %bit for amplitude 
Bp=2; %bit for phase
lenLUT=2^(Nw*(Ba+Bp));

ErrorList=zeros(lenLUT,1,'single');
FieldList=zeros(lenLUT,Nw,'single');
spList=zeros(lenLUT,Npx^2,'logical');

% starting parallel pool
gcp;
tic
parfor num=1:lenLUT
    utar=MaxA*DecU(num,Nw,Ba,Bp);
    [sp,err]=GaSearch(utar,Upx,Npx,100,1000,0.2);
    FieldList(num,:)=objfunc_(sp,Upx);
    ErrorList(num)=err;
    spList(num,:)=sp;
end
toc

%low error proportion:
disp(sum(ErrorList<0.2)/lenLUT)


function [sp,err]=GaSearch(Utar,Upx,Npx,NIND,MAXGEN,ERRtar)
%
%GA

NVAR=Npx^2;

GGAP=0.9;

Chrom=crtbp(NIND,NVAR);

gen=0;

ObjV=objfunc(Chrom,Upx,Utar); %这里没有bs2rv，因为优化变量本身是二进制

trace=zeros(1,MAXGEN);
while gen<MAXGEN
    FitnV=ranking(ObjV);
    SelCh=select('sus',Chrom,FitnV,GGAP);
    Rec=recombin('xovsp',SelCh,0.7);
    Mut=mut(Rec);
    ObjVM=objfunc(Mut,Upx,Utar);
    [Chrom,ObjV]=reins(Chrom,Mut,1,1,ObjV,ObjVM);
    gen=gen+1;
    trace(gen)=min(ObjV);

    if min(ObjV)<ERRtar
        break
    end

end


[err,imin]=min(ObjV);
sp=Chrom(imin,:);

end

% plot(trace,'-*');
% 
% disp(trace(gen))


function f=objfunc(Chrom,Upx,Utar)

Nind=size(Chrom,1);
tmp=zeros(Nind,1);

for i=1:Nind
U=objfunc_(Chrom(i,:),Upx);
tmp(i)=sum(abs(U-Utar).^2); %Error calculation
end

f=tmp;
end

function f=objfunc_(sp,Upx)
Nw=size(Upx,1);
U=zeros(1,Nw);
for w=1:Nw
    U(w)=sum(Upx(w,:).*sp);
end
f=U;
end
