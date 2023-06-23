%Holographic projection, chinese character 上海交通大学(shanghai jiaotong university)

% N=80;
% Utar=4*ones(N);
addpath("../util")
addpath("../MatDiffrac")

Img=imread("Picture/ChineseCharacters/1.png");
N=54;
MaxA=8;
Img=PreProcess(Img,MaxA,N);

mfLUT=matfile("../DLUT/LUT-3.mat");
mfKDT=matfile("../DLUT/KDTmodel-3.mat");

um=1e-6;
mm=1e-3;

px_size=7.56*um;

%first panel
Utar=zeros(2*N,3*N,3);
%
Utar(1:N,1:N,2)=Img(:,(1:N)+N*(1-1));
%
Utar(N+1:2*N,1:N,1)=Img(:,(1:N)+N*(2-1));
Utar(N+1:2*N,1:N,3)=Img(:,(1:N)+N*(2-1));

Utar1=FresnelIR(imresize(Utar,30),px_size/2.5/3,um*mfLUT.wavelength,40*mm,10);

%second panel
Utar=zeros(2*N,3*N,3);
%
Utar(1:N,N+1:2*N,1)=Img(:,(1:N)+N*(3-1));
Utar(1:N,N+1:2*N,2)=Img(:,(1:N)+N*(3-1));
Utar(1:N,N+1:2*N,3)=Img(:,(1:N)+N*(3-1));
%
Utar(N+1:2*N,N+1:2*N,2)=Img(:,(1:N)+N*(4-1));
Utar(N+1:2*N,N+1:2*N,3)=Img(:,(1:N)+N*(4-1));

Utar2=imresize(Utar,30);%FresnelIR(imresize(Utar,30),px_size/2.5/3,um*mfLUT.wavelength,50*mm,10);

%third panel
Utar=zeros(2*N,3*N,3);
%
Utar(1:N,2*N+1:3*N,1)=Img(:,(1:N)+N*(5-1));
%Utar(1:N,2*N+1:3*N,2)=Img(:,(1:N)+N*(5-1));
%Utar(1:N,2*N+1:3*N,2)=Img(:,(1:N)+N*(5-1));
%
Utar(N+1:2*N,2*N+1:3*N,1)=Img(:,(1:N)+N*(6-1));
Utar(N+1:2*N,2*N+1:3*N,2)=Img(:,(1:N)+N*(6-1));
%Utar(:,N+1:2*N,1)=Img1;
Utar3=FresnelIR(imresize(Utar,30),px_size/2.5/3,um*mfLUT.wavelength,-40*mm,10);

%  Utar=8*ones(54*2,54*3,3);
Utar=Utar1+Utar2+Utar3;

Utar=imresize(Utar,1/30);
 
%  Utar=Padding(Utar,20);

tic
[DMD_SP,Error]=DMDgeneration_SP(Utar,mfLUT.wavelength,mfLUT.wavelength0,mfLUT.Npx,mfKDT.model,mfLUT.spList);
toc
% tic
% %DMD_ED=DMDgeneration_ED(Utar,mfLUT.wavelength,mfLUT.wavelength0,mfLUT.Npx);
% toc

figure
% subplot(1,2,1)
imshow(DMD_SP)
%  subplot(1,2,2)
%  imshow(DMD_ED)

figure
% subplot(1,2,1)
imagesc(abs(fftshift(fft2(DMD_SP))));
%  subplot(1,2,2)
%  imagesc(abs(fftshift(fft2(DMD_ED))));

um=1e-6;
mm=1e-3;
px_size=7.56*um;
wavelength=[0.65,0.520,0.445]*um;
wavelength0=mfLUT.wavelength0*um;
Npx=mfLUT.Npx;
focal_len=100*mm;
a=wavelength0*focal_len/(Npx^2*px_size);
R=a*0.7;
npad=30*Npx;
DMD_SP=Padding(DMD_SP,npad);
Uout_SP=DMDforvard(DMD_SP,px_size,wavelength,Npx,a,R,focal_len);
Uout_SP=InvPadding(Uout_SP,npad);
DMD_SP=InvPadding(DMD_SP,npad);
[M,N,l]=size(Uout_SP);
Uout_SP=Uout_SP(M:-1:1,N:-1:1,:);
%  Uout_ED=DMDforvard(DMD_ED,px_size,wavelength,Npx,a,R,focal_len);
eta=0.023;
figure
imshow(abs(1/eta*1/4*Uout_SP).^2)
%  figure
%  imshow(abs(1/eta*1/4*Uout_ED).^2)

function A=PreProcess(Img,MaxA,N)
    Img=imresize(Img,[N,6*N]);
    Img=im2gray(Img);
    Img=255*double(Img>100);
    A=MaxA*(1-1/255*double(Img));
end