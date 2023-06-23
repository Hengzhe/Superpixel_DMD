%intensity modulation, chinese character 上海交通大学(shanghai jiaotong university)

% N=80;
% Utar=4*ones(N);
addpath("../util")

%%
%read image and process

Img=imread("Picture/ChineseCharacters/1.png");
N=54;
MaxA=8;
Img=PreProcess(Img,MaxA,N);
Utar=zeros(N,2*N,3);
%
Utar(1:N,1:N,2)=Img(:,(1:N)+N*(1-1));
%Utar(1:N,1:N,3)=Img(:,(1:N)+N*(1-1));
%
Utar(N+1:2*N,1:N,1)=Img(:,(1:N)+N*(2-1));
Utar(N+1:2*N,1:N,3)=Img(:,(1:N)+N*(2-1));
%
Utar(1:N,N+1:2*N,1)=Img(:,(1:N)+N*(3-1));
Utar(1:N,N+1:2*N,2)=Img(:,(1:N)+N*(3-1));
Utar(1:N,N+1:2*N,3)=Img(:,(1:N)+N*(3-1));
%
Utar(N+1:2*N,N+1:2*N,2)=Img(:,(1:N)+N*(4-1));
Utar(N+1:2*N,N+1:2*N,3)=Img(:,(1:N)+N*(4-1));
%
Utar(1:N,2*N+1:3*N,1)=Img(:,(1:N)+N*(5-1));
%Utar(1:N,2*N+1:3*N,2)=Img(:,(1:N)+N*(5-1));
%Utar(1:N,2*N+1:3*N,2)=Img(:,(1:N)+N*(5-1));
%
Utar(N+1:2*N,2*N+1:3*N,1)=Img(:,(1:N)+N*(6-1));
Utar(N+1:2*N,2*N+1:3*N,2)=Img(:,(1:N)+N*(6-1));
%Utar(:,N+1:2*N,1)=Img1;

%  Utar=8*ones(54*2,54*3,3);
 
%  Utar=Padding(Utar,20);

figure
imshow((1/MaxA*Utar).^2)
title("Target Image")

%%
%load DLUT
mfLUT=matfile("../DLUT/LUT-3.mat");
mfKDT=matfile("../DLUT/KDTmodel-3.mat");


%%
%generate the DMD pattern
tic
[DMD_SP,Error]=DMDgeneration_SP(Utar,mfLUT.wavelength,mfLUT.wavelength0,mfLUT.Npx,mfKDT.model,mfLUT.spList);
toc
% tic
% %DMD_ED=DMDgeneration_ED(Utar,mfLUT.wavelength,mfLUT.wavelength0,mfLUT.Npx);
% toc

figure
% subplot(1,2,1)
imshow(DMD_SP)
title("DMDPattern")
%  subplot(1,2,2)
%  imshow(DMD_ED)

figure
% subplot(1,2,1)
imagesc(abs(fftshift(fft2(DMD_SP))));
title("Fourier Plane")
%  subplot(1,2,2)
%  imagesc(abs(fftshift(fft2(DMD_ED))));


%%
%simulate the output image
um=1e-6;
mm=1e-3;
px_size=7.56*um;
wavelength=[0.650,0.520,0.445]*um;
wavelength0=mfLUT.wavelength0*um;
Npx=mfLUT.Npx;
focal_len=250*mm;
a=wavelength0*focal_len/(Npx^2*px_size);
R=a*0.6;
npad=30*Npx;
DMD_SP=Padding(DMD_SP,npad);
Uout_SP=DMDforvard(DMD_SP,px_size,wavelength,Npx,a,R,focal_len);
Uout_SP=InvPadding(Uout_SP,npad);
DMD_SP=InvPadding(DMD_SP,npad);
[M,N,l]=size(Uout_SP);
Uout_SP=Uout_SP(M:-1:1,N:-1:1,:);
%  Uout_ED=DMDforvard(DMD_ED,px_size,wavelength,Npx,a,R,focal_len);
eta=0.006;
figure
imshow(1/eta*abs(Uout_SP).^2)
title("Output Image")
%  figure
%  imshow(abs(1/eta*1/4*Uout_ED).^2)

%%
%preprocess of input image
function A=PreProcess(Img,MaxA,N)
    Img=imresize(Img,[N,6*N]);
    Img=im2gray(Img);
    Img=255*double(Img>100);
    A=MaxA*(1-1/255*double(Img));
end