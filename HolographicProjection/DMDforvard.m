function Uout=DMDforvard(DMD,px_size,wavelength,Npx,a,R,focal_len)
%==========================================================================
%Simulate the optical field after the whole process: DMD modulating, Forvard propagating, and filter, etc. 
%
%Uout=DMDforvard(DMD,wavelength,wavelength0,Npx,a,R,focal_len)
%input: DMD pattern, light wavelength, pin-hole position a,
%   superpixel size Npx, pin-hole radius R, focal length of lens. 
%output: output field, a matrix with size (M*N*Nw), (M*N) is size of DMD,
%   Nw the number of wavelength. 
%
%Written by Hengzhe Yan, 2022/12/25
%==========================================================================
Nw=length(wavelength);
[M,N]=size(DMD);
%single color
if Nw==1
    %DMD to focal plane. 
    Uf=fftshift(fft2(DMD));
    dfm=1/(M*px_size);
    dfn=1/(N*px_size);
    dxm=wavelength*focal_len*dfm;
    dxn=wavelength*focal_len*dfn;
    [nmat,mmat]=meshgrid(1:N,1:M);
    am=floor(Npx*a/dxm);
    an=floor(a/dxn);
    Rm=floor(R/dxm);
    Rn=floor(R/dxn);
    tmpf=zeros(M,N);
    tmpf(floor(M/2)-Rm:floor(M/2)+Rm,floor(N/2)-Rn:floor(N/2)+Rn)=Uf(floor(M/2)-Rm+am:floor(M/2)+Rm+am,floor(N/2+an)-Rn:floor(N/2)+Rn+an);
    tmpf=tmpf.*EllipseAper(mmat,nmat,M/2,N/2,Rm,Rn);
    
    %focal plane to image plane
    Uout=1/(M*N)*(fft2(fftshift(tmpf)));
    
%multi-color
else
    tmp=zeros(M,N,Nw);
    for w=1:Nw
       tmp(:,:,w)=DMDforvard(DMD,px_size,wavelength(w),Npx,a,R,focal_len); 
    end
    Uout=tmp;
end

end