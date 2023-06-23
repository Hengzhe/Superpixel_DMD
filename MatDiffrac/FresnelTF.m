function Uout=FresnelTF(Uin,dx,w,z,NPad)
%This function simulate the scalar diffraction of light, based on Fresnel
%diffraction formula, and sampled form of analytical transfer function H for calculation:
%Uout=IFFT{FFT{Uin}H}.
%
%Uout=FresnelTF(Uin,dx,w,z): 
%Uin the input field, Uout the output field. 
%dx the sample distance, w: wavelength, z: transport distance, NPad: number for padding.  

Uin=Padding(Uin,NPad);

[M,N]=size(Uin);
%evaluate the sampling criterion
L=min(dx*[M,N]);
if dx<w*z/L
   disp("WARNING: over-sampled regime for chirp function, result from FresnelTF can be invalid"); 
end

%transfer function
dfi=1/(M*dx);
dfj=1/(N*dx);
[J,I]=meshgrid(dfi*((1:N)-N/2),dfj*((1:M)-M/2));
H=exp(-1i*pi*w*z*(I.^2+J.^2));
H=fftshift(H);

%output
Uout=ifft2(fft2(Uin).*H);

end