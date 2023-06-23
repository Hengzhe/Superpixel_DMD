function Uout=FresnelIR(Uin,dx,w,z,NPad)
%This function simulate the scalar diffraction of light, based on Fresnel
%diffraction formula, and FFT of Impulse Response function for calculation:
%Uout=IFFT{FFT{Uin}FFT{h}}.
%
%Uout=FresnelIR(Uin,dx,w,z,Npad): 
%Uin the input field, Uout the output field. 
%dx the sample distance, w: wavelength, z: transport distance, NPad: number for padding.  
if length(w)==1

Uin=Padding(Uin,NPad);

[M,N]=size(Uin);
%evaluate the sampling criterion
L=max(dx*[M,N]);
if dx>w*z/L
   disp("WARNING: under-sampled regime for chirp function, result from FresnelIR can be invalid"); 
end

%Impulse Response function. 
j=dx*((1:N)-N/2);
i=dx*((1:M)-M/2);
[J,I]=meshgrid(j,i);
k=2*pi/w;
h=1/(1i*z*w)*exp(1i*k/(2*z)*(I.^2+J.^2));
h=h.*EllipseAper(I,J,0,0,M*dx/2,N*dx/2);
h=fftshift(h);
H=dx^2*fft2(h);


%output
Uout=ifft2(fft2(Uin).*H);
Uout=InvPadding(Uout,NPad);

else
   for iw=1:length(w)
      Uout(:,:,iw)=FresnelIR(Uin(:,:,iw),dx,w(iw),z,NPad); 
   end
end

end