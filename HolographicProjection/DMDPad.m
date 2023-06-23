function B=DMDPad(A)
%put the pattern A in the center of the DMD B, with size 1080x1920.
B=zeros(1080,1920);
[M,N]=size(A);
B(floor((1080-M)/2)+1:floor((1080-M)/2)+M,floor((1920-N)/2)+1:floor((1920-N)/2)+N)=A;
B=1-B;
end