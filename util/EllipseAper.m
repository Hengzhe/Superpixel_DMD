function A=EllipseAper(I,J,IC,JC,RI,RJ)
%CircAper(I,J,IC,JC,RI,RJ): Circle Aperture with Coordinate grid I,J,
%Circle center (IC,JC), radius (RI,RJ)
r2=1/RI^2*(I-IC).^2+1/RJ^2*(J-JC).^2;
A=double(r2<1);
end