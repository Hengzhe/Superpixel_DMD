function A=CircAper(I,J,IC,JC,R)
%CircAper(I,J,IC,JC,R): Circle Aperture with Coordinate grid I,J, Circle center (IC,JC), radius
%R
r2=(I-IC).^2+(J-JC).^2;
A=double(r2<R^2);
end