function [Dxc,Dyc,Lap,w0,w2,w4,w6] = matrizen(Nx,Ny,hx,hy)

Dyc = 1/(2*hy).* sten2mat([1 0 -1]',Nx,Ny);

Dxc = 1/(2*hx).* sten2mat([1 0 -1],Nx,Ny);

Lap = 1/hy^2.*(sten2mat([1 -2 1]',Nx,Ny)+sten2mat([1 -2 1],Nx,Ny));

w6 = sten2mat([-2 4 0],Nx,Ny)./hx^2 + sten2mat([-1 0 -1]',Nx,Ny)./hx^2;
w4 = sten2mat([-2 4 0]',Nx,Ny)./hy^2 + sten2mat([-1 0 -1],Nx,Ny)./hy^2;
w2 = sten2mat([0 4 -2],Nx,Ny)./hx^2 + sten2mat([-1 0 -1]',Nx,Ny)./hx^2;
w0 = sten2mat([0 4 -2]',Nx,Ny)./hy^2 + sten2mat([-1 0 -1],Nx,Ny)./hy^2;

end
