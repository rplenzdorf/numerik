function [u,v] = geschw(N,Psi,hx,hy,Nx,Ny)

%%
A0 = sten2mat([-1 0 -1],Nx,Ny)+sten2mat([0 2 -2]',Nx,Ny);
A2 = sten2mat([-1 0 -1]',Nx,Ny)+sten2mat([0 2 -2],Nx,Ny);
A4 = sten2mat([-1 0 -1],Nx,Ny)+sten2mat([-2 2 0]',Nx,Ny);
A6 = sten2mat([-1 0 -1]',Nx,Ny)+sten2mat([-2 2 0],Nx,Ny);

Dyc = 1/(2*hy).* sten2mat([1 0 -1]',Nx,Ny);
Dxc = 1/(2*hx).* sten2mat([1 0 -1],Nx,Ny);

%%
for k = 1:(Nx*Ny)
    
    if N(k) == 0
        u(k) = 0%A0(k,:)*Psi;
        v(k) = -Dxc(k,:)*Psi;
    elseif N(k) == 4
        u(k) = 0%A4(k,:)*Psi;
        v(k) = -Dxc(k,:)*Psi;
    elseif N(k) == 2
        u(k) = Dyc(k,:)*Psi;
        v(k) = 0%-A2(k,:)*Psi;
    elseif N(k) == 6
        u(k) = Dyc(k,:)*Psi;
        v(k) = 0%-A6(k,:)*Psi;
    else 
        u(k) = Dyc(k,:)*Psi;
        v(k) = -Dxc(k,:)*Psi;
    end
end

u = u';
v = v';
end