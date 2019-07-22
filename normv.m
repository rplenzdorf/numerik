function N = normv(Nx,Ny)

A = -ones(Nx,Ny);

A(1,:) = 0;
A(end,:) = 4;
A(:,1) = 2;
A(:,end) = 6;

% A(4:7,4) = 2;
% A(4:7,7) = 6;
% A(4,4:7) = 0;
% A(7,4:7) = 4;
% 
% A(5:6,5:6) = 9;

N=A;

N = reshape(N,Nx*Ny,1)
% imagesc(A)
end