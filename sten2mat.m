function D = sten2mat(s,Nx,Ny)
%%
close all
clc

%% Initialisierung

a = length(s);
k = size(s);
v = [((-a/2)+0.5):1:((a/2)-0.5)];
%%

% Zeilenvektor
if k(1,1) == 1 
    I = speye(Nx);
    A = spdiags(ones(Ny,1)*s,v,Ny,Ny);
    D = kron(A,I);
    
% Spaltenvektor
else
    s = s';
    I = speye(Ny);
    A = spdiags(ones(Nx,1)*s,v,Nx,Nx);
    D = kron(I,A);
end

%% Plot
% imagesc(full(D));
% colormap(jet(10));
% colorbar();