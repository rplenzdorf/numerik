function temperature(sv,Nx,Re)
%%
clear
clc
close all
tic
%% Input
CFL = .8;
DFL = .1;

global Nx;
global sv;
global Re;
%% Initialisierung

Ny = Nx/sv;
hx = sv/(Nx-1);
hy = 1/(Ny-1);

y_ = (0:(Ny-1))'*hy;
x_ = (0:(Nx-1))'*hx;

x = kron(ones(Ny,1),x_);
y = kron(y_,ones(Nx,1));

nu = 1/Re;

[Dyc,Dxc,Lap,w0,w2,w4,w6] = matrizen(Nx,Ny,hx,hy);
Dxx = 1/hy^2.*((sten2mat([1 -2 1]',Nx,Ny)));

dia = @(in) spdiags(in,0,Nx*Ny,Nx*Ny);
    
tv = 1; %Videoframes

dt = .01; %Zeitschritt sp�ter angepasst
u = zeros(Nx*Ny,1);
v = u;
%% Normvektor, Randvektor

N = normv(Nx,Ny);

R = N ~= -1;

Lap_p = diag(R==0)*Lap +diag(R);

%% Dirichlet
D = 0.*y;

%% Neumann Korrekturmatrix
korrx = zeros(Nx*Ny);
korry = zeros(Nx*Ny);
for i = 1:Nx*Ny
    if N(i) == 0
        korry(i,i) = -1;
    elseif N(i) == 4
        korry(i,i) = 1;
    elseif N(i) == 2
        korrx(i,i) = 1;
    elseif N(i) == 6
        korrx(i,i) = -1;
    end
end

%% Norm Matrizen
Uo = sten2mat([0 -1 1]',Nx,Ny)./hy;
Vo = sten2mat([0 -1 1],Nx,Ny)./hx;
Uw = sten2mat([-1 1 0]',Nx,Ny)./hy;
Vw = sten2mat([-1 1 0],Nx,Ny)./hx;

for i = 1:Nx*Ny
    if N(i) == 0
        wr(i,:) = w0(i,:);
    elseif N(i) ==  2
        wr(i,:) = w2(i,:);
    elseif N(i) ==  4
        wr(i,:) = w4(i,:);
    elseif N(i) ==  6
        wr(i,:) = w6(i,:);
    else
        wr(i,:) = 0;
    end
end
%% Anfangswirbel
w = zeros(Nx*Ny,1);
%% Anfangstemperatur
T = zeros(Nx*Ny,1);
%% Zeitschleife
for t = 0:dt:10
    w = ~R.*w;
    
    Psi = linsolve(Lap_p,((R==0).*(-w)) + R.*D);

    for k = 1:(Nx*Ny)
        if N(k) == -1
            u(k) = Dyc(k,:)*Psi;
            v(k) = -Dxc(k,:)*Psi;
        elseif N(k) == 2
            u(k) = 0;
            v(k) = 0;
        elseif N(k) == 6
            u(k) = 0;
            v(k) = 0;
        elseif N(k) == 0
            u(k) = Dyc(k,:)*Psi;
            v(k) = 0;
            T(k) = 1;
        elseif N(k) == 4
            u(k) = Dyc(k,:)*Psi;
            v(k) = 0;
        end
    end
    
    wkorr = (wr * Psi - korrx*u*2./hx-korry*v*2./hy);
    w = w + wkorr;
    w(Nx*Ny-Nx+1) = 0;
    w(Nx*Ny) = 0;
    
    %% Temperatur
    a = u > 0;
    b = v > 0;
    
    Tabl = -Uw*(a.*u.*T) -Uo*(~a.*u.*T) -Vw*(b.*v.*T) -Vo*(~b.*v.*T) - 1.*Dxx*T;
    T = T + dt*Tabl;
    wkorr = (wr * Psi - korrx*u*2./hx-korry*v*2./hy);
    T = T + wkorr;
    T(Nx*Ny-Nx+1) = 0;
    T(Nx*Ny) = 0;
    Pt = -10*Dxc*T;
    
    %% Downwind
    a = u > 0;
    b = v > 0;
    
    wabl = -Uw*(a.*u.*w) -Uo*(~a.*u.*w) -Vw*(b.*v.*w) -Vo*(~b.*v.*w) +nu.*Lap*w +Pt; %RHS

    %% Zeitanpassung
    dt1 = CFL*hy/max(sqrt(u.^2+v.^2));
    dt2 = DFL*hy.^2/nu;
    dt = min(dt1,dt2);
    
    %% Euler
    w = w + dt*wabl;
    
    
    %% Reshape
    X = reshape(x,Nx,Ny);
    Y = reshape(y,Nx,Ny);
    PSI = reshape(Psi,Nx,Ny);
    U = reshape(u,Nx,Ny);
    V = reshape(v,Nx,Ny);
    W = reshape(w,Nx,Ny);
    NORM = reshape(N,Nx,Ny);
    DIRI = reshape(D,Nx,Ny);
    C = sqrt(U.^2+V.^2);
    Tr = reshape(T,Nx,Ny);
                    
    
    %% Plot  
    colormap jet
    shading interp
    
    imagesc(x_,y_,Tr')
    hold on
    
% %     interpolation
%     s = pcolor(x_,y_,Tr')
%     s.FaceColor = 'interp'
%     s.EdgeColor = 'none'
%     hold on

%     quiver(x,y,u,v,'w')
%     hold on
    title(['t:',num2str(t),' Re:',num2str(Re)])
    colorbar()
    set(gca,'YDir','normal')
    axis equal
    x0=100;
    y0=100;
    width=1200;
    height=width/sv;
    set(gcf,'position',[x0,y0,width,height])
%     set(gca,'colorscale','log')
    drawnow
    
    F(tv) = getframe(gcf)
    tv = tv + 1;
    
    hold off
end
%% Video
vidname = input("Unter welchem Namen sooll das Video gespeichert werden?(Vollst�ndiger Dateipfad m�gich)\n",'s')
video = VideoWriter(vidname,'MPEG-4');
open(video)
writeVideo(video,F)
close(video)
t1 = toc/60
end