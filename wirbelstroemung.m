%%
clear
clc
close all
tic
% wie gut funktioniert github eigentlich?
% ziemlich gut
%% Input
sv = 4;

Nx = 160;

CFL = .8;
DFL = .1;

Re = 1000;

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

dia = @(in) spdiags(in,0,Nx*Ny,Nx*Ny);
    
tv = 1; 

dt = .01;
u = zeros(Nx*Ny,1);
v = u;
%% Normvektor, Randvektor

N = normv(Nx,Ny);

R = N ~= -1;

Lap_p = diag(R==0)*Lap +diag(R);

%% Dirichlet
D = -R.*y;

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
w = 1./((x-0.5).^2+(y-0.5).^2+0.001);

w = -w/10;
w = ~R.*w;

%% Zeitschleife
for t = 0:dt:4
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
        elseif N(k) == 4
            u(k) = Dyc(k,:)*Psi;
            v(k) = 0;
        end
        
%         if R(k) == 1
%             w(k) = Dxc(k,:)*v-Dyc(k,:)*u;
%         end
    end
    
    wkorr = (wr * Psi - korrx*u*2./hx-korry*v*2./hy);
    w = w + wkorr;
    w(Nx*Ny-Nx+1) = 0;
    w(Nx*Ny) = 0;
    %% Penalty
    obj1 = (x-.5).^2+((y -.5)).^2<.2^2;
    
    rx = .5;
    ry = .5;
    phi = 5*t;
    o1 = (x-rx)*cos(phi)+(y-ry)*sin(phi)+rx;
    o2 = (y-ry)*cos(phi)-(x-rx)*sin(phi)+ry;
    obj2 = (o1-rx).^10+((o2-ry)/4).^10<.03.^10;
    
    obj = obj1;
    
    c = 1/dt;
    
    fx = obj.*(0-u)*c;
    fy = obj.*(0-v)*c;
    
    P = Dxc*fy-Dyc*fx;
    
    %% Downwind
    a = u > 0;
    b = v > 0;
    
    wabl = -Uw*(a.*u.*w) -Uo*(~a.*u.*w) -Vw*(b.*v.*w) -Vo*(~b.*v.*w) +nu.*Lap*w+P;
    
    dt1 = CFL*hy/max(sqrt(u.^2+v.^2));
    dt2 = DFL*hy.^2/nu;
    dt = min(dt1,dt2);
    
    %% Euler
    w = w + dt*wabl;
%     w = ~R.*w;
    
    %% Kraft
    Fx = norm(fx);
    Fy = norm(fy);
    Fn = sqrt(Fx.^2+Fy.^2);
    
    %% Reshape
    X = reshape(x,Nx,Ny);
    Y = reshape(y,Nx,Ny);
    PSI = reshape(Psi,Nx,Ny);
    U = reshape(u,Nx,Ny);
    V = reshape(v,Nx,Ny);
    W = reshape(w,Nx,Ny);
    NORM = reshape(N,Nx,Ny);
    DIRI = reshape(D,Nx,Ny);
    OBJ = reshape(obj,Nx,Ny);
    C = sqrt(U.^2+V.^2);
                    
    
    %% Plot  
    colormap jet
    shading interp
    
    imagesc(x_,y_,C')
    hold on
    rectangle('Position',[.5-0.2 .5-0.2 .4 .4],'FaceColor',[.5 .5 .5],'Curvature',[1 1])
    hold on
%     quiver(x,y,u,v,'w')
%     hold on
%     quiver(.5,.5,Fx/Fn/2,Fy/Fn/2)
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

video = VideoWriter('E:\seife\Dokumente\Studium\Videos Numerik\kar_200-4_c.avi','Uncompressed avi');
open(video)
writeVideo(video,F)
close(video)
t1 = toc/60