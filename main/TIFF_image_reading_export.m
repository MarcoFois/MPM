clear all; clc; close all;
A = imread('f5.tif');
A = A(:,:,1);
C = imbinarize(A);
Bbin = bwboundaries(C);

for k = 1:length(Bbin)
   boundary = Bbin{k};
end

figure()
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)

Lx = 500; Ly = 350;
nEx = 0.4*Lx; 
nEy = 0.4*Ly;
hl = [Lx/nEx Ly/nEy];
nl=2;
[xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
NVx = size(xv,2);                                                     
NVy = size(yv,1);                                                      
NV  = NVx*NVy; 
NP = (NVx-1)*(NVy-1);
Nex = NVx-1;
Ney = NVy-1;
x = xv(:); y = yv(:);

xL      = 0+(0.5*hl/nl):hl/nl:Lx;
yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
[xpp,ypp] = meshgrid(xL,yL);
xp = xpp(:);
yp = ypp(:);
[in,on] = inpolygon(xp(:),yp(:),boundary(:,2),boundary(:,1));

xP = xp(in);
yP = yp(in);
hP = 15.*sin(xP(:)./150);
hi = hP;
%  figure()
% %  scatter(xp(:),yp(:),10,'r','filled');
% scatter(x(:),y(:),13,'r+')
% hold on
% 
% scatter3(xP,yP,hP,7,'b','filled');
% axis equal

% xp = xP(:);
% yp = yP(:);
% hp = hP(:);
% 
% Lx = 15; Ly = 15;
% nEx = 60; nEy = 60;
% hl = [Lx/nEx Ly/nEy];
% 
% [xv,yv] = meshgrid(0.0:hl(1):Lx,0:hl(2):Ly);
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy;                                                
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% np4e_x = 4;
% np4e_y = 4;
% xL      = 0+(0.5*hl/np4e_x):hl/np4e_x:4;
% yL      = 5.5+(0.5*hl/np4e_y):hl/np4e_y: 9.5;%-(0.5*hl/np4e_y);
% [xp,yp] = meshgrid(xL,yL);
% 
% Xp = [xp(:),yp(:)];
% nmp = length(xp(:));
% 
% hp = 2*ones(size(xp(:)));
%  figure()
% %  scatter(xp(:),yp(:),10,'r','filled');
% scatter(x(:),y(:),13,'r+')
% hold on
% 
% scatter3(xp(:),yp(:),hp(:),7,'b','filled');
% axis equal

xp = xP(:);
yp = yP(:);
hp = hP(:);
%% Setting initial area, mass, density

Xp = [xp(:) yp(:)];
 
%   figure()
 shp = alphaShape(Xp(:,1),Xp(:,2));
Fplot =  plot(shp);
 set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
 axis([0 size(C,2) 0 size(C,1)])
 Asys = area(shp);

nmp = length(xp(:));
rhosy = 1200;
Msys = (Asys.*mean(hp))*rhosy;
Mpar = Msys/nmp;
Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
Vpu = []; %(dam*15*hp(1))/nmp*ones(nmp,1);

%% Connectivity - p2e, e2p, p2v, v2e

NP = (NVx-1)*(NVy-1);
gn = reshape(1:NV,NVy,NVx);
e2v = zeros(Nex*Ney,4);

iele = 1;
% - e2v - connect each element to its vertices
for i = 1:Nex
    for j = 1:Ney
        e2v(iele,1) = gn(j,i);
        e2v(iele,2) = gn(j+1,i);
        e2v(iele,3) = gn(j,i+1);
        e2v(iele,4) = gn(j+1,i+1);
        iele = iele +1;
    end
end

p2e  = reshape((permute((floor((yp(:)-min(y(:)))./(hl(2)))+1)+...
               (Ney).*floor((xp(:)-min(x(:)))./(hl(1))),[2 1])),nmp,1);%
%%
neon = length(unique(p2e));
p2v  = reshape((e2v(p2e(:),:)'),4,nmp)';

               
%% Constants
g = 9.81;
t = 0;                  
T = 1.2;                  
%% Grid quantities initialization
Mv = zeros(NV,1);
MvR = zeros(2*NV,1);
vv = zeros(2*NV,1);
Uhv = zeros(2*NV,1);

dv = 0.0;
%% Boundary 
 
bnd = [min(x) max(x) min(y)  max(y)];
[left] = find(x<=bnd(1));
bcxl = left;
[right] = find(x>=bnd(2));
bcxr = right; %[bcx; r];

[up] = find(abs(y-bnd(4))<=1e-7);
bcyu = up;


[down] = find(abs(y-bnd(3))<=1e-7);
bcyd = down; %[bcy;r];
clear r;

bcxl = 2*bcxl(:,1)-1;
bcxl = [bcxl zeros(size(bcxl))];

bcxr = 2*bcxr(:,1)-1;
bcxr = [bcxr zeros(size(bcxr))];

bcyu = 2*bcyu(:,1);
bcyu = [bcyu zeros(size(bcyu))];

bcyd = 2*bcyd(:,1);
bcyd = [bcyd zeros(size(bcyd))];

%% Material point quantities initialization
yd = 4;
slo = -1/4;
Z = 10+slo*x;
hp = reshape(hp,[nmp 1]);
Fp   = zeros(nmp,4) + repmat([1 0 0 1],nmp,1);
rho = rhosy*ones(nmp,1);
vp  = zeros(nmp,2); 
Mp =  Mpar*ones(nmp,1);%Vp.*rho;
Vp = Vpu;
Ap = Apu;
vnp = sqrt(vp(:,1).^2+vp(:,2).^2);
% Ap = Mp./(rho.*hp);
hp = reshape(hp,[nmp 1]);
Uh = hp.*vp;
phi = atan(slo);
F_11 = Uh(:,1).^2./hp(:,1)+0.5*g*cos(phi).*hp(:,1).^2; 
F_12 = Uh(:,1).*Uh(:,2)./hp(:,1);
F_21 = Uh(:,1).*Uh(:,2)./hp(:,1);
F_22 = Uh(:,2).^2./hp(:,1)+0.5*g*cos(phi).*hp(:,1).^2; 
momp = zeros(nmp,2);
FL = zeros(3,nmp);
L = zeros(2,2,nmp);
Ep = zeros(2,2,nmp);
up  = zeros(nmp,2);
xi = 5000;
Fb(:,1) = zeros(nmp,1);
Fb(:,2) = zeros(nmp,1);

it = 1;
Xplot(:,:,:,it) = [Xp(:,1) Xp(:,2) hp(:,1)];
MP(:,it) = Mp(:,1);
AP(:,it) = Ap(:,1);
RHO(:,it) = rho(:,1);
% VP(:,it) = Vp(:,1);
VNP(:,it) = vnp(:,1);
%%

avpr = zeros(2*NV,1);
Ap0 = Ap;
%% Voellmy rheology + Bingham stress model initialization

vis = 50; %50; % Pa s
ty = 2000; %2000;
a = 3/2;
alpha = (6*vis*norm(vp))./(hp*ty);
b = -(114/32+alpha);
C = 65/32;

s_xx = zeros(nmp,1);
s_xy = zeros(nmp,1);
s_yy = zeros(nmp,1); 
D_bar(1:3,1:3,1:nmp) = zeros;
invII = 0.5*sum(dot(D_bar,D_bar));
sig_xx = s_xx;
sig_xy = s_xy;
sig_yy = s_yy;

SXX(:,it) = sig_xx(:,1);
SVX = zeros(NV,1);
svx(:,it) = SVX;
%%
DATA = struct("x",xp,"y",yp,"Mp",Mp,"Ap",Ap,"vpx",vp(:,1),"vpy",vp(:,2),"Nex",nEx,"Ney",nEy,"hx",hl(1),"hy",hl(2),"hp",hp,"mom_px",momp(:,1),...
       "mom_py",momp(:,2),"g",g,"T",T,"phi",phi,"xi",xi,"vis",vis,"ty",ty,"rho",rhosy);
json = jsonencode(DATA)

FID = fopen("DATA.json","w");
fprintf(FID,json);
fclose(FID);









