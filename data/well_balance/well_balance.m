% WELL-BALANCED


clear all;
clc;

Lx = 10; Ly = 10;
nEx = 50; %180
nEy = 50;
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

xL      = 0+(0.5*hl(1)/nl):hl(1)/nl:10;
yL      = 0+(0.5*hl(1)/nl):hl(1)/nl: Ly-(0.5*hl(1)/nl);
[xpl,ypl] = meshgrid(xL,yL);
L = Lx;

xpr = [];
ypr = xpr;

xp = [xpl(:); xpr(:)];
yp = [ypl(:); ypr(:)];

hp =  10.*ones(size(xp(:),1),1);

Xp = [xp(:),yp(:)];

Zplot =5.*exp(-2/5.*((xv-5).^2 + (yv-5).^2));  % choose this one for Z_1
%Zplot(xv>=4 & xv<=8 & yv>=4 & yv<=8) = 4;     % choose this one for Z_2
[gx,gy] = gradient(Zplot);

figure()
surf(xv,yv,Zplot);
title('normal')




Z = Zplot(:); %reshape(Zplot,NVy*NVx,1);
dZdx = gx(:); %reshape(gx,NV,1);
dZdy = gy(:); %reshape(gy,NV,1);
nmp = numel(xp(:));







 figure()
surf(xv,yv,Zplot)
hold on
scatter3(xp,yp,hp,7,'b','filled');
axis([0 10 0 10 0 10])
xlabel('x');
ylabel('y');



%% Constants
g     = 9.81;
xi    = 0;
vis   = 0;
ty    = 0;
T     = 1.5;

%% Material point quantities initialization
nmp   = numel(xp);
rhosy = 1000.0;
Msys  = sum (hp*hl(1)/nl.^2*rhosy);
Mp    = Msys/nmp * ones(nmp, 1);
Vp    = Mp./rhosy;
Ap    = Vp./hp;
vp    = zeros (nmp,2);

momp  = zeros (nmp,2);

Fb(:,1) = zeros (nmp,1);
Fb(:,2) = zeros (nmp,1);

%%
DATA = struct (
	   "x", xp, ...
	   "y", yp, ...
	   "Mp", Mp, ...
	   "Ap", Ap, ...
	   "vpx", vp(:,1), ...
	   "vpy", vp(:,2), ...
	   "Nex", nEx, ...
	   "Ney", nEy, ...
	   "hx", hl(1), ...
	   "hy", hl(2), ...
	   "hp", hp, ...
	   "mom_px", momp(:,1), ...
	   "mom_py", momp(:,2), ...
	   "g", g, ...
	   "T", T, ...
	   "xi", xi, ...
	   "vis", vis, ...
	   "ty", ty, ...
	   "rho", rhosy, ...
	   "Vp", Vp, ...
	   "Z", Z, ...
	   "dZdx", dZdx,
	   "dZdy", dZdy
	 );
json = jsonencode(DATA);

FID = fopen("DATA.json","w");
fprintf(FID,json);
fclose(FID);
