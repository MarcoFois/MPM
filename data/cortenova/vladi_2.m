
clear all;
close all;
clc;

load('mask_in_vladi.mat')
load('dem_real_vladi.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowstart = 1;
rowend   = 175;

colstart = 1;
colend   = 165;

ndivrows = rowend - rowstart;
ndivcols = colend - colstart;

Xmin = 180;
Xmax = 520;
Ymin = 180;
Ymax = 520;

hx = hy = 5;

y = linspace (0, hy*ndivrows, ndivrows+1);
x = linspace (0, hx*ndivcols, ndivcols+1);

%y = linspace (0, hy*(Ymax-Ymin), Ymax-Ymin );
%x = linspace (0, hx*(Xmax-Xmin), Xmax - Xmin );

[X, Y] = meshgrid (x, y);
%Z = double (immagine_int);
Z = dem; %double(dem(Xmin:Xmax,Ymin:Ymax));

zplot = dem; % double( dem(Xmin:Xmax,Ymin:Ymax));

[x_n, y_n] = meshgrid (0:1:hx*ndivcols, 0:1:hy*ndivrows);
dem_n = interp2 (X, Y, dem, x_n, y_n, 'spline');
Z = dem_n;
%h = 0*Z; %ymin 1400 ymax=1700   xmin = 700 xmax = 900
%{
%Ymin = 1350; Ymax = 1450; Xmin = 700; Xmax = 900;
select = (Y>Ymin & Y<Ymax & X >Xmin & X < Xmax);

h(select)= 38.*double(mask_in(select)) ; %(700.*double(mask_in(select))) - Z(select);
h(h < 0) = 0;
%}


%h = 38.*double(mask_in);

%h(h < 0) = 0;
DX = 1; DY = 1;
%[xp, yp] = meshgrid (Xmin:DX:Xmax, Ymin:DY:Ymax);
%hp = interp2 (X, Y, h, xp, yp, 'spline');
select = (Y>0 & Y<175 & X >0 & X < 165);
hp = (38).*double(mask_in);


yy = linspace (0, hy*ndivrows, ndivrows+1);
xx = linspace (0, hx*ndivcols, ndivcols+1);
[Xx, Yy] = meshgrid (xx, yy);
xp = Xx;
yp = Yy;

%[xp, yp] = meshgrid (0:1.5:hx*165, 0:1.5:hy*175);
%hp = interp2 (Xx, Yy, h, xp, yp, 'linear');

hp = hp(:);
xp = xp(:);
yp = yp(:);

hzero = find (hp < 38);
hp(hzero) = [];
xp(hzero) = [];
yp(hzero) = [];

%data = jsonencode (struct ('X', X(:), 'Y', Y(:), 'Z', Z(:), 'h', h(:)));
%immagine_out = uint16(immagine_int);
%imwrite (immagine_out, "sfondo.tif")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Zz        = double (Z);
Zz(Zz<=0) = -1.0;
[gx, gy]  = gradient (Zz);

Z    = Zz(:);
dZdx = gx(:);
dZdy = gy(:);
%{
figure(1)
surf(X(Ymin:Ymax,Xmin:Xmax),Y(Ymin:Ymax,Xmin:Xmax),dem(Ymin:Ymax,Xmin:Xmax));
hold on
scatter3(xp,yp,hp,10,'r')
colormap('hsv');
%}
figure(1)
surf(X,Y,zplot);
hold on
scatter3(xp,yp,hp,10,'r')
colormap('hsv');

figure(2)
surf(x_n,y_n,dem_n);
hold on
scatter3(xp,yp,hp,10,'r')
colormap('hsv');

%figure(2)
%surf(X(100:600,240:700),Y(100:600,240:700),dem(100:600,240:700));
%colormap('hsv');

figure(3)
scatter3(xp,yp,hp,10,'r')

%{
figure(3);

mesh(X,Y,Zz)
hold on
mesh(X,Y,h)
xlabel('x');
ylabel('y');
view([70 30])
%}


at = atan(sqrt(gx.^2 + gy.^2));
%{
figure(3)
mesh(X,Y,at)
xlabel('x');
ylabel('y');
colormap('hsv');


figure(4)
surf(X(100:600,240:700),Y(100:600,240:700),at(100:600,240:700));
xlabel('x');
ylabel('y');
colormap('hsv');
%}
figure(5)
scatter3(xp(:),yp(:),hp(:))


%% Constants
g     = 9.81;
xi    = 200;
vis   = 50;
ty    = 2000;
T     = 10;

%% Material point quantities initialization
nmp   = numel(xp);
rhosy = 1291.0;
Msys  = sum (hp*DX*DY*rhosy);
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
	   "Nex", ndivcols, ...
	   "Ney", ndivrows, ...
	   "hx", hx, ...
	   "hy", hy, ...
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









