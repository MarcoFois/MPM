
clear all;
close all;
clc;

ndivrows = 200;
ndivcols = 40;

hx = hy = .5;

y = linspace (0, hy*ndivrows, ndivrows+1);
x = linspace (0, hx*ndivcols, ndivcols+1);

[X, Y] = meshgrid (x, y);

Z = 10 - sqrt (10.^2 - (X - 10).^2) +  .3 * Y;


h = 0*Z;
Ymin = 70; Ymax = 100; Xmin = 0; Xmax = 20;
select = (Y>Ymin & Y<Ymax & X >Xmin & X < Xmax);

h(select)= 30 - Z(select);
h(h < 0) = 0;

surf (X, Y, Z)
hold all
surf (X, Y, Z+h)
axis equal


DX = .1; DY = .1;
[xp, yp] = meshgrid (Xmin:DX:Xmax, Ymin:DY:Ymax);

hp = interp2 (X, Y, h, xp, yp, 'nearest');
hp = hp(:);
xp = xp(:);
yp = yp(:);

hzero = find (hp <= 1);
hp(hzero) = [];
xp(hzero) = [];
yp(hzero) = [];

immagine_out = uint16(Z);
imwrite (immagine_out, "sfondo.tif")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Zz        = double (Z);
[gx, gy]  = gradient (Zz);

Z    = Zz(:);
dZdx = gx(:);
dZdy = gy(:);


figure()
mesh(X,Y,Zz)
axis equal

figure()
scatter3(xp(:),yp(:),hp(:))


%% Constants
g     = 9.81;
xi    = 200;
vis   = 50; 
ty    = 2000;
T     = 10;

%% Material point quantities initialization
nmp   = numel(xp);
rhosy = 1100.0;
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









