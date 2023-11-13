
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load elevation.octbin.gz
ndivcols = 135;
ndivrows = 165;
hx = (X(1, end) - X(1, 1))/ndivcols;
hy = (Y(end, 1) - Y(1, 1))/ndivrows;;

xptc = dlmread ("particles.csv", ',', 1, 0);
xp = xptc(:, 1);
yp = xptc(:, 2);
hp = xptc(:, 4);

hzero = find (hp <= 1);
hp(hzero) = [];
xp(hzero) = [];
yp(hzero) = [];

%data = jsonencode (struct ('X', X(:), 'Y', Y(:), 'Z', Z(:), 'h', h(:)));
immagine_out = uint16(Z);
imwrite (immagine_out, "sfondo.tif")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
mesh(X,Y,Z)
hold all
scatter3(xp(:),yp(:),hp(:)+interp2(X,Y,Z,xp(:), yp(:)))


[gx, gy]  = gradient (Z, hx, hy);

Z    = Z(:);
dZdx = gx(:);
dZdy = gy(:);

%% Constants
g     = 9.81;
xi    = 200;
vis   = 50;
ty    = 2000;
T     = 10;

%% Material point quantities initialization
nmp   = numel(xp);
rhosy = 1100.0;
Msys  = 634.606;
Mp    = Msys/nmp * ones(nmp, 1);
Vp    = Mp./rhosy;
Ap    = Vp./hp;
vp    = zeros (nmp,2);
BINGHAM = 1.0;
FRICTION = 1.0;
CFL = 0.1;
BC_FLAG = 0.0;
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
	   "dZdy", dZdy,
     "BINGHAM_ON", BINGHAM,
     "FRICTION_ON", FRICTION,
     "CFL", CFL,
     "BC_FLAG",BC_FLAG
	 );
json = jsonencode(DATA);

FID = fopen("DATA.json","w");
fprintf(FID,json);
fclose(FID);









