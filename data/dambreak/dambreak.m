
clear all;
close all;
clc;



ndivrows = 20;
ndivcols = 60;

hx = .5;
hy = .5;

y = linspace (0, hy*ndivrows, ndivrows+1);
x = linspace (0, hx*ndivcols, ndivcols+1);
[X, Y] = meshgrid (x, y);

Z = 0.0*X;


##
##h = 0*Z;
##Ymin = 0; Ymax = 10; Xmin = 0; Xmax = 15;
##select = (Y>=Ymin & Y<=Ymax & X >=Xmin & X < Xmax);
##
##h(select)= 10 - Z(select);
##h(h < 0) = 0;




DX = 0.6; DY = 0.6;
[xp, yp] = meshgrid (0:DX:15, 0:DY:10);

%hp = interp2 (X, Y, h, xp, yp, 'linear');
hp = 7.*ones(size(xp));
hp = hp(:);
xp = xp(:);
yp = yp(:);

surf (X, Y, Z)
hold all
scatter3(xp, yp, hp)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

##hzero = find (hp <= 1);
##hp(hzero) = [];
##xp(hzero) = [];
##yp(hzero) = [];

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
T     = 1.2;

%% Material point quantities initialization
nmp   = numel(xp);
rhosy = 1000.0;
Msys  = sum (hp*DX*DY*rhosy);
Mp    = Msys/nmp * ones(nmp, 1);
Vp    = Mp./rhosy;
Ap    = Vp./hp;
vp    = zeros (nmp,2);
BINGHAM = 0.0;
FRICTION = 0.0;
CFL = 0.01;
BC_FLAG = 1.0;

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









