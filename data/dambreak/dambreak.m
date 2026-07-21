clear all;
close all;
clc;

neley = 60;
nelex = 200;

y = linspace (0, 30, neley+1);
x = linspace (0, 100, nelex+1);

hx = 100/nelex;   % = 0.5
hy = 30/neley;    % = 0.5

[X, Y] = meshgrid (x, y);

Z = ones(size(X));  

% Regione iniziale occupata dall'acqua 
Ymin = 0; Ymax = 30; Xmin = 0; Xmax = 50;

ppc_x = 3;            % particelle per cella in x
ppc_y = 2;            % particelle per cella in y   (=> 4 ppc totali)
DX = hx / ppc_x;      % = 0.5/2 = 0.25
DY = hy / ppc_y;      % = 0.5/2 = 0.25

% posizioni ai centri delle sotto-celle: da Xmin+DX/2 a Xmax-DX/2
xp_vec = (Xmin + DX/2) : DX : (Xmax - DX/2);
yp_vec = (Ymin + DY/2) : DY : (Ymax - DY/2);
[xp, yp] = meshgrid (xp_vec, yp_vec);

% Altezza iniziale uniforme
hp = 5 .* ones(size(xp));

hp = hp(:);
xp = xp(:);
yp = yp(:);

% (inerte con hp = 5 ovunque; utile se in futuro hp e' variabile)
hzero = find (hp <= 1);
hp(hzero) = [];
xp(hzero) = [];
yp(hzero) = [];

immagine_out = uint16(Z);
imwrite (immagine_out, "sfondo.tif")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zz = double (Z);

[gx, gy] = gradient (Zz, hx, hy);

Z    = Zz(:);
dZdx = gx(:);
dZdy = gy(:);

figure()
mesh(X, Y, Zz)
axis equal
hold on
scatter3(xp(:), yp(:), hp(:))

%% Constants
g     = 9.81;
xi    = 200;
vis   = 50;
ty    = 2000;
T     = 5;

%% Material point quantities initialization
nmp   = numel(xp);
rhosy = 1000.0;
Msys  = sum (hp*DX*DY*rhosy);   % massa totale del sistema
Mp    = Msys/nmp * ones(nmp, 1);
Vp    = Mp./rhosy;
Ap    = Vp./hp;                 % = DX*DY = 0.0625 (con ppc=2)
vp    = zeros (nmp,2);
BINGHAM  = 0.0;
FRICTION = 0.0;
CFL = 0.09;
BC_FLAG = 1.0;
momp  = zeros (nmp,2);

Fb(:,1) = zeros (nmp,1);
Fb(:,2) = zeros (nmp,1);

%% Riepilogo a schermo
fprintf('hx = %g, hy = %g\n', hx, hy);
fprintf('DX = %g, DY = %g  -> %g x %g ppc (interi)\n', DX, DY, hx/DX, hy/DY);
fprintf('nmp = %d particelle\n', nmp);
fprintf('Ap iniziale = %g (=> aggiorna Ap0 nel main.cpp se serve)\n', Ap(1));
fprintf('Msys = %g kg\n', Msys);

%%
DATA = struct (...
	   "x", xp, ...
	   "y", yp, ...
	   "Mp", Mp, ...
	   "Ap", Ap, ...
	   "vpx", vp(:,1), ...
	   "vpy", vp(:,2), ...
	   "Nex", nelex, ...
	   "Ney", neley, ...
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
	   "dZdx", dZdx,...
	   "dZdy", dZdy,...
     "BINGHAM_ON", BINGHAM,...
     "FRICTION_ON", FRICTION,...
     "CFL", CFL,...
     "BC_FLAG", BC_FLAG...
	 );
json = jsonencode(DATA);

FID = fopen("DATA.json","w");
fprintf(FID,json);
fclose(FID);