
clear all;
close all;
clc;

load('mask_in_vladi.mat')
load('dem_real.mat')
% load('mask_out.tif');

PP = zeros(676,676);
pm = zeros(175,175);

for i = 1:175
    pm(176-i,i)=1;
end

mask_in = pm*mask_in;

for i = 1:676


    PP(677-i,i) = 1;

end

dem = PP*dem;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowstart = 1;
rowend   = 676;

colstart =1;
colend   =756;

ndivrows = rowend - rowstart;
ndivcols = colend - colstart;

Xmin = 180;
Xmax = 520;
Ymin = 180;
Ymax = 520;

hx =5;
hy = 5;

y = linspace (0, hy*ndivrows, ndivrows+1);
x = linspace (0, hx*ndivcols, ndivcols+1);

%y = linspace (0, hy*(Ymax-Ymin), Ymax-Ymin );
%x = linspace (0, hx*(Xmax-Xmin), Xmax - Xmin );

[X, Y] = meshgrid (x, y);
%Z = double (immagine_int);
Z = dem; %double(dem(Xmin:Xmax,Ymin:Ymax));

zplot = dem; % double( dem(Xmin:Xmax,Ymin:Ymax));



%h = 0*Z; %ymin 1400 ymax=1700   xmin = 700 xmax = 900
%{
%Ymin = 1350; Ymax = 1450; Xmin = 700; Xmax = 900;
select = (Y>Ymin & Y<Ymax & X >Xmin & X < Xmax);

h(select)= 38.*double(mask_in(select)) ; %(700.*double(mask_in(select))) - Z(select);
h(h < 0) = 0;
%}


%h = 38.*double(mask_in);

%h(h < 0) = 0;
% DX = 1; DY = 1;
% %[xp, yp] = meshgrid (Xmin:DX:Xmax, Ymin:DY:Ymax);
% %hp = interp2 (X, Y, h, xp, yp, 'spline');
% %select = (Y>Ymin & Y<Ymax & X >Xmin & X < Xmax);
% 
% 
% 
% 
% 
% hp = 38.*double(mask_in);
% 
% %yy = linspace (0, hy*ndivrows, ndivrows+1);
% %xx = linspace (0, hx*ndivcols, ndivcols+1);
% 
% yy = linspace (0, hy*size(hp,1), size(hp,1));
% xx = linspace (0, hx*size(hp,2), size(hp,2));
% 
% [Xx, Yy] = meshgrid (xx, yy);
% xp = Xx;
% yp = Yy;
% 
% hp = hp(:);
% xp = xp(:);
% yp = yp(:);
% 
% hzero = find (hp < 38);
% hp(hzero) = [];
% xp(hzero) = [];
% yp(hzero) = [];

%data = jsonencode (struct ('X', X(:), 'Y', Y(:), 'Z', Z(:), 'h', h(:)));
%immagine_out = uint16(immagine_int);
%imwrite (immagine_out, "sfondo.tif")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA = mask_in; %double(mask_in);
AA = AA(:,:,1);
CC = imbinarize(AA);
Bbin = bwboundaries(CC);

for k = 1:length(Bbin)
   boundary = Bbin{k};
end

figure()
plot(boundary(:,2), boundary(:,1), 'rx-', 'LineWidth', 1.0)

% xL      = 0+(0.5*hl/nl):hl/nl:Lx;
% yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
yy = linspace (0, 5*size(mask_in,1), 3650);
xx = linspace (0, 5*size(mask_in,2), 3650);
[xp,yp] = meshgrid(xx,yy);
xp = xp(:);
yp = yp(:);
[in,on] = inpolygon(xp(:),yp(:),boundary(:,2),boundary(:,1));
nbound = numel(boundary(:,2));
xp = 5.*[  xp(in)];
yp = 5.*[ yp(in)];
hp = 38.*ones(numel(xp),1);
DX = 5*size(mask_in,1)/6000;
DY = 5*size(mask_in,2)/6000;
xp = xp + 1530;
yp = yp + 1141;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%figure(2)
%surf(X(100:600,240:700),Y(100:600,240:700),dem(100:600,240:700));
%colormap('hsv');

figure(3)
scatter3(xp,yp,hp,10,'r')
hold on
surf(X,Y,dem,'LineStyle','none')
colormap("jet")

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
T     = 20;

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

BINGHAM = 1.0;
FRICTION = 1.0;
CFL = 0.1;
BC_FLAG = 0.0;

%%
DATA = struct( "x", xp, ...
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
	   "dZdx", dZdx,...
	   "dZdy", dZdy,...
            "BINGHAM_ON", BINGHAM,...
     "FRICTION_ON", FRICTION, ...
     "CFL", CFL,...
     "BC_FLAG",BC_FLAG);
json = jsonencode(DATA);

FID = fopen("DATA.json","w");
fprintf(FID,json);
fclose(FID);









