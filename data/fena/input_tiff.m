% clear all; clc; close all;
% A = imread('pr2.tif');
% A = A(:,:,1);
% C = imbinarize(A);
% Bbin = bwboundaries(C);
% 
% for k = 1:length(Bbin)
%    boundary = Bbin{k};
% end
% 
% figure()
% plot(boundary(:,2), boundary(:,1), 'rx-', 'LineWidth', 1.0)
% 
% Lx = 500; Ly = 400;
% nEx = 0.2*Lx; 
% nEy = 0.2*Ly;
% hl = [Lx/nEx Ly/nEy];
% nl=2;
% [xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy; 
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% xL      = 0+(0.5*hl/nl):hl/nl:Lx;
% yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
% [xp,yp] = meshgrid(xL,yL);
% xp = xp(:);
% yp = yp(:);
% [in,on] = inpolygon(xp(:),yp(:),boundary(:,2),boundary(:,1));
% nbound = numel(boundary(:,2));
% xP = [ boundary(:,2); xp(in)];
% yP = [boundary(:,1); yp(in)];
% % xP = xp(in);
% % yP = yp(in);
% hP = 15.*sin(xP(:)./150);
% hi = hP;
% %  figure()
% % %  scatter(xp(:),yp(:),10,'r','filled');
% % scatter(x(:),y(:),13,'r+')
% % hold on
% hold on
% scatter3(xP,yP,hP,7,'b','filled');
% axis equal
% hold on
% scatter(x(:),y(:),13,'r+')
% xp = xP(:);
% yp = yP(:);
% hp = hP(:);
% Xp = [xp(:) yp(:)];
%   figure()
%  shp = alphaShape(Xp(:,1),Xp(:,2));
% % Fplot =  plot(shp);
% %  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
% %  axis([0 size(C,2) 0 size(C,1)])
%  Asys = area(shp);
% 
% nmp = length(xp(:));
% rhosy = 1400.0;
% Msys = (Asys.*mean(hp))*rhosy;
% Mpar = Msys/nmp;
% Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
% %Vpu = (Apu.*mean(hp))./nmp;
% Vpu = (Apu.*mean(hp));
%  boundary_central = find(abs(x-10)<=1e-5 & (y>=8 | y<=2));
% 

%%




% % Domain
% clear; clc;
% % Lx = 15; Ly = 15;
% Lx = 60; Ly = 15;
% nEx = 96; nEy = 24;
% hl = [Lx/nEx Ly/nEy];
% 
% [xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy; 
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% nl = 4;
% nr =3;
% 
% dam = 10;
% % nmp = 4000;
% % hup = 3.5;
% % hdown = 2;
% % Msys = (7.5*15*hup)*1000+(7.5*15*hdown)*1000;
% % Mpar = Msys/nmp;
% % Apu = Mpar./(1e3*hup); Apd = Mpar./(1e3*hdown);
% % At = (Lx/2)*Ly;
% 
% 
% xL      = 0+(0.5*hl/nl):hl/nl:dam;
% yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
% [xpl,ypl] = meshgrid(xL,yL);
% 
% xR      = dam+(0.5*hl/nr):hl/nr:Lx;
% yR      = 0+(0.5*hl/nr):hl/nr: Ly-(0.5*hl/nr);
% [xpr,ypr] = meshgrid(xR,yR);
% 
% 
% xp = [xpl(:); xpr(:)];
% yp = [ypl(:); ypr(:)];
% 
% Xp = [xp(:),yp(:)];
% 
% hup = 5;
% hdown = 1.5;
% hp =hdown*ones(size(xp(:))); hp(xp(:)>=0 & xp(:)<=dam) =hup;
% % hp(sqrt((xp(:)-2.5).^2+(yp(:)-2.5).^2)<=0.5 ) = 2;
% % Ap = [Apl;Apr];
% rhosy = 1000.;
% nmp = length(xp(:));
% nmpup = length(xpl(:));
% nmpdo = length(xpr(:));
% Msys = (dam*Ly*hup)*1000+((Lx-dam)*Ly*hdown)*1000;
% Mpar = Msys/nmp;
% Apu = Mpar./(1e3*hup)*ones(nmpup,1); Apd = Mpar./(1e3*hdown)*ones(nmpdo,1);
% At = (Lx/2)*Ly;
% Vpu = (dam*15*hup)/nmpup*ones(nmpup,1);
% Vpd = ((Lx-dam)*Ly*hdown)/nmpdo*ones(nmpdo,1);
% 
% boundary_central = find(abs(x-10)<=1e-5 & (y<=5 | y>=10) );
% 
% 
% scatter(x(:),y(:),13,'r+')
% hold on
% scatter3(Xp(:,1),Xp(:,2),hp(:),5,'b','filled')
% hold on
% axis equal

%% GOCCIA CHE CADE
% clear all; clc; close all;
% 
% 
% Lx = 5; Ly = 5;
% nEx = 70;
% nEy = 70;
% hl = [Lx/nEx Ly/nEy];
% nl=6;
% [xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
% NVx = size(xv,2);
% NVy = size(yv,1);
% NV  = NVx*NVy;
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% xL      = 0+(0.5*hl(1)/nl):hl(1)/nl:Lx;
% yL      = 0+(0.5*hl(2)/nl):hl(2)/nl: Ly;
% [xpp,ypp] = meshgrid(xL,yL);
% xp = xpp(:);
% yp = ypp(:);
% 
% hp =  ones(size(xp(:),1),1);
%  hp(sqrt((xp-2.5).^2+(yp-2.5).^2)<0.5) = 2;
% % hp(xp<=10) = 2;
% xp(hp<2) = [];
% yp(hp<2) = [];
% hp(hp<2) = [];
% Xp = [xp(:) yp(:)];
%   figure()
% scatter(x(:),y(:),13,'r+')
% hold on
% 
% scatter3(xp,yp,hp,7,'b','filled');
% axis equal
% %   figure()
%  shp = alphaShape(Xp(:,1),Xp(:,2));
% % Fplot =  plot(shp);
% %  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
% %  axis([0 size(C,2) 0 size(C,1)])
%  Asys = area(shp);
% % boundary_central = find(abs(x-10)<=1e-5 & (y>=8 | y<=2));
% nmp = length(xp(:));
% rhosy = 1000;
% Msys = (Asys.*mean(hp))*rhosy;
% Mpar = Msys/nmp;
% Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
% Vpu = (Apu.*mean(hp)); Vpd = [];%(dam*15*hp(1))/nmp*ones(nmp,1);
%% TIFF EXAMPLE

% clear all; clc; close all;
% A = imread('f5.tif');
% A = A(:,:,1);
% C = imbinarize(A);
% Bbin = bwboundaries(C);
% 
% for k = 1:length(Bbin)
%    boundary = Bbin{k};
% end
% 
% figure()
% plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
% 
% Lx = 1000; Ly = 350;
% nEx = 0.4*Lx; 
% nEy = 0.4*Ly;
% hl = [Lx/nEx Ly/nEy];
% nl=2;
% [xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy; 
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% xL      = 0+(0.5*hl/nl):hl/nl:Lx;
% yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
% [xp,yp] = meshgrid(xL,yL);
% xp = xp(:);
% yp = yp(:);
% [in,on] = inpolygon(xp(:),yp(:),boundary(:,2),boundary(:,1));
% 
% % xP = [xp(in); xp(on)];
% % yP = [yp(in); yp(on)];
% xP = xp(in);
% yP = yp(in);
% hP = 15.*sin(xP(:)./150);
% hi = hP;
% %  figure()
% % %  scatter(xp(:),yp(:),10,'r','filled');
% % scatter(x(:),y(:),13,'r+')
% % hold on
% % 
% % scatter3(xP,yP,hP,7,'b','filled');
% % axis equal
% 
% xp = xP(:);
% yp = yP(:);
% hp = hP(:);
% % Setting initial area, mass, density
% 
% Xp = [xp(:) yp(:)];
%  
%   figure()
%  shp = alphaShape(Xp(:,1),Xp(:,2));
% % Fplot =  plot(shp);
% %  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
% %  axis([0 size(C,2) 0 size(C,1)])
%  Asys = area(shp);
% 
% nmp = length(xp(:));
% rhosy = 1200;
% Msys = (Asys.*mean(hp))*rhosy;
% Mpar = Msys/nmp;
% Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
% %Vpu = (Apu.*mean(hp))./nmp;
% Vpu = (Apu.*mean(hp));

%% SLIDING BENCHMARK with exponential obstl
% clear all;
% clc;
% 
% Lx = 150; Ly = 20;
% nEx = 300; %180
% nEy = 40;
% hl = [Lx/nEx Ly/nEy];
% nl=3;
% [xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy; 
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% xL      = 0+(0.5*hl/nl):hl/nl:Lx;
% yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
% [xpl,ypl] = meshgrid(xL,yL);
% L = Lx;
% % omega = (((xpl-L/2).^2 + (ypl-L/2).^2)/L/L) <= (.2 + .01 * sin(10*pi*(ypl-L/2)/L)).^2;
% %omega = (((x-L/2).^2 + (y-L/2).^2)) <= (.2 + .01 * sin(10*pi*(y-L/2))).^2;
% %hp = max (0, min ((60-(100 - 100 * xpl/L)), 2)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega ./5; 
% 
% % xR      = dam+(0.5*hl/nr):hl/nr:Lx;
% % yR      = 0+(0.5*hl/nr):hl/nr: Ly-(0.5*hl/nr);
% % [xpr,ypr] = meshgrid(xR,yR);
% xpr = [];
% ypr = xpr;
% 
% xp = [xpl(:); xpr(:)];
% yp = [ypl(:); ypr(:)];
%  hp =  1.*ones(size(xp(:),1),1);
% %   Z = 3+2.*sin(x).*cos(0.5.*y);
% % dZdx = 2.*cos(x).*cos(0.5.*y);
% % dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);
%  Z =60 -0.30.*x + 7.*exp(-(x-50).^2./7 - (y-10).^2./7); % 3+2.*sin(x).*cos(0.5.*y);
% dZdx = -0.30 - 7/7.*2.*(x-50).*exp(-(x-50).^2./7 - (y-10).^2./7); %2.*cos(x).*cos(0.5.*y);
% dZdy = - 7/7.*2.*(y-10).*exp(-(x-50).^2./7 - (y-10).^2./7);
% %    hp(sqrt((xp-3).^2+(yp-5).^2)<1) = 1.3;
% %   hp(((xp-10).^2./3+(yp-10).^2./5 - sqrt((xp-12).^2 + (yp-12).^2)./4) <0.00001) = 1.5;
%  hp(xp<=30 & yp>=0) = 2;
% % xp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% % yp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% % hp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% xp(hp<=1)=[];
% yp(hp<=1)=[];
% hp(hp<=1)=[];
% % xp(hp<=2)=[];
% % yp(hp<=2)=[];
% % hp(hp<=2)=[];
% 
% Xp = [xp(:),yp(:)];
% Zplot =60 -0.30.*xv + 7.*exp(-(xv-50).^2./7 - (yv-10).^2./7); %2.*sin(xv).*cos(0.5.*yv);
% 
% 
% 
%  figure()
% % scatter(x(:),y(:),13,'r+')
% % hold on
% surf(xv,yv,Zplot)
% hold on
% scatter3(xp,yp,hp,7,'b','filled');
% axis equal
% xlabel('x');
% ylabel('y');
%   figure()
%  shp = alphaShape(Xp(:,1),Xp(:,2),1);
% % Fplot =  plot(shp);
% %  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
% %  axis([0 size(C,2) 0 size(C,1)])
%  Asys = area(shp);
% 
% nmp = length(xp(:));
% rhosy = 1000.0;
% Msys = (Asys.*mean(hp))*rhosy;
% Mpar = Msys/nmp;
% Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
% %Vpu = (Apu.*mean(hp))./nmp;
% Vpu = (Apu.*mean(hp));
%  boundary_central = find(abs(x-10)<=1e-5 & (y>=8 | y<=2));
% 
% X = x(:);
% Y = y(:);
%  Z = 3+2.*sin(x).*cos(0.5.*y);
% dZdx = 2.*cos(x).*cos(0.5.*y);
% dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);

%% SIN/COS DOMAIN

% clear all;
% clc;
% 
% Lx = 40; Ly = 40;
% nEx = 160; %180
% nEy = 160;
% hl = [Lx/nEx Ly/nEy];
% nl=3;
% [xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy; 
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% xL      = 0+(0.5*hl/nl):hl/nl:Lx;
% yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
% [xpl,ypl] = meshgrid(xL,yL);
% L = Lx;
% % omega = (((xpl-L/2).^2 + (ypl-L/2).^2)/L/L) <= (.2 + .01 * sin(10*pi*(ypl-L/2)/L)).^2;
% %omega = (((x-L/2).^2 + (y-L/2).^2)) <= (.2 + .01 * sin(10*pi*(y-L/2))).^2;
% %hp = max (0, min ((60-(100 - 100 * xpl/L)), 2)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega ./5; 
% 
% % xR      = dam+(0.5*hl/nr):hl/nr:Lx;
% % yR      = 0+(0.5*hl/nr):hl/nr: Ly-(0.5*hl/nr);
% % [xpr,ypr] = meshgrid(xR,yR);
% xpr = [];
% ypr = xpr;
% 
% xp = [xpl(:); xpr(:)];
% yp = [ypl(:); ypr(:)];
%  hp =  1.*ones(size(xp(:),1),1);
% %    Z = 3+2.*sin(x).*cos(0.5.*y);
% %  dZdx = 2.*cos(x).*cos(0.5.*y);
% %  dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);
% 
% %Z =60 -0.30.*x + 7.*exp(-(x-50).^2./7 - (y-10).^2./7); % 3+2.*sin(x).*cos(0.5.*y);
% %dZdx = -0.30 - 7/7.*2.*(x-50).*exp(-(x-50).^2./7 - (y-10).^2./7); %2.*cos(x).*cos(0.5.*y);
% %dZdy = - 7/7.*2.*(y-10).*exp(-(x-50).^2./7 - (y-10).^2./7);
% %    hp(sqrt((xp-3).^2+(yp-5).^2)<1) = 1.3;
%    hp(((xp-15).^2./3+(yp-15).^2./5 - sqrt((xp-20).^2 + (yp-20).^2)./4) <0.00001) = 2;
% % hp(xp<=30 & yp>=0) = 2;
% % xp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% % yp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% % hp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% xp(hp<=1)=[];
% yp(hp<=1)=[];
% hp(hp<=1)=[];
% % xp(hp<=2)=[];
% % yp(hp<=2)=[];
% % hp(hp<=2)=[];
% 
% Xp = [xp(:),yp(:)];
% Zplot =3+2.*sin(xv).*cos(0.5.*yv);
% [gx,gy] = gradient(Zplot);
% 
% Z = reshape(Zplot,NVy*NVx,1);
% dZdx = reshape(gx,NV,1);
% dZdy = reshape(gy,NV,1);
% 
% 
%  figure()
% % scatter(x(:),y(:),13,'r+')
% % hold on
% surf(xv,yv,Zplot)
% hold on
% scatter3(xp,yp,hp,7,'b','filled');
% axis equal
% xlabel('x');
% ylabel('y');
%   figure()
%  shp = alphaShape(Xp(:,1),Xp(:,2),1);
% % Fplot =  plot(shp);
% %  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
% %  axis([0 size(C,2) 0 size(C,1)])
%  Asys = area(shp);
% 
% nmp = length(xp(:));
% rhosy = 1000.0;
% Msys = (Asys.*mean(hp))*rhosy;
% Mpar = Msys/nmp;
% Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
% %Vpu = (Apu.*mean(hp))./nmp;
% Vpu = (Apu.*mean(hp));
%  boundary_central = find(abs(x-10)<=1e-5 & (y>=8 | y<=2));
% 
% X = x(:);
% Y = y(:);
%  Z = 3+2.*sin(x).*cos(0.5.*y);
% dZdx = 2.*cos(x).*cos(0.5.*y);
% dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);

 %% VOLCANO
% 
% clear all;
% clc;
% 
% Lx = 250; Ly = 250;
% nEx = 149; %180
% nEy = 149;
% hl = [Lx/nEx Ly/nEy];
% nl=5;
% [xv,yv] = meshgrid(0:hl(1):Lx,0:hl(2):Ly);
% %%
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy; 
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% xL      = 0+(0.5*hl/nl):hl/nl:Lx;
% yL      = 0+(0.5*hl/nl):hl/nl: Ly-(0.5*hl/nl);
% [xpl,ypl] = meshgrid(xL,yL);
% L = Lx;
% % omega = (((xpl-L/2).^2 + (ypl-L/2).^2)/L/L) <= (.2 + .01 * sin(10*pi*(ypl-L/2)/L)).^2;
% %omega = (((x-L/2).^2 + (y-L/2).^2)) <= (.2 + .01 * sin(10*pi*(y-L/2))).^2;
% %hp = max (0, min ((60-(100 - 100 * xpl/L)), 2)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega ./5; 
% 
% % xR      = dam+(0.5*hl/nr):hl/nr:Lx;
% % yR      = 0+(0.5*hl/nr):hl/nr: Ly-(0.5*hl/nr);
% % [xpr,ypr] = meshgrid(xR,yR);
% xpr = [];
% ypr = xpr;
% 
% xp = [xpl(:); xpr(:)];
% yp = [ypl(:); ypr(:)];
%  hp =  1.*ones(size(xp(:),1),1);
% %    Z = 3+2.*sin(x).*cos(0.5.*y);
% %  dZdx = 2.*cos(x).*cos(0.5.*y);
% %  dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);
% 
% %Z =60 -0.30.*x + 7.*exp(-(x-50).^2./7 - (y-10).^2./7); % 3+2.*sin(x).*cos(0.5.*y);
% %dZdx = -0.30 - 7/7.*2.*(x-50).*exp(-(x-50).^2./7 - (y-10).^2./7); %2.*cos(x).*cos(0.5.*y);
% %dZdy = - 7/7.*2.*(y-10).*exp(-(x-50).^2./7 - (y-10).^2./7);
% %    hp(sqrt((xp-3).^2+(yp-5).^2)<1) = 1.3;
% 
%    hp(((xp-92).^2./3+(yp-83).^2./5) <200*(.2 + .01 * sin(37*pi*(yp-L/2)/L))) = 2;
% 
% %    hp(((xp-.55).^2./3+(yp-.7).^2./5) <0.001*(.2 + .01 * sin(150*pi*(yp-L/2)/L)).^2) = 1.1;
% % hp(((xp-800).^2./3+(yp-1020).^2./5) <4000*(.2 + .01 * sin(150*pi*(yp))).^2) = 17;
% % hp(xp<=30 & yp>=0) = 2;
% % xp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% % yp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% % hp(hp<=min ((( 500 * xpl/L-200)), 30))=[];
% xp(hp<=1)=[];
% yp(hp<=1)=[];
% hp(hp<=1)=[];
% % xp(hp<=2)=[];
% % yp(hp<=2)=[];
% % hp(hp<=2)=[];
% 
% Xp = [xp(:),yp(:)];
% %load volcano_data.mat
% load HH.mat
% heights = 0.65.*HH(151:300,151:300);
% H = reshape(heights,NVy,NVx);
% Zplot =H; %3+2.*sin(xv).*cos(0.5.*yv);
% Zplot2 = 3+2.*sin(xv).*cos(0.5.*yv);
% [gx,gy] = gradient(Zplot);
% 
% Z = reshape(Zplot,NVy*NVx,1);
% dZdx = reshape(gx,NV,1);
% dZdy = reshape(gy,NV,1);
% 
% 
%  figure()
% % scatter(x(:),y(:),13,'r+')
% % hold on
% surf(xv,yv,Zplot)
% hold on
% scatter3(xp,yp,hp,7,'b','filled');
% axis([0 Lx 0 Ly])
% xlabel('x');
% ylabel('y');
%   figure()
%  shp = alphaShape(Xp(:,1),Xp(:,2),1);
% % Fplot =  plot(shp);
% %  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
% %  axis([0 size(C,2) 0 size(C,1)])
%  Asys = area(shp);
% 
% nmp = length(xp(:));
% rhosy = 1291.0;
% Msys = (Asys.*mean(hp))*rhosy;
% Mpar = Msys/nmp;
% Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
% %Vpu = (Apu.*mean(hp))./nmp;
% Vpu = (Apu.*mean(hp));
%  boundary_central = find(abs(x-10)<=1e-5 & (y>=8 | y<=2));
% 
% X = x(:);
% Y = y(:);
% %  Z = 3+2.*sin(x).*cos(0.5.*y);
% % dZdx = 2.*cos(x).*cos(0.5.*y);
% % dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);
% 

%% LIGURIA


clear all;
clc;

Lx = 259; Ly = 250;
nEx = 259; %180
nEy = 250;
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

xL      = 0+(0.5*hl/nl):hl/nl:259;
yL      = 0+(0.5*hl/nl):hl/nl: 250-(0.5*hl/nl);
[xpl,ypl] = meshgrid(xL,yL);
L = Lx;
% omega = (((xpl-L/2).^2 + (ypl-L/2).^2)/L/L) <= (.2 + .01 * sin(10*pi*(ypl-L/2)/L)).^2;
%omega = (((x-L/2).^2 + (y-L/2).^2)) <= (.2 + .01 * sin(10*pi*(y-L/2))).^2;
%hp = max (0, min ((60-(100 - 100 * xpl/L)), 2)) .* omega; 
%hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega; 
%hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega ./5; 

% xR      = dam+(0.5*hl/nr):hl/nr:Lx;
% yR      = 0+(0.5*hl/nr):hl/nr: Ly-(0.5*hl/nr);
% [xpr,ypr] = meshgrid(xR,yR);
xpr = [];
ypr = xpr;

xp = [xpl(:); xpr(:)];
yp = [ypl(:); ypr(:)];
 hp =  1.*ones(size(xp(:),1),1);
%    Z = 3+2.*sin(x).*cos(0.5.*y);
%  dZdx = 2.*cos(x).*cos(0.5.*y);
%  dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);

%Z =60 -0.30.*x + 7.*exp(-(x-50).^2./7 - (y-10).^2./7); % 3+2.*sin(x).*cos(0.5.*y);
%dZdx = -0.30 - 7/7.*2.*(x-50).*exp(-(x-50).^2./7 - (y-10).^2./7); %2.*cos(x).*cos(0.5.*y);
%dZdy = - 7/7.*2.*(y-10).*exp(-(x-50).^2./7 - (y-10).^2./7);
%    hp(sqrt((xp-3).^2+(yp-5).^2)<1) = 1.3;

%    hp(((xp-92).^2./3+(yp-83).^2./5) <200*(.2 + .01 * sin(37*pi*(yp-L/2)/L))) = 2;

  % hp(((xp-.55).^2./3+(yp-.7).^2./5) <0.001*(.2 + .01 * sin(150*pi*(yp-L/2)/L)).^2) = 17;
 hp(((xp-154).^2./3+(yp-210).^2./5) <350*(.2 + .01 * sin(150*pi*(yp))./100).^2) = 17;

xp(hp<=1)=[];
yp(hp<=1)=[];
hp(hp<=1)=[];


Xp = [xp(:),yp(:)];
rowstart = 550;
rowend   = 800;

colstart = 171;
colend   = 430;

ndivrows = rowend - rowstart;
ndivcols = colend - colstart;
%load volcano_data.mat
immagine = imread ("dtm_liguria_2017.tif", "PixelRegion", {[rowstart rowend], [colstart, colend]});
IMM = double(immagine(end:-1:1, :));
% immagine_int = reshape (typecast (immagine(end:-1:1, :), "int16"), size (immagine));
immagine_int = reshape (IMM, size (immagine));
immagine_int(immagine_int<0) = -1;

Zz = double (immagine_int);
Zz(Zz<=0) = 0.0;


[gx, gy] = gradient(Zz);

Z  = Zz(:);
dZdx = gx(:); %reshape(gx,NV,1);
dZdy = gy(:); %reshape(gy,NV,1);




 figure()
% scatter(x(:),y(:),13,'r+')
% hold on
mesh(xv,yv,Zz)
hold on
scatter3(xp,yp,hp,7,'b','filled');
axis([0 Lx 0 Ly])
xlabel('x');
ylabel('y');
  figure()
 shp = alphaShape(Xp(:,1),Xp(:,2),1);
% Fplot =  plot(shp);
%  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
%  axis([0 size(C,2) 0 size(C,1)])
 Asys = area(shp);

nmp = length(xp(:));
rhosy = 1200.0;
Msys = (Asys.*mean(hp))*rhosy;
Mpar = Msys/nmp;
Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
%Vpu = (Apu.*mean(hp))./nmp;
Vpu = (Apu.*mean(hp));
 boundary_central = ones(nmp,1); %find(abs(x-10)<=1e-5 & (y>=8 | y<=2));

X = x(:);
Y = y(:);
%  Z = 3+2.*sin(x).*cos(0.5.*y);
% dZdx = 2.*cos(x).*cos(0.5.*y);
% dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);


%% LIGURIA 2

% 
% clear all;
% clc;
% rowstart = 550; %550;
% rowend   = 800; %800;
% 
% colstart = 171; %171;
% colend   = 430; %430;
% Lx = colend-1; Ly = rowend-1;
% nEx = Lx; %180
% nEy = Ly;
% hl = [Lx/nEx Ly/nEy];
% nl=2;
% [xv,yv] = meshgrid(colstart:hl(1):Lx+1,rowstart:hl(2):Ly+1);
% 
% NVx = size(xv,2);                                                     
% NVy = size(yv,1);                                                      
% NV  = NVx*NVy; 
% NP = (NVx-1)*(NVy-1);
% Nex = NVx-1;
% Ney = NVy-1;
% x = xv(:); y = yv(:);
% 
% xL      = colstart+(0.5*hl/nl):hl/nl:Lx+1;
% yL      = rowstart+(0.5*hl/nl):hl/nl: Ly+1-(0.5*hl/nl);
% [xpl,ypl] = meshgrid(xL,yL);
% L = Lx;
% % omega = (((xpl-L/2).^2 + (ypl-L/2).^2)/L/L) <= (.2 + .01 * sin(10*pi*(ypl-L/2)/L)).^2;
% %omega = (((x-L/2).^2 + (y-L/2).^2)) <= (.2 + .01 * sin(10*pi*(y-L/2))).^2;
% %hp = max (0, min ((60-(100 - 100 * xpl/L)), 2)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega; 
% %hp = max (0, min ((( 500 * xpl/L-200)), 30)) .* omega ./5; 
% 
% % xR      = dam+(0.5*hl/nr):hl/nr:Lx;
% % yR      = 0+(0.5*hl/nr):hl/nr: Ly-(0.5*hl/nr);
% % [xpr,ypr] = meshgrid(xR,yR);
% xpr = [];
% ypr = xpr;
% 
% xp = [xpl(:); xpr(:)];
% yp = [ypl(:); ypr(:)];
%  hp =  1.*ones(size(xp(:),1),1);
% %    Z = 3+2.*sin(x).*cos(0.5.*y);
% %  dZdx = 2.*cos(x).*cos(0.5.*y);
% %  dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);
% 
% %Z =60 -0.30.*x + 7.*exp(-(x-50).^2./7 - (y-10).^2./7); % 3+2.*sin(x).*cos(0.5.*y);
% %dZdx = -0.30 - 7/7.*2.*(x-50).*exp(-(x-50).^2./7 - (y-10).^2./7); %2.*cos(x).*cos(0.5.*y);
% %dZdy = - 7/7.*2.*(y-10).*exp(-(x-50).^2./7 - (y-10).^2./7);
% %    hp(sqrt((xp-3).^2+(yp-5).^2)<1) = 1.3;
% 
% %    hp(((xp-92).^2./3+(yp-83).^2./5) <200*(.2 + .01 * sin(37*pi*(yp-L/2)/L))) = 2;
% 
%   % hp(((xp-.55).^2./3+(yp-.7).^2./5) <0.001*(.2 + .01 * sin(150*pi*(yp-L/2)/L)).^2) = 17;
%  hp(((xp-330).^2./3+(yp-750).^2./5) <100*(.2 + .01 * sin(150*pi*(yp))).^2) = 17;
% % hp(xp<330 & xp>325 & yp>736 & yp<760) = 17;
% xp(hp<=1)=[];
% yp(hp<=1)=[];
% hp(hp<=1)=[];
% 
% 
% Xp = [xp(:),yp(:)];
% 
% 
% ndivrows = rowend - rowstart;
% ndivcols = colend - colstart;
% %load volcano_data.mat
% immagine = imread ("dtm_liguria_2017.tif", "PixelRegion", {[rowstart rowend], [colstart, colend]});
% IMM = double(immagine(end:-1:1, :));
% % immagine_int = reshape (typecast (immagine(end:-1:1, :), "int16"), size (immagine));
% immagine_int = reshape (IMM, size (immagine));
% immagine_int(immagine_int<0) = -1;
% 
% Zz = double (immagine_int);
% Zz(Zz<=0) = 0.0;
% 
% 
% [gx, gy] = gradient(Zz);
% 
% Z  = Zz(:);
% dZdx = gx(:); %reshape(gx,NV,1);
% dZdy = gy(:); %reshape(gy,NV,1);
% 
% 
% 
% 
%  figure()
% % scatter(x(:),y(:),13,'r+')
% % hold on
% mesh(xv,yv,Zz)
% hold on
% scatter3(xp,yp,hp,7,'b','filled');
% axis([colstart Lx rowstart Ly])
% xlabel('x');
% ylabel('y');
%   figure()
%  shp = alphaShape(Xp(:,1),Xp(:,2),1);
% % Fplot =  plot(shp);
% %  set(Fplot,'FaceColor','red','FaceLighting','gouraud','EdgeColor','none')
% %  axis([0 size(C,2) 0 size(C,1)])
%  Asys = area(shp);
% 
% nmp = length(xp(:));
% rhosy = 1200.0;
% Msys = (Asys.*mean(hp))*rhosy;
% Mpar = Msys/nmp;
% Apu =Mpar./(rhosy*mean(hp))*ones(nmp,1); Apd =[];
% %Vpu = (Apu.*mean(hp))./nmp;
% Vpu = (Apu.*mean(hp));
%  boundary_central = ones(nmp,1); %find(abs(x-10)<=1e-5 & (y>=8 | y<=2));
% 
% X = x(:);
% Y = y(:);
% %  Z = 3+2.*sin(x).*cos(0.5.*y);
% % dZdx = 2.*cos(x).*cos(0.5.*y);
% % dZdy = -2.*0.5.*sin(x).*sin(0.5.*y);




               
%% Constants
g = 9.81;
t = 0;                  
T = 1.5;                  
%% Grid quantities initialization
Mv = zeros(NV,1);
MvR = zeros(2*NV,1);
vv = zeros(2*NV,1);
Uhv = zeros(2*NV,1);

dv = 0.0;

%% Material point quantities initialization
yd = 4;
% rhosy = 1400.0;
slo = 1/sqrt(3);
% Z = 500-slo*x;
hp = reshape(hp,[nmp 1]);
Fp   = zeros(nmp,4) + repmat([1 0 0 1],nmp,1);
rho = rhosy*ones(nmp,1);
vp  = zeros(nmp,2); 
Mp =  Mpar*ones(nmp,1);%Vp.*rho;
Vp = [Vpu];%;Vpd];
Ap = [Apu];%Apd];
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
xi = 200;
Fb(:,1) = zeros(nmp,1);
Fb(:,2) = zeros(nmp,1);

it = 1;
% Xplot(:,:,:,it) = [Xp(:,1) Xp(:,2) hp(:,1)];
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

meanZ = ones(nmp,1);


%%
DATA = struct("x",xp,"y",yp,"Mp",Mp,"Ap",Ap,"vpx",vp(:,1),"vpy",vp(:,2),"Nex",nEx,"Ney",nEy,"hx",hl(1),"hy",hl(2),"hp",hp,"mom_px",momp(:,1),...
       "mom_py",momp(:,2),"g",g,"T",T,"phi",phi,"xi",xi,"vis",vis,"ty",ty,"rho",rhosy,"Vp",Vp,"Z",Z,"dZdx",dZdx,"dZdy",dZdy);
json = jsonencode(DATA)

FID = fopen("DATA.json","w");
fprintf(FID,json);
fclose(FID);









