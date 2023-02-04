
/*
   MPM 2D - SWE

		While t<T
(0)   Set to zero grid quantities:
(1)		Particles -> Grid (P2G)
			  (1.0)	 Mass on vertices:
				(1.1)  Momenta on vertices:
				(1.2)  External Forces on vertices:
      End P2G
(2)		Internal Forces on vertices:
(3)		Compute momenta:
(4)		Get nodal accelerations and velocities:
(5)		Enforce boundary conditions (e.g. Dirichlet)
(6)		Particles <- Grid (G2P)
				(6.0)   Update particle velocities:
				(6.1)   Update particle positions:
				(6.2)   Update particle momenta:
      End G2P
(7)   Update particle heights:
(8)   Update particle area and stresses:
		    (8.0)   Update particle fluxes  using the quantities just computed
			  (8.1)   Update stresses (USL approach):
(9)	Check that Masses  and densities  must be invariant throughout the simulation
(10)  Advance time:

		EndWhile
*/
#include <iostream>
#include <quadgrid_cpp.h>
#include <particles.h>
#include "json.hpp"
#include <fstream>
#include <cmath>
#include <map>

struct DATA
{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> Mp;
  std::vector<double> Ap;
  std::vector<double> vpx;
  std::vector<double> vpy;
  int Nex;
  int Ney;
  double hx;
  double hy;
  std::vector<double> hp;
  std::vector<double> mom_px;
  std::vector<double> mom_py;
//  std::vector<double> F_ext_px;//serve divergenza per calcolarle ---> come fare su particella?
//  std::vector<double> F_ext_py;
//  std::vector<double> F_int_px;
//  std::vector<double> F_int_py;
  double g;
  double T;
  double phi;
  double xi;
  double vis;
  double ty;
  double rho;
};

void
from_json(const nlohmann::json &j, DATA &d)
{
  j.at("x").get_to(d.x);
  j.at("y").get_to(d.y);
  j.at("Mp").get_to(d.Mp);
  j.at("Ap").get_to(d.Ap);
  j.at("vpx").get_to(d.vpx);
  j.at("vpy").get_to(d.vpy);
  j.at("Nex").get_to(d.Nex);
  j.at("Ney").get_to(d.Ney);
  j.at("hx").get_to(d.hx);
  j.at("hy").get_to(d.hy);
  j.at("hp").get_to(d.hp);
  j.at("mom_px").get_to(d.mom_px);
  j.at("mom_py").get_to(d.mom_py);
  j.at("rho").get_to(d.rho);
  j.at("xi").get_to(d.xi);
  j.at("ty").get_to(d.ty);
  j.at("vis").get_to(d.vis);
  j.at("g").get_to(d.g);
  j.at("T").get_to(d.T);
  j.at("phi").get_to(d.phi);
  //j.at("F_ext_px").get_to(d.F_ext_px);
//  j.at("F_ext_py").get_to(d.F_ext_py);
//  j.at("F_int_px").get_to(d.F_int_px);
//  j.at("F_int_py").get_to(d.F_int_py);
}


using idx_t = quadgrid_t<std::vector<double>>::idx_t;

int main ()
{
  nlohmann::json json;
  std::ifstream json_file("DATA.json");
  json_file>>json;
  json_file.close();
  DATA data = json.get<DATA>();

  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (data.Ney, data.Nex, data.hx, data.hy);

  idx_t num_particles = data.x.size();
  particles_t ptcls (num_particles, {"label"}, {"Mp", "Ap","vpx","vpy","apx","apy","mom_px","mom_py","hp","F_ext_px","F_ext_py","F_int_px","F_int_py","F_11","F_12","F_21","F_22"}, grid, data.x, data.y);
  ptcls.dprops["Mp"] = data.Mp;
  ptcls.dprops["Ap"] = data.Ap;
  ptcls.dprops["vpx"] = data.vpx; // da matlab
  ptcls.dprops["vpy"] = data.vpy;
  ptcls.dprops["mom_px"] = data.mom_px; // da matlab
  ptcls.dprops["mom_py"] = data.mom_py;
  ptcls.dprops["hp"] = data.hp;
  ptcls.dprops["F_ext_px"].assign (num_particles,0.0);
  ptcls.dprops["F_ext_py"].assign (num_particles,0.0);
  ptcls.dprops["F_int_px"].assign (num_particles,0.0);
  ptcls.dprops["F_int_py"].assign (num_particles,0.0);
  ptcls.dprops["apx"].assign (num_particles,0.0);
  ptcls.dprops["apy"].assign (num_particles,0.0);

  std::iota(ptcls.iprops["label"].begin(),ptcls.iprops["label"].end(),0);


  std::vector<double> Fux(num_particles,0.);
  std::vector<double> Fuy(num_particles,0.);
  std::vector<double> Fvx(num_particles,0.);
  std::vector<double> Fvy(num_particles,0.);

  std::vector<double> sig_xx(num_particles,0.);
  std::vector<double> sig_xy(num_particles,0.);
  std::vector<double> sig_yy(num_particles,0.);

  std::vector<double> Fbx(num_particles,0.);
  std::vector<double> Fby(num_particles,0.);

  std::vector<double> Ftot_vx (grid.num_global_nodes (), 0.0);
  std::vector<double> Ftot_vy (grid.num_global_nodes (), 0.0);

  std::vector<double> avx (grid.num_global_nodes (), 0.0);
  std::vector<double> avy (grid.num_global_nodes (), 0.0);
  std::vector<double> vvx (grid.num_global_nodes (), 0.0);
  std::vector<double> vvy (grid.num_global_nodes (), 0.0);

double t = 0;
double dt = 0.01;



std::map<std::string, std::vector<double>>
  vars{{"Mv", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"mom_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"mom_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"F_ext_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"F_ext_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"F_int_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"F_int_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"avx", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"avy", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"vvx", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"vvy", std::vector<double>(grid.num_global_nodes (), 0.0)}};



//while (t < T)
{


//  (0)  CONNECTIVITY and BASIS FUNCTIONS

/*
p2e  = reshape((permute((floor((Xp(:,2)-min(y))./(hl(2)))+1)+(Ney).*floor((Xp(:,1)-min(x))./(hl(2))),[2 1])),nmp,1);
neon = length(unique(p2e));
p2v  = reshape((e2v(p2e,:)'),4,nmp)';
[S,Sx,Sy,B] = Shptest(Xp(:,1),Xp(:,2),x,y,hl(1),nmp,p2v);
*/

 ptcls.init_particle_mesh ();


/*
% Reset grid quantities
Mv = zeros(NV,1);
MvR = zeros(2*NV,1);
F  = zeros(2*NV,1);
Flux = zeros(2*NV,1);
av = zeros(2*NV,1);
vv = zeros(2*NV,1);
mom = zeros(2*NV,1);


*/

for (auto &v : vars)
{
  v.second.assign (v.second.size (),0.0);
}

// (1) PROJECTION FROM MP TO NODES (P2G)

/*
Projection to nodes
for m = 1:nmp
    for nd = 1:4
        l  = p2v(m,nd); iDy = 2*l;  iDx = iDy-1;

         Mv(l) = Mv(l) + S(m,nd).*Mp(m);

         mom(iDx) = mom(iDx) + S(m,nd).*(momp(m,1));
         mom(iDy) = mom(iDy) + S(m,nd).*(momp(m,2));
  end
end
*/
ptcls.p2g (vars,std::vector<std::string>{"Mp","mom_px","mom_py"},
           std::vector<std::string>{"Mv","mom_vx","mom_vy"});


// (2) INTERNAL AND EXTERNAL FORCES ON VERTICES
/*
for m = 1:nmp
    for nd = 1:4
        l  = p2v(m,nd); iDy = 2*l;  iDx = iDy-1;
        F_ext(iDx) = F_ext(iDx)  -S(m,nd).*rho(m)*Ap(m)*hp(m)*g*(phi) -
                 Ap(m)*S(m,nd)*Fb(m,1);
        F_ext(iDy) = F(iDy)   - Ap(m)*S(m,nd)*Fb(m,2);%*cos(atan(1));
        end
    end





*/

for (idx_t ip = 0; ip < num_particles; ++ip) // FRICTION
{
  double norm = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip]);

  Fbx[ip] = - ((1e5 + data.rho * ptcls.dprops["hp"][ip] * data.g) * std::tan(75 * M_PI / 180) / (norm +1e-4) +
            data.rho * data.g * norm / data.xi) * ptcls.dprops["vpx"][ip];

  Fby[ip] = - ((1e5 + data.rho * ptcls.dprops["hp"][ip] * data.g) * std::tan(75 * M_PI / 180) / (norm +1e-4) +
            data.rho * data.g * norm / data.xi) * ptcls.dprops["vpy"][ip];

/*  std::cout <<std::setprecision(8) << " hp = " << ptcls.dprops["hp"][ip] << " "
            << "rho = " << data.rho               << " "
            << "g = " << data.g               << " "
            << "xi = " << data.xi               << " "
            << "norm = " << norm               <<  " "
            << "tan = " << std::tan(75 * M_PI / 180)  << " "
            << " vpx = " << ptcls.dprops["vpx"][ip] << " "
            << " vpy = " << ptcls.dprops["vpy"][ip] << " "
            << " Fbx = " << Fbx[ip] << " "
            << " Fby = " << Fby[ip] << std::endl;
*/
  }

for (idx_t ip = 0; ip < num_particles; ++ip) // Calcolo forze esterne su ogni particella
{
  /*
  Fb(mp,1) = -((1e5+rho(mp)*hp(mp)*g)*tan(75*pi/180)/(norm(vp(mp,:))+1e-4) +...
            rho(mp)*g*norm(vp(mp,:))/xi)*vp(mp,1);
Fb(mp,2) = -((1e5+rho(mp)*hp(mp)*g)*tan(75*pi/180)/(norm(vp(mp,:))+1e-4) +...
            rho(mp)*g*norm(vp(mp,:))/xi)*vp(mp,2);
            */
  ptcls.dprops.at("F_ext_px")[ip] = -data.rho * ptcls.dprops["Ap"][ip] * ptcls.dprops["hp"][ip] * data.g * data.phi - ptcls.dprops["Ap"][ip] * Fbx[ip];
  ptcls.dprops.at("F_ext_py")[ip] = -data.rho * ptcls.dprops["Ap"][ip] * ptcls.dprops["hp"][ip] * data.g * data.phi - ptcls.dprops["Ap"][ip] * Fby[ip];

/*  std::cout <<std::setprecision(8) << " Ap = " << ptcls.dprops["Ap"][ip] << " "
            << " hp = " << ptcls.dprops["hp"][ip] << " "
            << "rho = " << data.rho               << " "
            << "g = " << data.g               << " "
            << "phi = " << data.phi              << " "
            << "F_ext_px = " << ptcls.dprops.at("F_ext_px")[ip] << " "
            << "F_ext_py = " << ptcls.dprops.at("F_ext_py")[ip] << " "
            << " Fbx = " << Fbx[ip] << " "
            << " Fby = " << Fby[ip] << std::endl;
*/
}

ptcls.p2g (vars,std::vector<std::string>{"F_ext_px","F_ext_py"}, // External Forces on vertices
           std::vector<std::string>{"F_ext_vx","F_ext_vy"}); //   apply mass matrix for confronting the two fields

//std::cout << vars.at("F_ext_vx").size ()<< std::endl;





/* (3) FORZE INTERNE
for m = 1:nmp
for nd = 1:4
    l  = p2v(m,nd); iDy = 2*l;  iDx = iDy-1;
F_int(iDx) = F_int(iDx) + Ap(m)*(Sx(m,nd)*F_11(m,1)+Sy(m,nd)*F_12(m,1))+Ap(m)*(Sx(m,nd)*sig_xx(m,1)+Sy(m,nd)*sig_xy(m,1));
F_int(iDy) = F_int(iDy) + Ap(m)*(Sx(m,nd)*F_21(m,1)+Sy(m,nd)*F_22(m,1))+Ap(m)*(Sx(m,nd)*sig_xy(m,1)+Sy(m,nd)*sig_yy(m,1));
end
end
*/

for (idx_t ip = 0; ip < num_particles; ++ip)
{

  ptcls.dprops["F_11"][ip] = ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip]  +
          .5 * data.g * std::cos(std::atan(data.phi)) * (ptcls.dprops["hp"][ip]  * ptcls.dprops["hp"][ip] );
  ptcls.dprops["F_12"][ip] = ptcls.dprops["hp"][ip]  * ptcls.dprops["vpx"][ip] * ptcls.dprops["vpy"][ip];
  ptcls.dprops["F_21"][ip] = ptcls.dprops["hp"][ip]  * ptcls.dprops["vpx"][ip] * ptcls.dprops["vpy"][ip];
  ptcls.dprops["F_22"][ip] = ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip]  +
          .5 * data.g * std::cos(std::atan(data.phi)) * (ptcls.dprops["hp"][ip]  * ptcls.dprops["hp"][ip] );
}

ptcls.p2gd (vars, std::vector<std::string>{"F_11","F_21"}, std::vector<std::string>{"F_12","F_22"},
            "Ap",std::vector<std::string>{"F_int_vx","F_int_vy"});

// (3) MOMENTUM BALANCE
/*
mom = mom + dt.*F;
*/
for (idx_t iv = 0; iv < grid.num_global_nodes (); ++iv) // SOMMA FORZE INTERNE ED ESTERNE
{
  Ftot_vx[iv] = vars["F_ext_vx"][iv] + vars["F_int_vx"][iv];
  Ftot_vy[iv] = vars["F_ext_vy"][iv] + vars["F_int_vy"][iv];
}

for (idx_t iv = 0; iv < grid.num_global_nodes (); ++iv)   // UPDATE MOMENTI ---> mom = mom + dt*F
{
  vars["mom_vx"][iv] += dt * Ftot_vx[iv];
  vars["mom_vy"][iv] += dt * Ftot_vy[iv];
}

// (4)  COMPUTE NODAL ACCELERATIONS AND VELOCITIES
/*
MvR = reshape(repmat(Mv',2,1),2*NV,1);

av  = F./(MvR);  av(isnan(av) | MvR <=1e-7) = 0.0;
vv  = mom./(MvR);  vv(isnan(vv) | MvR<=1e-7) = 0.0;
*/
for (idx_t iv = 0; iv < grid.num_global_nodes (); ++iv)
{
  avx[iv] = Ftot_vx[iv] / vars["Mv"][iv];
  avy[iv] = Ftot_vy[iv] / vars["Mv"][iv];

  vvx[iv] = vars["mom_vx"][iv] / vars["Mv"][iv];
  vvy[iv] = vars["mom_vy"][iv] / vars["Mv"][iv];
}
// (5) BOUNDARY CONDITIONS
/*
av(bcxl(:,1)) = bcxl(:,2);
av(bcxr(:,1)) = bcxr(:,2);
av(bcyu(:,1)) = bcyu(:,2);
av(bcyd(:,1)) = bcyd(:,2);

vv(bcxl(:,1)) = bcxl(:,2);
vv(bcxr(:,1)) = bcxr(:,2);
vv(bcyu(:,1)) = bcyu(:,2);
vv(bcyd(:,1)) = bcyd(:,2);
*/

// (6) RETURN TO POINTS (G2P)
/*
for mp = 1:nmp
    iDx = 2*p2v(mp,:)-1;  iDy = iDx+1;

    Xp(mp,:) = Xp(mp,:) + dt*[S(mp,:)*vv(iDx) S(mp,:)*vv(iDy)];
    vp(mp,:) = vp(mp,:) + dt*[S(mp,:)*av(iDx) S(mp,:)*av(iDy)];
end
*/
ptcls.g2p (vars,std::vector<std::string>{"vvx","vvy"}, // External Forces on vertices
           std::vector<std::string>{"vpx","vpy"});
ptcls.g2p (vars,std::vector<std::string>{"avx","avy"}, // External Forces on vertices
           std::vector<std::string>{"apx","apy"});

for (idx_t ip = 0; ip < num_particles; ++ip)
{
  data.x[ip] += dt * ptcls.dprops["vpx"][ip];
  data.y[ip] += dt * ptcls.dprops["vpy"][ip];

  ptcls.dprops["vpx"][ip] += dt * ptcls.dprops["apx"][ip];
  ptcls.dprops["vpy"][ip] += dt * ptcls.dprops["apy"][ip];
}


// (7) COMPUTE HEIGHT WITH STRAIN
/*
for mp = 1:nmp
iDx = 2*p2v(mp,:)-1;  iDy = iDx+1;
    momp(mp,:) = vp(mp,:).*repmat(Mp(mp),1,2);
    hp(mp,:) = hp(mp,:)./(1+(Sx(mp,:)*vv(iDx)+Sy(mp,:)*vv(iDy)).*dt);
end
*/

ptcls.p2gd (vars, std::vector<std::string>{"vpx"}, std::vector<std::string>{"vpy"},
            "Ap",std::vector<std::string>{"div_v"});
for (idx_t ip = 0; ip < num_particles; ++ip)
{
  ptcls.dprops["hp"][ip] /= dt * ptcls.dprops["apx"];
}

// (8) UPDATE PARTICLE AREA and STRESS (USL)
/*
for mp = 1:nmp

      Ap(mp) = Mp(mp)/(rho(mp)*(hp(mp)));

     Uh(mp,1) = vp(mp,1).*hp(mp,1);
     Uh(mp,2) = vp(mp,2).*hp(mp,1);

     F_11(mp,1) = vp(mp,1).^2.*hp(mp,1)+0.5*g*cos(phi).*hp(mp,1).^2;
     F_12(mp,1) = vp(mp,1)*vp(mp,2).*hp(mp,1);
     F_21(mp,1) = vp(mp,1)*vp(mp,2).*hp(mp,1);
     F_22(mp,1) = vp(mp,2).^2.*hp(mp,1)+0.5*g.*cos(phi)*hp(mp,1).^2;


    Fb(mp,1) = -((1e5+rho(mp)*hp(mp)*g)*tan(75*pi/180)/(norm(vp(mp,:))+1e-4) +...
                rho(mp)*g*norm(vp(mp,:))/xi)*vp(mp,1);
    Fb(mp,2) = -((1e5+rho(mp)*hp(mp)*g)*tan(75*pi/180)/(norm(vp(mp,:))+1e-4) +...
                rho(mp)*g*norm(vp(mp,:))/xi)*vp(mp,2);



    alpha(mp) = (6*vis*norm(vp(mp,:)))./(hp(mp)*ty);
    b(mp) = -(114/32+alpha(mp));


    z(mp) = (-b(mp)-sqrt(b(mp).^2-4*a*C))/(2*a);

    s_xx(mp) = Sx(mp,:)*vv(iDx);
    s_xy(mp) = 0.5*(Sy(mp,:)*vv(iDx)+Sx(mp,:)*vv(iDy));
    s_yy(mp) = Sy(mp,:)*vv(iDy);


    D_bar(1,1,mp) = s_xx(mp);  D_bar(1,2,mp) = s_xy(mp); D_bar(1,3,mp) = 0.5*(3/(2+z(mp)))*(vp(mp,1)/(hp(mp)+1e-4));
    D_bar(2,1,mp) = s_xy(mp);  D_bar(2,2,mp) = s_yy(mp); D_bar(2,3,mp) = 0.5*(3/(2+z(mp)))*(vp(mp,2)/(hp(mp)+1e-4));
    D_bar(3,1,mp) = D_bar(1,3,mp); D_bar(3,2,mp) = D_bar(2,3,mp); D_bar(3,3,mp) = -(D_bar(1,1,mp)+D_bar(2,2,mp));

    invII(mp) = 0.5*sum(dot(D_bar(:,:,mp),D_bar(:,:,mp)));

    sig_xx(mp) = (ty/sqrt(invII(mp))+2*vis).*s_xx(mp)*hp(mp);
    sig_xy(mp) = (ty/sqrt(invII(mp))+2*vis).*s_xy(mp)*hp(mp);
    sig_yy(mp) = (ty/sqrt(invII(mp))+2*vis).*s_yy(mp)*hp(mp);

end
*/

  ptcls.print<particles_t::output_format::csv>(std::cout);


  grid.vtk_export("GRID.vts", vars);

  t += 5; //dt;
}


/*  for (auto ii : vars["Mv"])
  std::cout << ii << std::endl;




/*
  // (0)
  while (t < T)
  {
    std::cout << "(0)" << std::endl;
    std::map<std::string, std::vector<double>>
      vars{{"m", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"vy", std::vector<double>(grid.num_global_nodes (), 0.0)}};

    // (1)
    ptcls.p2g (vars);

    for (auto ii : vars["m"])
      std::cout << ii << std::endl;


    // (2)
    // (3)
    // (4)
    // (5)
    // (6)
    // (7)
    // (8)
    // (9)
    // (10)
    // (11)
    // (12)
    // (13)
    // (14)
    // (15)
    t += 5;
  }

*/


  return 0;
}
