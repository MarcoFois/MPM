
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
  double g;
  double T;
  double rho;
  std::vector<double> Vp;
  std::vector<double> bound;
  double xi;
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
  j.at("g").get_to(d.g);
  j.at("T").get_to(d.T);
  j.at("Vp").get_to(d.Vp);
  j.at("bound").get_to(d.bound);
  j.at("xi").get_to(d.xi);

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
    particles_t ptcls (num_particles, {"label"}, {"Mp", "Ap","vpx","vpy","mom_px","mom_py","hp",
                                                  "Vp","F_ext_px","F_ext_py","apx","apy","VpR",
                                                  "F_11","F_12","F_21","F_22","vpx_dx","vpx_dy",
                                                  "vpy_dx","vpy_dy","Fb_x","Fb_y"}, grid, data.x, data.y);
    ptcls.dprops["Mp"] = data.Mp;
    ptcls.dprops["Ap"] = data.Ap;
    ptcls.dprops["vpx"] = data.vpx;
    ptcls.dprops["vpy"] = data.vpy;
    ptcls.dprops["mom_px"] = data.mom_px;
    ptcls.dprops["mom_py"] = data.mom_py;
    ptcls.dprops["hp"] = data.hp;
    ptcls.dprops["Vp"] = data.Vp;


    for (idx_t ip = 0; ip < num_particles; ++ip)
    {
        ptcls.dprops["vpx_dx"][ip] = 0.0;
        ptcls.dprops["vpx_dy"][ip] = 0.0;
        ptcls.dprops["vpy_dx"][ip] = 0.0;
        ptcls.dprops["vpy_dy"][ip] = 0.0;
        ptcls.dprops["apx"][ip] = 0.0;
        ptcls.dprops["apy"][ip] = 0.0;
        ptcls.dprops["VpR"][ip] = ptcls.dprops["Vp"][ip] *  data.rho ;
        ptcls.dprops["Fb_x"][ip] = 0.0;
        ptcls.dprops["Fb_y"][ip] = 0.0;
    }

    std::iota(ptcls.iprops["label"].begin(),ptcls.iprops["label"].end(),0);

    double t = 0.0;
    double dt  ;
    double cel;
    double phi = 1./1.732050;
    double atm = 100000.;
    double fric_ang = 25. * M_PI / 180.;

    for (idx_t ip = 0; ip < num_particles; ++ip)
    {
      ptcls.dprops["F_11"][ip] =   data.g  * std::cos(std::atan(phi)) *  (ptcls.dprops["hp"][ip]  ) ;
      ptcls.dprops["F_12"][ip] = 0.0;
      ptcls.dprops["F_21"][ip] = 0.0;
      ptcls.dprops["F_22"][ip] =  data.g  * std::cos(std::atan(phi)) *  (ptcls.dprops["hp"][ip]  );
    }

    std::vector<double> Ftot_vx (grid.num_global_nodes (), 0.0);
    std::vector<double> Ftot_vy (grid.num_global_nodes (), 0.0);
    std::vector<double> norm_v (num_particles, 0.0);








    std::map<std::string, std::vector<double>>
      vars{{"Mv", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"mom_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"mom_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_ext_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_ext_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_int_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_int_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"div_v", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"avx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"avy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"vvx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"vvy", std::vector<double>(grid.num_global_nodes (), 0.0)}},
      Plotvars{{"rho_v", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"avx", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"avy", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"vvx", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"vvy", std::vector<double>(grid.num_global_nodes (), 0.0)}};

    int it = 0;
    ptcls.build_mass();

    while (t < 1)
    {

        double max_vel_x = *std::max_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());
        double max_vel_y = *std::max_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());
        double max_vel = std::max(1+max_vel_x,1+ max_vel_y);
        cel = std::abs( max_vel);
        dt = 0.01 * data.hx / cel;

        std::string filename = "nc_particles_";
        filename = filename + std::to_string (it++);
        filename = filename + ".csv";
        std::ofstream OF (filename.c_str ());
        ptcls.print<particles_t::output_format::csv>(OF);
        OF.close ();



      //  (0)  CONNECTIVITY and BASIS FUNCTIONS

        ptcls.init_particle_mesh ();


        for (auto &v : vars)
        {
            v.second.assign (v.second.size (),0.0);
        }

        for (auto &v : Plotvars)
        {
            v.second.assign (v.second.size (),0.0);
        }




      // (1) PROJECTION FROM MP TO NODES (P2G)

        ptcls.p2g (vars,std::vector<std::string>{"Mp","mom_px","mom_py"},
                 std::vector<std::string>{"Mv","mom_vx","mom_vy"});


      // (2)  EXTERNAL FORCES ON VERTICES (P2G)


        for (idx_t ip = 0; ip < num_particles; ++ip)
        {
            ptcls.dprops["F_ext_px"][ip] += ptcls.dprops["Mp"][ip] * data.g * phi + ptcls.dprops["Ap"][ip] * ptcls.dprops["Fb_x"][ip];
            ptcls.dprops["F_ext_py"][ip] +=  ptcls.dprops["Ap"][ip] * ptcls.dprops["Fb_y"][ip];
        }



        ptcls.p2g (vars,std::vector<std::string>{"F_ext_px","F_ext_py"},
                 std::vector<std::string>{"F_ext_vx","F_ext_vy"});


      // (3) INTERNAL FORCES (p2gd) and MOMENTUM BALANCE



        ptcls.p2gd (vars, std::vector<std::string>{"F_11","F_21"},
                  std::vector<std::string>{"F_12","F_22"},
                  "VpR",std::vector<std::string>{"F_int_vx","F_int_vy"});



        for (idx_t iv = 0; iv < grid.num_global_nodes (); ++iv)
        {
            Ftot_vx[iv] = vars["F_ext_vx"][iv] + vars["F_int_vx"][iv];
            Ftot_vy[iv] = vars["F_ext_vy"][iv] + vars["F_int_vy"][iv];
        }

        for (idx_t iv = 0; iv < grid.num_global_nodes (); ++iv)
        {
            vars["mom_vx"][iv] += dt * Ftot_vx[iv];
            vars["mom_vy"][iv] += dt * Ftot_vy[iv];
        }

      // (4)  COMPUTE NODAL ACCELERATIONS AND VELOCITIES

        for (idx_t iv = 0; iv < grid.num_global_nodes (); ++iv)
        {
            vars["avx"][iv] = vars["Mv"][iv] > 1e-4 ?  Ftot_vx[iv] / vars["Mv"][iv] : 0.0;
            vars["avy"][iv] = vars["Mv"][iv] > 1e-4 ?  Ftot_vy[iv] / vars["Mv"][iv] : 0.0;

            vars["vvx"][iv] = vars["Mv"][iv] > 1e-4 ?  vars["mom_vx"][iv] / vars["Mv"][iv] : 0.0;
            vars["vvy"][iv]= vars["Mv"][iv] > 1e-4 ?  vars["mom_vy"][iv] / vars["Mv"][iv] : 0.0;
        }

      // (5) BOUNDARY CONDITIONS - TO DO
       for (auto icell = grid.begin_cell_sweep ();
	      icell != grid.end_cell_sweep (); ++icell)
        {
           if ( (icell->e (2) == 2) || (icell->e(1)==3)  )
           {
              for (idx_t inode = 0; inode < 4; ++inode)
              {
                 vars["avx"][icell->gt(inode)] = 0.0;
                 vars["vvx"][icell->gt(inode)] = 0.0;
	            }
           }
	         if ( (icell->e(0) == 0) || (icell->e(1)==1) )
  	       {
  	          for (idx_t inode = 0; inode < 4; ++inode)
  	          {
                 vars["avy"][icell->gt(inode)] = 0.0;
                 vars["vvy"][icell->gt(inode)] = 0.0;
  	          }
	         }
        }

      /*  for (auto iv : data.bound)
        {
          vars["avx"][iv] = 0.0;
          vars["vvx"][iv] = 0.0;
          vars["avy"][iv] = 0.0;
          vars["vvy"][iv] = 0.0;
        } */




      // (6) RETURN TO POINTS (G2P) and UPDATE POS AND VEL ON PARTICLES

        ptcls.dprops.at("vpx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("vpy").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apy").assign(ptcls.num_particles, 0.0);

        ptcls.g2p (vars,std::vector<std::string>{"vvx","vvy","avx","avy"},
                 std::vector<std::string>{"vpx","vpy","apx","apy"});


        for (idx_t ip = 0; ip < num_particles; ++ip)
        {
            ptcls.dprops["vpx"][ip] += dt * ptcls.dprops["apx"][ip];
            ptcls.dprops["vpy"][ip] += dt * ptcls.dprops["apy"][ip];
        }

        for (idx_t ip = 0; ip < num_particles; ++ip)
        { //ptcls.x[ip]
            ptcls.x[ip] += dt * ptcls.dprops["vpx"][ip];
            ptcls.y[ip] += dt * ptcls.dprops["vpy"][ip];
        }



      // (7) COMPUTE HEIGHT WITH STRAIN (divergence of velocities)

        ptcls.g2pd (vars,std::vector<std::string>{"vvx","vvy"},
                 std::vector<std::string>{"vpx_dx","vpy_dx"},
                 std::vector<std::string>{"vpx_dy","vpy_dy"});

        for (idx_t ip = 0; ip < num_particles; ++ip)
        {
            ptcls.dprops["hp"][ip] /= (1+dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip]));

            ptcls.dprops["mom_px"][ip] = ptcls.dprops["vpx"][ip] * ptcls.dprops["Mp"][ip];
            ptcls.dprops["mom_py"][ip] = ptcls.dprops["vpy"][ip] * ptcls.dprops["Mp"][ip];

            ptcls.dprops["Vp"][ip]  /=(1+dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip]));
            ptcls.dprops["VpR"][ip] = ptcls.dprops["Vp"][ip]  * data.rho ;
            ptcls.dprops["Ap"][ip] = ptcls.dprops["Vp"][ip] / (ptcls.dprops["hp"][ip] + 0.0001);
        }

      // (8) UPDATE PARTICLE STRESS (USL)

        for (idx_t ip = 0; ip < num_particles; ++ip)
        {
            ptcls.dprops["F_11"][ip] =   data.g * std::cos(std::atan(phi)) *  (ptcls.dprops["hp"][ip]   ) ;
            ptcls.dprops["F_12"][ip] = 0.0;
            ptcls.dprops["F_21"][ip] = 0.0;
            ptcls.dprops["F_22"][ip] =  data.g * std::cos(std::atan(phi)) *  (ptcls.dprops["hp"][ip]  );
        }

        for (idx_t ip = 0; ip < num_particles; ++ip)
        {
            norm_v[ip] = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip] );

            ptcls.dprops["Fb_x"][ip] =  ((atm + data.rho * ptcls.dprops["hp"][ip] * data.g) * std::tan(fric_ang) / (norm_v[ip] + 0.0001) +
                                        data.rho * data.g * norm_v[ip] / data.xi) * ptcls.dprops["vpx"][ip]  ;

            ptcls.dprops["Fb_y"][ip] =  ((atm + data.rho * ptcls.dprops["hp"][ip] * data.g) * std::tan(fric_ang) /( norm_v[ip] + 0.0001) +
                                        data.rho * data.g * norm_v[ip] / data.xi) * ptcls.dprops["vpy"][ip] ;
        }

        ptcls.p2g (Plotvars,std::vector<std::string>{"Mp","vpx","vpy","apx","apy"},
                std::vector<std::string>{"rho_v","vvx","vvy","avx","avy"},true);



        filename = "nc_grid_";
        filename = filename + std::to_string (it);
        filename = filename + ".vts";
        grid.vtk_export(filename.c_str(), Plotvars);

        t +=dt;


    }

//  ptcls.print<particles_t::output_format::csv>(std::cout);

  return 0;
}
