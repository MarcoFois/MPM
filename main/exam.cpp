
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
  std::vector<double> meanZ;
  std::vector<double> Z;
  std::vector<double> dZdx;
  std::vector<double> dZdy;
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
  j.at("meanZ").get_to(d.meanZ);
  j.at("Z").get_to(d.Z);
  j.at("dZdx").get_to(d.dZdx);
  j.at("dZdy").get_to(d.dZdy);

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
                                                  "Vp","F_ext_px","F_ext_py","apx","apy",
                                                  "F_11","F_12","F_21","F_22","vpx_dx","vpx_dy",
                                                  "vpy_dx","vpy_dy","Fb_x","Fb_y","hpZ","dZxp","dZyp","Zp","xp","yp","Fric_px","Fric_py"}, grid, data.x, data.y);
    ptcls.dprops["Mp"] = data.Mp;
    ptcls.dprops["Ap"] = data.Ap;
    ptcls.dprops["vpx"] = data.vpx;
    ptcls.dprops["vpy"] = data.vpy;
    ptcls.dprops["mom_px"] = data.mom_px;
    ptcls.dprops["mom_py"] = data.mom_py;
    ptcls.dprops["hp"] = data.hp;
    ptcls.dprops["Vp"] = data.Vp;
    ptcls.dprops["xp"] = data.x;
    ptcls.dprops["yp"] = data.y;



    for (idx_t ip = 0; ip < num_particles; ++ip)
    {
        ptcls.dprops["vpx_dx"][ip] = 0.0;
        ptcls.dprops["vpx_dy"][ip] = 0.0;
        ptcls.dprops["vpy_dx"][ip] = 0.0;
        ptcls.dprops["vpy_dy"][ip] = 0.0;
        ptcls.dprops["apx"][ip] = 0.0;
        ptcls.dprops["apy"][ip] = 0.0;
        ptcls.dprops["Fb_x"][ip] = 0.0;
        ptcls.dprops["Fb_y"][ip] = 0.0;
        ptcls.dprops["dZxp"][ip] = 0.0;
        ptcls.dprops["dZyp"][ip] = 0.0;
        ptcls.dprops["Zp"][ip] = 0.0;
        ptcls.dprops["hpZ"][ip] = 0.0;
    }


    std::iota(ptcls.iprops["label"].begin(),ptcls.iprops["label"].end(),0);

    double t = 0.0;
    double dt  ;
    double cel;
    double phi = 0.0;
    double atm = 100000.;
    double fric_ang = 1. * M_PI / 40.;
    double atan_grad_z;
    std::vector<double> norm_v (num_particles, 0.0);




    std::map<std::string, std::vector<double>>
      vars{{"Mv", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"mom_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"mom_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_ext_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_ext_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_int_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"F_int_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"Fric_x", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"Fric_y", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"div_v", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"avx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"avy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"vvx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"vvy", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"Z", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"dZdx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"dZdy", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"Ftot_vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
    {"Ftot_vy", std::vector<double>(grid.num_global_nodes (), 0.0)}},
      Plotvars{{"rho_v", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"avx", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"avy", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"vvx", std::vector<double>(grid.num_global_nodes (), 0.0)},
          {"vvy", std::vector<double>(grid.num_global_nodes (), 0.0)}};


         vars["Z"] = data.Z;
         vars["dZdx"] = data.dZdx;
         vars["dZdy"] = data.dZdy;

    /*      for (auto icell = grid.begin_cell_sweep ();
                   icell != grid.end_cell_sweep (); ++icell)
               {
                   for (auto iv = 0; iv < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++iv)
                   {
                     auto gv = icell -> gt (iv);

              //     vars["dZdx"][gv] = - 0.30 - 2. * (icell -> p(0,gv) -50.) * std::exp(-(icell -> p(0,gv)-50.) * (icell -> p(0,gv)-50.)/7. - (icell -> p(1,gv)-10.) * (icell -> p(1,gv)-10.)/7.);//2. * std::cos(icell->p (0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
                //   vars["dZdy"][gv] = - 2. * (icell -> p(1,gv)-10.) * std::exp(-(icell -> p(0,gv)-50.) * (icell -> p(0,gv)-50.)/7. - (icell -> p(1,gv)-10.) * (icell -> p(1,gv)-10.)/7.);            // -2. * 0.5 * std::sin(icell -> p(0,gv)) * std::sin(0.5 * (icell-> p(1,gv)));
                  //   vars["Z"][gv] = 50. - 0.30 * icell -> p(0,gv) + 7. * std::exp(-(icell -> p(0,gv)-50.) * (icell -> p(0,gv)-50.)/7. - (icell -> p(1,gv)-10.) * (icell -> p(1,gv)-10.)/7.); //3. + 2. * std::sin(icell -> p(0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
                  vars["dZdx"][gv] = 2. * std::cos(icell->p (0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
                  vars["dZdy"][gv] =  -2. * 0.5 * std::sin(icell -> p(0,gv)) * std::sin(0.5 * (icell-> p(1,gv)));
                    vars["Z"][gv] = 3. + 2. * std::sin(icell -> p(0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
                   }
               } */

          ptcls.g2p (vars,std::vector<std::string>{"Z"},
                               std::vector<std::string>{"Zp"});
          ptcls.dprops.at("dZxp").assign(ptcls.num_particles, 0.0);
          ptcls.dprops.at("dZyp").assign(ptcls.num_particles, 0.0);
          ptcls.g2p (vars,std::vector<std::string>{"dZdx","dZdy"},
                    std::vector<std::string>{"dZxp","dZyp"});

          for (idx_t ip = 0; ip < num_particles; ++ip)
          {
            ptcls.dprops["hpZ"][ip] = ptcls.dprops["hp"][ip] + ptcls.dprops["Zp"][ip];
          }

    int it = 0;
    ptcls.build_mass();

   grid.vtk_export("GRID_forZ.vts", vars);


    while (t < 4)
    {

        double max_vel_x = *std::max_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());
        double max_vel_y = *std::max_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());
        double max_vel = std::max(1+max_vel_x,1+max_vel_y);
        cel = std::abs( max_vel);
        dt = 0.01 *  data.hx / (1e-4 + cel); //0.2 *  data.hx / (1e-4 + cel);

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

        for (idx_t ip = 0; ip<num_particles; ++ip)
        {
          ptcls.dprops["F_11"][ip] =   .5 * data.rho * data.g *   (ptcls.dprops["hp"][ip]   ) ;
          ptcls.dprops["F_12"][ip] = 0.0;
          ptcls.dprops["F_21"][ip] = 0.0;
          ptcls.dprops["F_22"][ip] =  .5 * data.rho * data.g *   (ptcls.dprops["hp"][ip] );
        }


        ptcls.dprops.at("dZxp").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("dZyp").assign(ptcls.num_particles, 0.0);
        vars["Z"] = data.Z;
        vars["dZdx"] = data.dZdx;
        vars["dZdy"] = data.dZdy;

  /*      for (auto icell = grid.begin_cell_sweep ();
                 icell != grid.end_cell_sweep (); ++icell)
             {
                 for (auto iv = 0; iv < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++iv)
                 {
                   auto gv = icell -> gt (iv);

            //     vars["dZdx"][gv] = - 0.30 - 2. * 17./7. * (icell -> p(0,gv) -50.) * std::exp(-(icell -> p(0,gv)-50.) * (icell -> p(0,gv)-50.)/7. - (icell -> p(1,gv)-10.) * (icell -> p(1,gv)-10.)/7.);//2. * std::cos(icell->p (0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
            //     vars["dZdy"][gv] = - 2. * 17./7. * (icell -> p(1,gv)-10.) * std::exp(-(icell -> p(0,gv)-50.) * (icell -> p(0,gv)-50.)/7. - (icell -> p(1,gv)-10.) * (icell -> p(1,gv)-10.)/7.);            // -2. * 0.5 * std::sin(icell -> p(0,gv)) * std::sin(0.5 * (icell-> p(1,gv)));
            //       vars["Z"][gv] = 50. - 0.30 * icell -> p(0,gv) + 17. * std::exp(-(icell -> p(0,gv)-50.) * (icell -> p(0,gv)-50.)/7. - (icell -> p(1,gv)-10.) * (icell -> p(1,gv)-10.)/7.); //3. + 2. * std::sin(icell -> p(0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
            vars["dZdx"][gv] = 2. * std::cos(icell->p (0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
            vars["dZdy"][gv] =  -2. * 0.5 * std::sin(icell -> p(0,gv)) * std::sin(0.5 * (icell-> p(1,gv)));
              vars["Z"][gv] = 3. + 2. * std::sin(icell -> p(0,gv)) * std::cos(0.5 * (icell -> p(1,gv)));
                 }
             } */

      // (1) PROJECTION FROM MP TO NODES (P2G)

        ptcls.p2g (vars,std::vector<std::string>{"Mp","mom_px","mom_py"},
                 std::vector<std::string>{"Mv","mom_vx","mom_vy"});
                 ptcls.g2pd (vars,std::vector<std::string>{"Z"},
                          std::vector<std::string>{"dZxp"},
                          std::vector<std::string>{"dZyp"});
      // (2)  EXTERNAL FORCES ON VERTICES (P2G)
     for (auto icell = grid.begin_cell_sweep ();
           icell != grid.end_cell_sweep (); ++icell)
       {
         for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
         {
           auto gv = icell -> gt (inode);

            if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
             for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip)
             {
               auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];


                 ptcls.dprops["Fric_px"][gp] +=  -  ptcls.dprops["Ap"][gp] * ptcls.dprops["Fb_x"][gp];

                 ptcls.dprops["Fric_py"][gp] +=    -   ptcls.dprops["Ap"][gp] *  ptcls.dprops["Fb_y"][gp];
            }
         }

       }

       ptcls.p2g (vars,std::vector<std::string>{"Fric_px","Fric_py"},
                  std::vector<std::string>{"Fric_x","Fric_y"});


      for (auto icell = grid.begin_cell_sweep ();
             icell != grid.end_cell_sweep (); ++icell)
         {
           for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
           {
             auto gv = icell -> gt (inode);

                   vars["F_ext_vx"][gv] -=  vars["Mv"][gv] *  9.81 * vars["dZdx"][gv] - vars["Fric_x"][gv]; // vars["dZdx"][gv]; // vars["dZdx"][gv]; //  dZdx[gv]  ; //-  ptcls.dprops["Ap"][gp] * ptcls.dprops["Fb_x"][gp];

                   vars["F_ext_vy"][gv] -=  vars["Mv"][gv] *  9.81 *  vars["dZdy"][gv] - vars["Fric_y"][gv]; // vars["dZdy"][gv]; //  -   ptcls.dprops["Ap"][gp] *  ptcls.dprops["Fb_y"][gp];

           }

         }


    //   ptcls.p2g (vars,std::vector<std::string>{"F_ext_px","F_ext_py"},
          //        std::vector<std::string>{"F_ext_vx","F_ext_vy"});



      // (3) INTERNAL FORCES (p2gd) and MOMENTUM BALANCE


        ptcls.p2gd (vars, std::vector<std::string>{"F_11","F_21"},
                  std::vector<std::string>{"F_12","F_22"},
                  "Vp",std::vector<std::string>{"F_int_vx","F_int_vy"});



       for (auto icell = grid.begin_cell_sweep ();
            icell != grid.end_cell_sweep (); ++icell)
        {
            for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                {
                  auto iv = icell -> gt (inode);
                  vars["Ftot_vx"][iv] = vars["F_ext_vx"][iv] + vars["F_int_vx"][iv];
                  vars["Ftot_vy"][iv] = vars["F_ext_vy"][iv] + vars["F_int_vy"][iv];
                }

        }

        for (auto icell = grid.begin_cell_sweep ();
             icell != grid.end_cell_sweep (); ++icell)
         {
             for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                 {
                   auto iv = icell -> gt (inode);
                   vars["mom_vx"][iv] += dt * vars["Ftot_vx"][iv];
                   vars["mom_vy"][iv] += dt * vars["Ftot_vy"][iv];
                //   vars["sum"][iv] += 1.;
                 }

         }

      // (4)  COMPUTE NODAL ACCELERATIONS AND VELOCITIES

     for (auto icell = grid.begin_cell_sweep ();
           icell != grid.end_cell_sweep (); ++icell)
       {
           for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
               {
                 auto iv = icell -> gt (inode);
                 vars["avx"][iv] = vars["Mv"][iv] > 1e-7 ?  vars["Ftot_vx"][iv] / vars["Mv"][iv] : 0.0;
                 vars["avy"][iv] = vars["Mv"][iv] > 1e-7 ? vars[ "Ftot_vy"][iv] / vars["Mv"][iv] : 0.0;

                 vars["vvx"][iv] = vars["Mv"][iv] > 1e-7 ?  vars["mom_vx"][iv] / vars["Mv"][iv] : 0.0;
                 vars["vvy"][iv]= vars["Mv"][iv] > 1e-7 ?  vars["mom_vy"][iv] / vars["Mv"][iv] : 0.0;
               }

       }


      // (5) BOUNDARY CONDITIONS - TO DO
       for (auto icell = grid.begin_cell_sweep ();
	      icell != grid.end_cell_sweep (); ++icell)
        {
          /* if ( (icell->e (2) == 2) || (icell->e(1)==3)  )
           {
              for (idx_t inode = 0; inode < 4; ++inode)
              {
                 vars["avx"][icell->gt(inode)] = 0.0;
                 vars["vvx"][icell->gt(inode)] = 0.0;
	            }
           } */
	         if ( (icell->e(0) == 0) || (icell->e(1)==1) )
  	       {
  	          for (idx_t inode = 0; inode < 4; ++inode)
  	          {
                 vars["avy"][icell->gt(inode)] = 0.0;
                 vars["vvy"][icell->gt(inode)] = 0.0;
  	          }
	         }
        }


      // (6) RETURN TO POINTS (G2P) and UPDATE POS AND VEL ON PARTICLES

        ptcls.dprops.at("vpx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("vpy").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apy").assign(ptcls.num_particles, 0.0);

        ptcls.g2p (vars,std::vector<std::string>{"vvx","vvy","avx","avy"},
                 std::vector<std::string>{"vpx","vpy","apx","apy"});


      for (auto icell = grid.begin_cell_sweep ();
                         icell != grid.end_cell_sweep (); ++icell)
            {
               for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                  {
                     if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
                     for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip)
                        {
                         auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];
                          ptcls.dprops["vpx"][gp] += dt * ptcls.dprops["apx"][gp] ;
                          ptcls.dprops["vpy"][gp] += dt * ptcls.dprops["apy"][gp] ;
                        }
                   }

            }



           for (auto icell = grid.begin_cell_sweep ();
                             icell != grid.end_cell_sweep (); ++icell)
                {
                   for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                      {
                         if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
                         for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip)
                            {
                             auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];
                              ptcls.x[gp] += dt * ptcls.dprops["vpx"][gp];
                              ptcls.y[gp] += dt * ptcls.dprops["vpy"][gp];
                            //  ptcls.dprops["xp"][gp] += dt * ptcls.dprops["vpx"][gp];
                          //    ptcls.dprops["yp"][gp] += dt * ptcls.dprops["vpy"][gp];
                            }
                       }

                }


      // (7) COMPUTE HEIGHT WITH STRAIN (divergence of velocities)
      ptcls.dprops.at("vpx_dx").assign(ptcls.num_particles, 0.0);
      ptcls.dprops.at("vpy_dy").assign(ptcls.num_particles, 0.0);
              ptcls.g2pd (vars,std::vector<std::string>{"vvx","vvy"},
                 std::vector<std::string>{"vpx_dx","vpy_dx"},
                 std::vector<std::string>{"vpx_dy","vpy_dy"});

      for (auto icell = grid.begin_cell_sweep ();
            icell != grid.end_cell_sweep (); ++icell)
            {
                for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                   {
                     if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
                     for (auto gp = 0; gp < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++gp)
                        {
                          auto ip = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[gp];
                           ptcls.dprops["hp"][ip] /= (1+dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip]));
                           ptcls.dprops["mom_px"][ip] = ptcls.dprops["vpx"][ip] * ptcls.dprops["Mp"][ip];
                           ptcls.dprops["mom_py"][ip] = ptcls.dprops["vpy"][ip] * ptcls.dprops["Mp"][ip];
                           ptcls.dprops["Vp"][ip]  /=(1+dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip]));

                         }
                   }

            }



      // (8) UPDATE PARTICLE STRESS (USL)
      for (auto icell = grid.begin_cell_sweep ();
          icell != grid.end_cell_sweep (); ++icell)
          {
              for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                 {
                   auto gv = icell -> gt (inode);
                   if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
                   for (auto gp = 0; gp < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++gp)
                      {


                         auto ip = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[gp];

                         ptcls.dprops["F_11"][ip] =  .5 * 1000. * 9.81 *  (ptcls.dprops["hp"][ip] ) ;
                         ptcls.dprops["F_12"][ip] = 0.0;
                         ptcls.dprops["F_21"][ip] = 0.0;
                         ptcls.dprops["F_22"][ip] =  .5 * 1000. * 9.81 *   (ptcls.dprops["hp"][ip]  );
                       }
                 }

          }

          ptcls.dprops.at("hpZ").assign(ptcls.num_particles, 0.0);
          ptcls.dprops.at("Zp").assign(ptcls.num_particles, 0.0);
          ptcls.g2p (vars,std::vector<std::string>{"Z"},
                               std::vector<std::string>{"Zp"});

          for (auto icell = grid.begin_cell_sweep ();
              icell != grid.end_cell_sweep (); ++icell)
                {
                  for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                    {
                      auto gv = icell -> gt (inode);
                      if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
                      for (auto gp = 0; gp < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++gp)
                         {
                        auto ip = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[gp];
                        ptcls.dprops["hpZ"][ip] = ptcls.dprops["hp"][ip] + ptcls.dprops["Zp"][ip];

                     }
                   }

                }
      /*    for (idx_t ip = 0; ip < num_particles; ++ip)
          {
            ptcls.dprops["hpZ"][ip] = ptcls.dprops["hp"][ip] + ptcls.dprops["Zp"][ip];
          } */

     for (idx_t ip = 0; ip < num_particles; ++ip)
        {
            norm_v[ip] = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip] );

            ptcls.dprops["Fb_x"][ip] = - 0.5 * ((atm + data.rho * ptcls.dprops["hp"][ip] * data.g) * std::tan(fric_ang) / (norm_v[ip] + 0.0001) +
                                        data.rho * data.g * norm_v[ip] / data.xi) * ptcls.dprops["vpx"][ip]  ;

            ptcls.dprops["Fb_y"][ip] = - 0.5 * ((atm + data.rho * ptcls.dprops["hp"][ip] * data.g) * std::tan(fric_ang) /( norm_v[ip] + 0.0001) +
                                        data.rho * data.g * norm_v[ip] / data.xi) * ptcls.dprops["vpy"][ip] ;
        }



  /*      for (auto icell = grid.begin_cell_sweep ();
                icell != grid.end_cell_sweep (); ++icell)
            {
              for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
              {

                  if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
                  for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip)
                  {
                    auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];

                      ptcls.dprops["hpZ"][gp] = ptcls.dprops["hp"][gp] + meanZ[icell];

                  }
              }

         } */

//        ptcls.p2g (Plotvars,std::vector<std::string>{"Mp","vpx","vpy","apx","apy"},
        //        std::vector<std::string>{"rho_v","vvx","vvy","avx","avy"});



        filename = "nc_grid_";
        filename = filename + std::to_string (it);
        filename = filename + ".vts";
      //  grid.vtk_export(filename.c_str(), Plotvars);
grid.vtk_export(filename.c_str(), vars);
        t +=dt;


    }

//  ptcls.print<particles_t::output_format::csv>(std::cout);

  return 0;
}
