#include <fstream>
#include <functional>
#include <map>
#include <cmath>
#include <iostream>
#include <particles.h>
#include <quadgrid_cpp.h>
#include <timer.h>
#include "mpm_data.h"

using idx_t = quadgrid_t<std::vector<double>>::idx_t;
cdf::timer::timer_t my_timer{};

int main ()
{

  DATA data ("DATA.json");

  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (data.Ney, data.Nex, data.hx, data.hy);

  idx_t num_particles = data.x.size();
  particles_t ptcls (num_particles, {"label"}, {"Mp", "Ap","vpx","vpy","mom_px","mom_py","hp",
        "Vp","F_ext_px","F_ext_py","apx","apy",
        "F_11","F_12","F_21","F_22","vpx_dx","vpx_dy",
        "vpy_dx","vpy_dy","Fb_x","Fb_y","hpZ","dZxp","dZyp","Zp","xp","yp","Fric_px","Fric_py","Fpx","Fpy"}, grid, data.x, data.y);
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
      ptcls.dprops["Fpx"][ip] = 0.0;
      ptcls.dprops["Fpy"][ip] = 0.0;
    }


  std::iota(ptcls.iprops["label"].begin(),ptcls.iprops["label"].end(),0);

  double t = 0.0;
  double dt  ;
  double cel;
  double phi = 0.0;
  double atm = 100000.;
  double fric_ang = 0.5 * M_PI / 6.;
  double atan_grad_z;
  std::vector<double> norm_v (num_particles, 0.0);

  // double A = 3./2.;
  // std::vector<double> ALF (num_particles, 0.0);
  // std::vector<double> B (num_particles, -114./32.);
  // double C = 65./32.;
  // std::vector<double> s_xx (num_particles, 0.0);
  // std::vector<double> s_xy (num_particles, 0.0);
  // std::vector<double> s_yy (num_particles, 0.0);
  // std::vector<double> D_xx (num_particles, 0.0);
  // std::vector<double> D_xy (num_particles, 0.0);
  // std::vector<double> D_xz (num_particles, 0.0);
  // std::vector<double> D_yx (num_particles, 0.0);
  // std::vector<double> D_yy (num_particles, 0.0);
  // std::vector<double> D_yz (num_particles, 0.0);
  // std::vector<double> D_zx (num_particles, 0.0);
  // std::vector<double> D_zy (num_particles, 0.0);
  // std::vector<double> D_zz (num_particles, 0.0);
  // std::vector<double> invII (num_particles, 0.0);
  // std::vector<double> Z1 (num_particles, 0.0);
  // std::vector<double> Z2 (num_particles, 0.0);
  // std::vector<double> ZZ (num_particles, 0.0);
  // std::vector<double> sig_xx (num_particles, 0.0);
  // std::vector<double> sig_xy (num_particles, 0.0);
  // std::vector<double> sig_yy (num_particles, 0.0);
  double A{3./2.}, B{-114./32.}, C{65./32.};
  double nrm{0.0}, ALF{0.0}, Z1{0.0}, Z2{0.0}, ZZ{0.0};
  double s_xx{0.0}, s_xy{0.0}, s_yy{0.0};
  double D_xx{0.0}, D_yx{0.0}, D_zx{0.0}, D_xy{0.0}, D_yy{0.0}, D_zy{0.0}, D_xz{0.0}, D_yz{0.0}, D_zz{0.0};
  double invII{0.0}, sig_yy{0.0}, sig_xy{0.0}, sig_xx{0.0}, cc{1.0};

  std::map<std::string, std::vector<double>>
    vars{
      {"Mv", std::vector<double>(grid.num_global_nodes (), 0.0)},
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
      {"Ftot_vy", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"HV", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"FPxv", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"FPyv", std::vector<double>(grid.num_global_nodes (), 0.0)}
    },
    Plotvars{
      {"rho_v", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"avx", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"avy", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"vvx", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"vvy", std::vector<double>(grid.num_global_nodes (), 0.0)}
    };


    vars["Z"] = data.Z;
    vars["dZdx"] = data.dZdx;
    vars["dZdy"] = data.dZdy;


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

      for (idx_t ip = 0; ip<num_particles; ++ip)
      {
          ptcls.dprops["F_11"][ip] =   .5 * data.rho * data.g *   (ptcls.dprops["hp"][ip]   ) ;
          ptcls.dprops["F_12"][ip] = 0.0;
          ptcls.dprops["F_21"][ip] = 0.0;
          ptcls.dprops["F_22"][ip] =  .5 * data.rho * data.g *   (ptcls.dprops["hp"][ip] );
      }

    int it = 0;

    ptcls.build_mass();


    grid.vtk_export("GRID_forZ.vts", vars);


    dt = 1.0e-5;
    std::vector<idx_t> ordering (ptcls.num_particles);
    while (t < 2.6) //data.T
      {

        my_timer.tic ("update dt");

	double max_vel_x = *std::max_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());
	double min_vel_x = *std::min_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());
	
	max_vel_x = std::max (std::abs(max_vel_x), std::abs(min_vel_x));
	
        double max_vel_y = *std::max_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());
	double min_vel_y = *std::min_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());
	
	max_vel_y = std::max (std::abs(max_vel_y), std::abs(min_vel_y));
	
        double hmean = *std::max_element (ptcls.dprops["hp"].begin(), ptcls.dprops["hp"].end());
        double max_vel = std::max(std::sqrt(data.g * hmean) + max_vel_x, std::sqrt(data.g * hmean) + max_vel_y);
        //  double max_vel = std::max(1+max_vel_x,1+max_vel_y);
        cel = std::abs(max_vel);
	
        // double max_vel_x = *std::max_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());
        // double max_vel_y = *std::max_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());
        // double hmean = accumulate(ptcls.dprops["hp"].begin(), ptcls.dprops["hp"].end(), 0.0/ptcls.dprops["hp"].size());
        // double max_vel = std::max(std::sqrt(data.g * hmean)+max_vel_x,-std::sqrt(data.g * hmean)+max_vel_y);
        // //  double max_vel = std::max(1+max_vel_x,1+max_vel_y);
        // cel = std::abs(max_vel);
        if (it > 0)
          dt = 0.1 * data.hx / (1e-2 + cel); //0.2 *  data.hx / (1e-4 + cel);
        std::cout << "time = " << t << "  " << " dt = " <<  dt << std::endl;
	std::cout << "cel = " << cel << std::endl;
        my_timer.toc ("update dt");

	my_timer.tic ("save csv");
        std::string filename = "nc_particles_";
        filename = filename + std::to_string (it++);
        filename = filename + ".csv";
        if (it % 250 == 0)
        {
          std::ofstream OF (filename.c_str ());
          ptcls.print<particles_t::output_format::csv>(OF);
          OF.close ();
        }
	my_timer.toc ("save csv");

        //  (0)  CONNECTIVITY and BASIS FUNCTIONS
        my_timer.tic ("reorder");
        ptcls.init_particle_mesh ();

        ordering.resize (ptcls.num_particles);
        idx_t iordering = 0;
        for (auto const & ii : ptcls.grd_to_ptcl) {
          for (auto const & jj : ii.second) {
            ordering[iordering++] = jj;
          }
        }

        ptcls.reorder (ordering);
	iordering = 0;
        for (auto & ii : ptcls.grd_to_ptcl) {
          for (auto & jj : ii.second) {
            jj = iordering++;
          }
        }
	my_timer.toc ("reorder");

	my_timer.tic ("step 0");
        for (auto &v : vars) {
	  v.second.assign (v.second.size (), 0.0);
	}

        for (auto &v : Plotvars) {
	  v.second.assign (v.second.size (), 0.0);
	}




        ptcls.dprops.at("dZxp").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("dZyp").assign(ptcls.num_particles, 0.0);
        vars["Z"] = data.Z;
        vars["dZdx"] = data.dZdx;
        vars["dZdy"] = data.dZdy;

        my_timer.toc ("step 0");

        // (1) PROJECTION FROM MP TO NODES (P2G)
        my_timer.tic ("step 1");
        ptcls.p2g (vars,std::vector<std::string>{"Mp","mom_px","mom_py"},
                   std::vector<std::string>{"Mv","mom_vx","mom_vy"});
        ptcls.g2pd (vars,std::vector<std::string>{"Z"},
                    std::vector<std::string>{"dZxp"},
                    std::vector<std::string>{"dZyp"});

        my_timer.toc ("step 1");

        // (2)  EXTERNAL FORCES ON VERTICES (P2G)
        my_timer.tic ("step 2");

        std::transform (ptcls.dprops["Mp"].begin (), ptcls.dprops["Mp"].end (), ptcls.dprops["dZxp"].begin (),ptcls.dprops["Fpx"].begin (),
        [&] (double mp, double gzx) { return - data.g * mp * gzx; });
        std::transform (ptcls.dprops["Mp"].begin (), ptcls.dprops["Mp"].end (), ptcls.dprops["dZyp"].begin (),ptcls.dprops["Fpy"].begin (),
        [&] (double mp, double gzy) { return - data.g * mp * gzy; });

        std::transform (ptcls.dprops["Ap"].begin (), ptcls.dprops["Ap"].end (), ptcls.dprops["Fb_x"].begin (),
        ptcls.dprops["Fric_px"].begin (), std::multiplies<double> ());
        std::transform (ptcls.dprops["Ap"].begin (), ptcls.dprops["Ap"].end (), ptcls.dprops["Fb_y"].begin (),
        ptcls.dprops["Fric_py"].begin (), std::multiplies<double> ());

        ptcls.p2g (vars,std::vector<std::string>{"Fpx","Fpy","Fric_px","Fric_py"},
        std::vector<std::string>{"FPxv","FPyv","Fric_x","Fric_y"});

       std::transform (vars["FPxv"].begin (), vars["FPxv"].end (), vars["Fric_x"].begin (),
       vars["F_ext_vx"].begin (), std::plus<double> ());
       std::transform (vars["FPyv"].begin (), vars["FPyv"].end (), vars["Fric_y"].begin (),
       vars["F_ext_vy"].begin (), std::plus<double> ());


    /*    for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell) {
	  if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
	    for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
	      auto gv = icell -> gt (inode);
	      for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip) {
		auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];
		ptcls.dprops["Fric_px"][gp] =  - ptcls.dprops["Mp"][gp] *    9.81 * ptcls.dprops["dZxp"][gp]; // vars["dZdx"][gv]; // - ptcls.dprops["Ap"][gp] * ptcls.dprops["Fb_x"][gp] * 0.0 ;
		ptcls.dprops["Fric_py"][gp] =    - ptcls.dprops["Mp"][gp] *    9.81 * ptcls.dprops["dZyp"][gp]; // -  ptcls.dprops["Ap"][gp] *  ptcls.dprops["Fb_y"][gp] * 0.0;
	      }
	    }
	  }
	} */

    //  ptcls.p2g (vars,std::vector<std::string>{"Fric_px","Fric_py"},std::vector<std::string>{"Fric_x","Fric_y"});
      //  ptcls.p2g (vars,std::vector<std::string>{"Fric_px","Fric_py"},
                //   std::vector<std::string>{"F_ext_vx","F_ext_vy"});
  /*    std::transform (vars["Mv"].begin (), vars["Mv"].end (), vars["dZdx"].begin (), vars["Fric_x"].begin (), vars["F_ext_vx"].begin (),
      [&] (double mv, double gzx, double frx) { return - data.g * mv * gzx - frx; });
      std::transform (vars["Mv"].begin (), vars["Mv"].end (), vars["dZdy"].begin (), vars["Fric_y"].begin (), vars["F_ext_vy"].begin (),
      [&] (double mv, double gzy, double fry) { return - data.g * mv * gzy - fry; }); */

          /*   for (auto icell = grid.begin_cell_sweep ();
                icell != grid.end_cell_sweep (); ++icell)
                {
                for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
                {
                auto gv = icell -> gt (inode);

                vars["F_ext_vx"][gv] = -  vars["Mv"][gv] *  9.81 * vars["dZdx"][gv] - vars["Fric_x"][gv]; // vars["dZdx"][gv]; // vars["dZdx"][gv]; //  dZdx[gv]  ; //-  ptcls.dprops["Ap"][gp] * ptcls.dprops["Fb_x"][gp];

                vars["F_ext_vy"][gv] = -  vars["Mv"][gv] *  9.81 *  vars["dZdy"][gv] - vars["Fric_y"][gv]; // vars["dZdy"][gv]; //  -   ptcls.dprops["Ap"][gp] *  ptcls.dprops["Fb_y"][gp];

                }

              } */

        my_timer.toc ("step 2");

        // (3) INTERNAL FORCES (p2gd) and MOMENTUM BALANCE
        my_timer.tic ("step 3");
        ptcls.p2gd (vars, std::vector<std::string>{"F_11","F_21"},
                    std::vector<std::string>{"F_12","F_22"},
                    "Vp",std::vector<std::string>{"F_int_vx","F_int_vy"});



	std::transform (vars["F_ext_vx"].begin (), vars["F_ext_vx"].end (), vars["F_int_vx"].begin (), vars["Ftot_vx"].begin (), std::minus<double> ());
	std::transform (vars["F_ext_vy"].begin (), vars["F_ext_vy"].end (), vars["F_int_vy"].begin (), vars["Ftot_vy"].begin (), std::minus<double> ());

        // for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell) {
	//   if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
	//       auto iv = icell -> gt (inode);
	//       vars["Ftot_vx"][iv] = vars["F_ext_vx"][iv] + vars["F_int_vx"][iv];
	//       vars["Ftot_vy"][iv] = vars["F_ext_vy"][iv] + vars["F_int_vy"][iv];
	//     }
	//   }
	// }

	std::transform (vars["Ftot_vx"].begin (), vars["Ftot_vx"].end (), vars["Ftot_vx"].begin (), vars["mom_vx"].begin (), [=] (double x, double y) { return dt*x + y; });
	std::transform (vars["Ftot_vy"].begin (), vars["Ftot_vy"].end (), vars["Ftot_vy"].begin (), vars["mom_vy"].begin (), [=] (double x, double y) { return dt*x + y; });

        // for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell) {
	//   if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
	//       auto iv = icell -> gt (inode);
	//       vars["mom_vx"][iv] += dt * vars["Ftot_vx"][iv];
	//       vars["mom_vy"][iv] += dt * vars["Ftot_vy"][iv];
	//     }
	//   }
	// }

        my_timer.toc ("step 3");

        // (4)  COMPUTE NODAL ACCELERATIONS AND VELOCITIES
        my_timer.tic ("step 4");

	std::transform (vars["Ftot_vx"].begin (), vars["Ftot_vx"].end (), vars["Mv"].begin (), vars["avx"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	std::transform (vars["Ftot_vy"].begin (), vars["Ftot_vy"].end (), vars["Mv"].begin (), vars["avy"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	std::transform (vars["mom_vx"].begin (), vars["mom_vx"].end (), vars["Mv"].begin (), vars["vvx"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	std::transform (vars["mom_vy"].begin (), vars["mom_vy"].end (), vars["Mv"].begin (), vars["vvy"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	
        // for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell) {
	//   if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
	//       auto iv = icell -> gt (inode);
	//       vars["avx"][iv] = vars["Mv"][iv] > 1e-6 ?  vars["Ftot_vx"][iv] / vars["Mv"][iv] : 0.0;
	//       vars["avy"][iv] = vars["Mv"][iv] > 1e-6 ? vars[ "Ftot_vy"][iv] / vars["Mv"][iv] : 0.0;

	//       vars["vvx"][iv] = vars["Mv"][iv] > 1e-6 ?  vars["mom_vx"][iv] / vars["Mv"][iv] : 0.0;
	//       vars["vvy"][iv] = vars["Mv"][iv] > 1e-6 ?  vars["mom_vy"][iv] / vars["Mv"][iv] : 0.0;
	//     }
	//   }
	// }
        my_timer.toc ("step 4");

        // (5) BOUNDARY CONDITIONS - TO DO
        my_timer.tic ("step 5");
        /*      for (auto icell = grid.begin_cell_sweep ();
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
                } */
        my_timer.toc ("step 5");

        // (6) RETURN TO POINTS (G2P) and UPDATE POS AND VEL ON PARTICLES
        my_timer.tic ("step 6");
        ptcls.dprops.at("vpx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("vpy").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apy").assign(ptcls.num_particles, 0.0);

        ptcls.g2p (vars,std::vector<std::string>{"vvx","vvy","avx","avy"},
                   std::vector<std::string>{"vpx","vpy","apx","apy"});


	std::transform (ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (),  ptcls.dprops["apx"].begin (), ptcls.dprops["vpx"].begin (), [=] (double x, double y) { return x + dt * y; } );
	std::transform (ptcls.dprops["vpy"].begin (), ptcls.dprops["vpy"].end (),  ptcls.dprops["apy"].begin (), ptcls.dprops["vpy"].begin (), [=] (double x, double y) { return x + dt * y; } );
        // for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell) {
	//   if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
	//       for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip) {
	// 	auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];
	// 	ptcls.dprops["vpx"][gp] += dt * ptcls.dprops["apx"][gp] ;
	// 	ptcls.dprops["vpy"][gp] += dt * ptcls.dprops["apy"][gp] ;
	//       }
	//     }
	//   }
	// }


	std::transform (ptcls.x.begin (), ptcls.x.end (),  ptcls.dprops["vpx"].begin (), ptcls.x.begin (), [=] (double x, double y) { return x + dt * y; } );
	std::transform (ptcls.y.begin (), ptcls.y.end (),  ptcls.dprops["vpy"].begin (), ptcls.y.begin (), [=] (double x, double y) { return x + dt * y; } );
        // for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell) {
	//   if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
	//       for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip) {
	// 	auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];
	// 	ptcls.x[gp] += dt * ptcls.dprops["vpx"][gp];
	// 	ptcls.y[gp] += dt * ptcls.dprops["vpy"][gp];
	// 	//  ptcls.dprops["xp"][gp] += dt * ptcls.dprops["vpx"][gp];
	// 	//    ptcls.dprops["yp"][gp] += dt * ptcls.dprops["vpy"][gp];
	//       }
	//     }
	//   }
	// }

        my_timer.toc ("step 6");

        // (7) COMPUTE HEIGHT WITH STRAIN (divergence of velocities)
        my_timer.tic ("step 7");
        ptcls.dprops.at("vpx_dx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("vpy_dy").assign(ptcls.num_particles, 0.0);
        ptcls.g2pd (vars,std::vector<std::string>{"vvx","vvy"},
                    std::vector<std::string>{"vpx_dx","vpy_dx"},
                    std::vector<std::string>{"vpx_dy","vpy_dy"});


	std::transform (ptcls.dprops["vpx_dx"].begin (), ptcls.dprops["vpx_dx"].end (), ptcls.dprops["vpy_dy"].begin (), norm_v.begin (), std::plus<double> ());
	std::transform (ptcls.dprops["hp"].begin (), ptcls.dprops["hp"].end (), norm_v.begin (), ptcls.dprops["hp"].begin (), [=] (double x, double y) { return x / (1 + dt * y); } );
	std::transform (ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (), ptcls.dprops["Mp"].begin (), ptcls.dprops["mom_px"].begin (), std::multiplies<double> () );
	std::transform (ptcls.dprops["vpy"].begin (), ptcls.dprops["vpy"].end (), ptcls.dprops["Mp"].begin (), ptcls.dprops["mom_py"].begin (), std::multiplies<double> () );
	std::transform (ptcls.dprops["Vp"].begin (), ptcls.dprops["Vp"].end (), norm_v.begin (), ptcls.dprops["Vp"].begin (), [=] (double x, double y) { return x / (1 + dt * y); } );
  std::transform (ptcls.dprops["Vp"].begin (), ptcls.dprops["Vp"].end (), ptcls.dprops["hp"].begin (), ptcls.dprops["Ap"].begin (), std::divides<double> () );
//	std::transform (ptcls.dprops["Mp"].begin (), ptcls.dprops["Mp"].end (), ptcls.dprops["hp"].begin (), ptcls.dprops["Ap"].begin (), [&] (double x, double y) { return x / (data.rho * y); } );
        // for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell) {
	//   if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
	//       for (auto gp = 0; gp < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++gp) {
	// 	auto ip = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[gp];
	// 	ptcls.dprops["hp"][ip] /= (1 + dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip]));
	// 	ptcls.dprops["mom_px"][ip] = ptcls.dprops["vpx"][ip] * ptcls.dprops["Mp"][ip];
	// 	ptcls.dprops["mom_py"][ip] = ptcls.dprops["vpy"][ip] * ptcls.dprops["Mp"][ip];
	// 	ptcls.dprops["Vp"][ip] /= (1 + dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip]));
	// 	ptcls.dprops["Ap"][ip] = ptcls.dprops["Mp"][ip] / (data.rho * ptcls.dprops["hp"][ip]);
	//       }
	//     }
	//   }
	// }

        my_timer.toc ("step 7");

// (7b) UPDATE FRICTIONS
my_timer.tic ("step 7b");
	std::transform (ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (), ptcls.dprops["vpy"].begin (), norm_v.begin (), [] (double x, double y) { return std::sqrt (x*x + y*y); });
	std::transform (ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (), ptcls.dprops["hp"].begin (), ptcls.dprops["Fb_x"].begin (),
			[&] (double x, double y) { return - data.rho * data.g * ( y * std::tan(fric_ang) + x * x  / data.xi) * x; });
	std::transform (ptcls.dprops["Fb_x"].begin (), ptcls.dprops["Fb_x"].end (), norm_v.begin (), ptcls.dprops["Fb_x"].begin (), std::divides<double> ());
	std::transform (ptcls.dprops["vpy"].begin (), ptcls.dprops["vpy"].end (), ptcls.dprops["hp"].begin (), ptcls.dprops["Fb_y"].begin (),
			[&] (double x, double y) { return - data.rho * data.g * ( y * std::tan(fric_ang) + x * x  / data.xi) * x; });
	std::transform (ptcls.dprops["Fb_y"].begin (), ptcls.dprops["Fb_y"].end (), norm_v.begin (), ptcls.dprops["Fb_y"].begin (), std::divides<double> ());
        // for (idx_t ip = 0; ip < num_particles; ++ip)
        //   {
        //     norm_v[ip] = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip] );
	//
        //     ptcls.dprops["Fb_x"][ip] = - ( data.rho * ptcls.dprops["hp"][ip] * data.g * std::tan(fric_ang)  +
        //                                    data.rho * data.g * ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip]  / data.xi) * ptcls.dprops["vpx"][ip] / (norm_v[ip] + 0.001) ;
	//
        //     ptcls.dprops["Fb_y"][ip] = - ( data.rho * ptcls.dprops["hp"][ip] * data.g * std::tan(fric_ang)  +
        //                                    data.rho * data.g * ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip]  / data.xi)  * ptcls.dprops["vpy"][ip] / (norm_v[ip] + 0.001) ;
        //   }


        //        ptcls.p2g (Plotvars,std::vector<std::string>{"Mp","vpx","vpy","apx","apy"},
        //        std::vector<std::string>{"rho_v","vvx","vvy","avx","avy"});
   my_timer.toc ("step 7b");


   // (8) UPDATE PARTICLE STRESS (USL)
        my_timer.tic ("step 8");

	// std::transform (ptcls.dprops["hp"].begin (), ptcls.dprops["hp"].end (), ptcls.dprops["F_11"].begin (), [&] (double x) { return .5 * data.rho * data.g * x; });
	// ptcls.dprops["F_12"].assign (ptcls.num_particles, 0.0);
	// ptcls.dprops["F_21"].assign (ptcls.num_particles, 0.0);
	// std::transform (ptcls.dprops["hp"].begin (), ptcls.dprops["hp"].end (), ptcls.dprops["F_22"].begin (), [&] (double x) { return .5 * data.rho * data.g * x; });


        // ptcls.dprops.at("hpZ").assign(ptcls.num_particles, 0.0);
        // ptcls.dprops.at("Zp").assign(ptcls.num_particles, 0.0);

        // ptcls.g2p (vars,std::vector<std::string>{"Z"},
        //            std::vector<std::string>{"Zp"});

	// std::transform (ptcls.dprops["hp"].begin (), ptcls.dprops["hp"].end (), ptcls.dprops["Zp"].begin (), ptcls.dprops["hpZ"].begin (), std::plus<double> ());


//============= DA OTTIMIZZARE TENSORE SIGMA!! *************=====================
//==============================================================================

	// //std::cout << "\n\n\n***********************************\n";
	for (idx_t ip = 0; ip < ptcls.num_particles; ++ip) {
	  nrm = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip] );
	  //std::cout << nrm << " ";
	  
	  ALF = ptcls.dprops["hp"][ip] > 1.e-3
	    ? (6. * 50. * nrm)/((ptcls.dprops["hp"][ip]+0.001) * 2000.)
	    : 0.0;
	  B = -114./32. - ALF;
	  Z1 = (-B + std::sqrt(B * B - 4. * A * C))/(2. * A);
	  Z2 = (-B - std::sqrt(B * B - 4. * A * C))/(2. * A);
	  ZZ = std::abs (Z1 - .5) <= .5 ? Z1 : Z2;
	  //std::cout << ZZ << " ";
 
	  s_xx = ptcls.dprops["vpx_dx"][ip];
	  s_xy = 0.5 * (ptcls.dprops["vpx_dy"][ip] + ptcls.dprops["vpy_dx"][ip]);
	  s_yy = ptcls.dprops["vpy_dy"][ip];

	  D_xx = s_xx;
	  D_yx = s_xy;
	  D_zx = ptcls.dprops["hp"][ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ)) * (ptcls.dprops["vpx"][ip] / (ptcls.dprops["hp"][ip] + 0.001) ) : 0.0;

	  D_xy = s_xy;
	  D_yy = s_yy;
	  D_zy = ptcls.dprops["hp"][ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ))  * (ptcls.dprops["vpy"][ip] / (ptcls.dprops["hp"][ip]+ 0.001) ) : 0.0;

	  D_xz =  0.0;//0.0; // 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpx"][ip] / (ptcls.dprops["hp"][ip] + 0.001) );
	  D_yz =  0.0;// 0.0; // 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpy"][ip] / (ptcls.dprops["hp"][ip]+ 0.001) );
	  D_zz = - (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip] );

	  invII = 0.5 * (D_xx * D_xx + D_yy * D_yy + D_zz * D_zz +
			 D_xz * D_xz + D_yz * D_yz + D_xy * D_xy);
	  //std::cout << invII << " ";
	  
	  sig_xx = invII != 0 ? (2000./std::sqrt(invII) + 2. * 50.) * D_xx : 0.0;
	  sig_xy = invII != 0 ? (2000./std::sqrt(invII) + 2. * 50.) * D_xy : 0.0;
	  sig_yy = invII != 0 ? (2000./std::sqrt(invII) + 2. * 50.) * D_yy : 0.0;
	  // std::cout << sig_xx << " ";
	  // std::cout << sig_xy << " ";
	  // std::cout << sig_yy << " ";
	  // std::cout << "\n";
	  
	  ptcls.dprops["F_11"][ip] =  cc * sig_xx - .5 * data.rho * data.g *  (ptcls.dprops["hp"][ip]);
	  ptcls.dprops["F_12"][ip] =  cc * sig_xy;
	  ptcls.dprops["F_21"][ip] =  cc * sig_xy;
	  ptcls.dprops["F_22"][ip] =  cc * sig_yy - .5 * data.rho * data.g *  (ptcls.dprops["hp"][ip]);
	  
	}

	ptcls.dprops.at("hpZ").assign(ptcls.num_particles, 0.0);
	ptcls.dprops.at("Zp").assign(ptcls.num_particles, 0.0);
	ptcls.g2p (vars,std::vector<std::string>{"Z"},
		   std::vector<std::string>{"Zp"});

	std::transform (ptcls.dprops["hp"].begin (), ptcls.dprops["hp"].end (), ptcls.dprops["Zp"].begin (), ptcls.dprops["hpZ"].begin (), std::plus<double> ());
	
	// for (idx_t ip = 0; ip < num_particles; ++ip)
	//   {
	//     sig_xx[ip] = 0.0;
	//     sig_xy[ip] = 0.0;
	//     sig_yy[ip] =  0.0;
	//   }

	// double double_dot = 0.0;

	// for (auto icell = grid.begin_cell_sweep ();
	//      icell != grid.end_cell_sweep (); ++icell)
        //   {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
	//       {
	// 	auto gv = icell -> gt (inode);
	// 	if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
	// 	  for (auto gp = 0; gp < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++gp)
	// 	    {


	// 	      auto ip = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[gp];
	// 	      norm_v[ip] = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip] );

	// 	      ALF[ip] =ptcls.dprops["hp"][ip] > 1.e-3 ? (6. * 50. * norm_v[ip])/((ptcls.dprops["hp"][ip]+0.001) * 2000.):0.0;
	// 	      B[ip] = -114./32. - ALF[ip];
	// 	      Z1[ip] = (-B[ip] + std::sqrt(B[ip] * B[ip] - 4. * A * C))/(2. * A);
	// 	      Z2[ip] = (-B[ip] - std::sqrt(B[ip] * B[ip] - 4. * A * C))/(2. * A);
	// 	      ZZ[ip] = std::abs(Z1[ip] - .5)<=.5 ? Z1[ip] : Z2[ip];

	// 	      s_xx[ip] = ptcls.dprops["vpx_dx"][ip];
	// 	      s_xy[ip] = 0.5 * (ptcls.dprops["vpx_dy"][ip] + ptcls.dprops["vpy_dx"][ip]);
	// 	      s_yy[ip] = ptcls.dprops["vpy_dy"][ip];

	// 	      D_xx[ip] = s_xx[ip];
	// 	      D_yx[ip] = s_xy[ip];
	// 	      D_zx[ip] = ptcls.dprops["hp"][ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpx"][ip] / (ptcls.dprops["hp"][ip] + 0.001) ) : 0.0;

	// 	      D_xy[ip] = s_xy[ip];
	// 	      D_yy[ip] = s_yy[ip];
	// 	      D_zy[ip] = ptcls.dprops["hp"][ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ[ip]))  * (ptcls.dprops["vpy"][ip] / (ptcls.dprops["hp"][ip]+ 0.001) ) : 0.0;

	// 	      D_xz[ip] =  0.0;//0.0; // 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpx"][ip] / (ptcls.dprops["hp"][ip] + 0.001) );
	// 	      D_yz[ip] =  0.0;// 0.0; // 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpy"][ip] / (ptcls.dprops["hp"][ip]+ 0.001) );
	// 	      D_zz[ip] = - (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip] );

	// 	      invII[ip] = 0.5 * (D_xx[ip] * D_xx[ip] + D_yy[ip] * D_yy[ip] + D_zz[ip] * D_zz[ip] +
	// 				 D_xz[ip] * D_xz[ip] + D_yz[ip] * D_yz[ip] + D_xy[ip] * D_xy[ip]);

	// 	      sig_xx[ip] = invII[ip] != 0 ? (2000./std::sqrt(invII[ip]) + 2. * 50.) * D_xx[ip] : 0.0;
	// 	      sig_xy[ip] = invII[ip] != 0 ? (2000./std::sqrt(invII[ip]) + 2. * 50.) * D_xy[ip] : 0.0;
	// 	      sig_yy[ip] = invII[ip] != 0 ? (2000./std::sqrt(invII[ip]) + 2. * 50.) * D_yy[ip] : 0.0;

	// 	      double cc = 1.0;

	// 	      ptcls.dprops["F_11"][ip] =  cc * sig_xx[ip] -  .5 * data.rho * data.g *  (ptcls.dprops["hp"][ip] ) ;
	// 	      ptcls.dprops["F_12"][ip] =  cc * sig_xy[ip];
	// 	      ptcls.dprops["F_21"][ip] =  cc * sig_xy[ip];
	// 	      ptcls.dprops["F_22"][ip] =  cc * sig_yy[ip] - .5 * data.rho * data.g *   (ptcls.dprops["hp"][ip]  );
	// 	    }
	//       }

        //   }

	// ptcls.dprops.at("hpZ").assign(ptcls.num_particles, 0.0);
	// ptcls.dprops.at("Zp").assign(ptcls.num_particles, 0.0);
	// ptcls.g2p (vars,std::vector<std::string>{"Z"},
	// 	   std::vector<std::string>{"Zp"});

	// for (auto icell = grid.begin_cell_sweep ();
	//      icell != grid.end_cell_sweep (); ++icell)
	//   {
	//     for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
	//       {
	// 	auto gv = icell -> gt (inode);
	// 	if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
	// 	  for (auto gp = 0; gp < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++gp)
	// 	    {
	// 	      auto ip = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[gp];
	// 	      ptcls.dprops["hpZ"][ip] = ptcls.dprops["hp"][ip] + ptcls.dprops["Zp"][ip];

	// 	    }
	//       }

	//   }  

    //******====================================================================
    //******====================================================================

        my_timer.toc ("step 8");

	my_timer.tic ("step 8b");
        ptcls.p2g (vars,std::vector<std::string>{"hp"}, std::vector<std::string>{"HV"});
	my_timer.toc ("step 8b");


	//        my_timer.tic ("save vts");
	// filename = "nc_grid_";
	// filename = filename + std::to_string (it);
	// filename = filename + ".vts";
	// if (it % 250 == 0)
	// {
        //  grid.vtk_export(filename.c_str(), Plotvars);
        //  grid.vtk_export(filename.c_str(), vars);
	// }
        t +=dt;
	// my_timer.toc ("save vts");

        my_timer.print_report ();

      }

    my_timer.print_report ();
    //  ptcls.print<particles_t::output_format::csv>(std::cout);

    return 0;
}
