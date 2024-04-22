#include <fstream>
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
  double fric_ang = 34. * M_PI / 180.;
  double atan_grad_z;
  std::vector<double> norm_v (num_particles, 0.0);
  double A = 3./2.;
  std::vector<double> ALF (num_particles, 0.0);
  std::vector<double> B (num_particles, -114./32.);
  double C = 65./32.;
  std::vector<double> s_xx (num_particles, 0.0);
  std::vector<double> s_xy (num_particles, 0.0);
  std::vector<double> s_yy (num_particles, 0.0);
  std::vector<double> D_xx (num_particles, 0.0);
  std::vector<double> D_xy (num_particles, 0.0);
  std::vector<double> D_xz (num_particles, 0.0);
  std::vector<double> D_yx (num_particles, 0.0);
  std::vector<double> D_yy (num_particles, 0.0);
  std::vector<double> D_yz (num_particles, 0.0);
  std::vector<double> D_zx (num_particles, 0.0);
  std::vector<double> D_zy (num_particles, 0.0);
  std::vector<double> D_zz (num_particles, 0.0);
  std::vector<double> invII (num_particles, 0.0);
  std::vector<double> Z1 (num_particles, 0.0);
  std::vector<double> Z2 (num_particles, 0.0);
  std::vector<double> ZZ (num_particles, 0.0);
  std::vector<double> sig_xx (num_particles, 0.0);
  std::vector<double> sig_xy (num_particles, 0.0);
  std::vector<double> sig_yy (num_particles, 0.0);




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
      {"HV", std::vector<double>(grid.num_global_nodes (), 0.0)}
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
    ptcls.dprops["F_11"][ip] = -  0.5 * data.rho * data.g *   (ptcls.dprops["hp"][ip]  ) ;
    ptcls.dprops["F_12"][ip] = 0.0;
    ptcls.dprops["F_21"][ip] = 0.0;
    ptcls.dprops["F_22"][ip] = - 0.5 * data.rho * data.g *   (ptcls.dprops["hp"][ip] );
  }

    int it = 0;

    ptcls.build_mass();


    grid.vtk_export("GRID_forZ.vts", vars);


    dt = 1.0e-5;
    std::vector<idx_t> ordering (ptcls.num_particles);
    while (t < 1.0) //data.T;
        {

  	my_timer.tic ("update dt");
    /*      double max_vel_x = *std::max_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());
            double max_vel_y = *std::max_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());
            double hmean = accumulate(ptcls.dprops["hp"].begin(), ptcls.dprops["hp"].end(), 0.0/ptcls.dprops["hp"].size());
            double max_vel = std::max(std::sqrt(data.g * hmean)+max_vel_x,-std::sqrt(data.g * hmean)+max_vel_y);
          //  double max_vel = std::max(1+max_vel_x,1+max_vel_y);
            cel = std::abs(max_vel); */

           double max_vel_x = *std::max_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());
            double min_vel_x = *std::min_element(ptcls.dprops["vpx"].begin(), ptcls.dprops["vpx"].end());

            max_vel_x = std::max (std::abs(max_vel_x), std::abs(min_vel_x));

            double max_vel_y = *std::max_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());
            double min_vel_y = *std::min_element(ptcls.dprops["vpy"].begin(), ptcls.dprops["vpy"].end());

            max_vel_y = std::max (std::abs(max_vel_y), std::abs(min_vel_y));

            double hmean = *std::max_element (ptcls.dprops["hp"].begin(), ptcls.dprops["hp"].end());
            double max_vel = std::max(std::sqrt(data.g * hmean) + max_vel_x, std::sqrt(data.g * hmean) + max_vel_y);
            cel = std::abs(max_vel);

	    if (it > 0)
	      dt = 0.2* data.hx / (1e-2 + cel); // 0.7      0.2 *  data.hx / (1e-4 + cel);
	      std::cout << "time = " << t << "  " << " dt = " <<  dt << std::endl;
    my_timer.toc ("update dt");

    my_timer.tic ("save csv");

            std::string filename = "nc_particles_";
            filename = filename + std::to_string (it++);
            filename = filename + ".csv";
            if (t>=0.0)//(it % 50 == 0)
            {
              std::ofstream OF (filename.c_str ());
              ptcls.print<particles_t::output_format::csv>(OF);
              OF.close ();
            }
    my_timer.toc("save csv");

	//  (0)  CONNECTIVITY and BASIS FUNCTIONS
	my_timer.tic ("reorder");
        ptcls.init_particle_mesh ();

  /*      ordering.resize (ptcls.num_particles);
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
    my_timer.toc ("reorder"); */

    my_timer.tic("step 0");
        for (auto &v : vars)
	  {
            v.second.assign (v.second.size (),0.0);
	  }

        for (auto &v : Plotvars)
	  {
            v.second.assign (v.second.size (),0.0);
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
	for (auto icell = grid.begin_cell_sweep (); icell != grid.end_cell_sweep (); ++icell)
	  {
	    for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
	      {
		auto gv = icell -> gt (inode);
   	if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
		  for (auto ip = 0; ip < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ip)
		     {
		      auto gp = ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ip];
		      ptcls.dprops["Fric_px"][gp] = 0.0 * ptcls.dprops["Ap"][gp] *  ptcls.dprops["Fb_x"][gp]  ; // - ptcls.dprops["Mp"][gp] *    9.81 * ptcls.dprops["dZxp"][gp] + ptcls.dprops["Ap"][gp] * ptcls.dprops["Fb_x"][gp]  ;
		      ptcls.dprops["Fric_py"][gp] = 0.0 * ptcls.dprops["Ap"][gp] * ptcls.dprops["Fb_y"][gp]  ;//  - ptcls.dprops["Mp"][gp] *    9.81 * ptcls.dprops["dZyp"][gp]  + ptcls.dprops["Ap"][gp] *  ptcls.dprops["Fb_y"][gp]  ;
		     }
	     }

	  }

	ptcls.p2g (vars,std::vector<std::string>{"Fric_px","Fric_py"},
		  std::vector<std::string>{"Fric_x","Fric_y"});
//    ptcls.p2g (vars,std::vector<std::string>{"Fric_px","Fric_py"},
//  		   std::vector<std::string>{"F_ext_vx","F_ext_vy"});

	for (auto icell = grid.begin_cell_sweep ();
             icell != grid.end_cell_sweep (); ++icell)
	  {
	    for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
	      {
		auto gv = icell -> gt (inode);

		vars["F_ext_vx"][gv] = -  vars["Mv"][gv] *  9.81 * 0.0 * vars["dZdx"][gv] + vars["Fric_x"][gv]; // vars["dZdx"][gv]; // vars["dZdx"][gv]; //  dZdx[gv]  ; //-  ptcls.dprops["Ap"][gp] * ptcls.dprops["Fb_x"][gp];

		vars["F_ext_vy"][gv] = -  vars["Mv"][gv] *  9.81 * 0.0 *  vars["dZdy"][gv] + vars["Fric_y"][gv]; // vars["dZdy"][gv]; //  -   ptcls.dprops["Ap"][gp] *  ptcls.dprops["Fb_y"][gp];

	      }

	  }

	my_timer.toc ("step 2");

	// (3) INTERNAL FORCES (p2gd) and MOMENTUM BALANCE
	my_timer.tic ("step 3");
        ptcls.p2gd (vars, std::vector<std::string>{"F_11","F_21"},
		    std::vector<std::string>{"F_12","F_22"},
		    "Vp",std::vector<std::string>{"F_int_vx","F_int_vy"});



	for (auto icell = grid.begin_cell_sweep ();
	     icell != grid.end_cell_sweep (); ++icell)
	  {
            for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
	      {
		auto iv = icell -> gt (inode);
		vars["Ftot_vx"][iv] = vars["F_ext_vx"][iv] - vars["F_int_vx"][iv];
		vars["Ftot_vy"][iv] = vars["F_ext_vy"][iv] - vars["F_int_vy"][iv];
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

	      }

	  }
	my_timer.toc ("step 3");



	// (4)  COMPUTE NODAL ACCELERATIONS AND VELOCITIES
	my_timer.tic ("step 4");
	for (auto icell = grid.begin_cell_sweep ();
	     icell != grid.end_cell_sweep (); ++icell)
	  {
	    for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode)
	      {
		auto iv = icell -> gt (inode);
		vars["avx"][iv] = vars["Mv"][iv] > 1e-8 ?  vars["Ftot_vx"][iv] / vars["Mv"][iv] : 0.0;
		vars["avy"][iv] = vars["Mv"][iv] > 1e-8 ? vars[ "Ftot_vy"][iv] / vars["Mv"][iv] : 0.0;

		vars["vvx"][iv] = vars["Mv"][iv] > 1e-8 ?  vars["mom_vx"][iv] / vars["Mv"][iv] : 0.0;
		vars["vvy"][iv]= vars["Mv"][iv] > 1e-8 ?  vars["mom_vy"][iv] / vars["Mv"][iv] : 0.0;
	      }

	  }
	my_timer.toc ("step 4");

  // (5) BOUNDARY CONDITIONS - TO DO
 my_timer.tic ("step 5");
  for (auto icell = grid.begin_cell_sweep ();
      icell != grid.end_cell_sweep (); ++icell)
   {
        if ( (icell->e (2) == 2) || (icell->e(3)==3)  )
          {
          for (idx_t inode = 0; inode < 4; ++inode)
          {
                vars["avx"][icell->gt(inode)] = 0.0;
              vars["vvx"][icell->gt(inode)] = 0.0;
              //  vars["mom_vx"][icell->gt(inode)] = 0.0;
          }
          }
       if ( (icell->e(0) == 0) || (icell->e(1)==1) )
         {
          for (idx_t inode = 0; inode < 4; ++inode)
               {
                 vars["avy"][icell->gt(inode)] = 0.0;
                vars["vvy"][icell->gt(inode)] = 0.0;
                //   vars["mom_vy"][icell->gt(inode)] = 0.0;
               }
         }
      }


 my_timer.toc ("step 5");



	// (6) RETURN TO POINTS (G2P) and UPDATE POS AND VEL ON PARTICLES
	my_timer.tic ("step 6");
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
	my_timer.toc ("step 6");

	// (7) COMPUTE HEIGHT WITH STRAIN (divergence of velocities)
	my_timer.tic ("step 7");
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
		      ptcls.dprops["Ap"][ip]   /=(1+dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip]));
ptcls.dprops["Vp"][ip] = ptcls.dprops["hp"][ip] * ptcls.dprops["Ap"][ip];
        // ptcls.dprops["Ap"][ip] =  ptcls.dprops["hp"][ip] >= 1.e-2 ? ptcls.dprops["Vp"][ip]/ptcls.dprops["hp"][ip] : ptcls.dprops["hp"][ip]; //     (1+dt * (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip])); // ptcls.dprops["Vp"][ip]/ptcls.dprops["hp"][ip]; //ptcls.dprops["Mp"][ip] / (data.rho * ptcls.dprops["hp"][ip]);
		    }
	      }

	  }




	my_timer.toc ("step 7");

  for (idx_t ip = 0; ip < num_particles; ++ip)
    {
            norm_v[ip] = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip] );

            ptcls.dprops["Fb_x"][ip] = ptcls.dprops["hp"][ip] > 1.e-3 ?  - 0.0 * 2.5 * (  data.rho * ptcls.dprops["hp"][ip] * data.g * std::tan(fric_ang)  +
                                           data.rho * data.g * ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip]  / 100.)* ptcls.dprops["vpx"][ip] / (norm_v[ip] + 0.001) : 0.0 ;

            ptcls.dprops["Fb_y"][ip] = ptcls.dprops["hp"][ip] > 1.e-3 ? - 0.0 * 2.5 * (data.rho * ptcls.dprops["hp"][ip] * data.g * std::tan(fric_ang)  +
                                           data.rho * data.g * ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip]  / 100.)  * ptcls.dprops["vpy"][ip] / (norm_v[ip] + 0.001) : 0.0 ;
    }

	// (8) UPDATE PARTICLE STRESS (USL)
	my_timer.tic ("step 8");

  for (idx_t ip = 0; ip < num_particles; ++ip)
    {
      sig_xx[ip] = 0.0;
      sig_xy[ip] = 0.0;
      sig_yy[ip] =  0.0;
    }

 double double_dot = 0.0;

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
          norm_v[ip] = std::sqrt(ptcls.dprops["vpx"][ip] * ptcls.dprops["vpx"][ip] + ptcls.dprops["vpy"][ip] * ptcls.dprops["vpy"][ip] );

          ALF[ip] =ptcls.dprops["hp"][ip] > 1.e-3 ? (6. * 50. * norm_v[ip])/((ptcls.dprops["hp"][ip]+0.001) * 2000.):0.0;
          B[ip] = -114./32. - ALF[ip];
          Z1[ip] = (-B[ip] + std::sqrt(B[ip] * B[ip] - 4. * A * C))/(2. * A);
          Z2[ip] = (-B[ip] - std::sqrt(B[ip] * B[ip] - 4. * A * C))/(2. * A);
          Z2[ip] = std::abs(Z1[ip] - .5)<=.5 ? Z1[ip] : Z2[ip];

          s_xx[ip] = ptcls.dprops["vpx_dx"][ip];
          s_xy[ip] = 0.5 * (ptcls.dprops["vpx_dy"][ip] + ptcls.dprops["vpy_dx"][ip]);
          s_yy[ip] = ptcls.dprops["vpy_dy"][ip];

          D_xx[ip] = s_xx[ip];
          D_yx[ip] = s_xy[ip];
          D_zx[ip] = ptcls.dprops["hp"][ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpx"][ip] / (ptcls.dprops["hp"][ip] + 0.001) ) : 0.0;

          D_xy[ip] = s_xy[ip];
          D_yy[ip] = s_yy[ip];
          D_zy[ip] = ptcls.dprops["hp"][ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ[ip]))  * (ptcls.dprops["vpy"][ip] / (ptcls.dprops["hp"][ip]+ 0.001) ) : 0.0;

          D_xz[ip] =  0.0;//0.0; // 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpx"][ip] / (ptcls.dprops["hp"][ip] + 0.001) );
          D_yz[ip] =  0.0;// 0.0; // 0.5 * (3. / (2. + ZZ[ip])) * (ptcls.dprops["vpy"][ip] / (ptcls.dprops["hp"][ip]+ 0.001) );
          D_zz[ip] = - (ptcls.dprops["vpx_dx"][ip] + ptcls.dprops["vpy_dy"][ip] );

          invII[ip] = 0.5 * (D_xx[ip] * D_xx[ip] + D_yy[ip] * D_yy[ip] + D_zz[ip] * D_zz[ip] +
                      D_xz[ip] * D_xz[ip] + D_yz[ip] * D_yz[ip] + D_xy[ip] * D_xy[ip]);

          sig_xx[ip] = invII[ip] != 0 ? (2000./std::sqrt(invII[ip]) + 2. * 50.) * D_xx[ip] : 0.0;
          sig_xy[ip] = invII[ip] != 0 ? (2000./std::sqrt(invII[ip]) + 2. * 50.) * D_xy[ip] : 0.0;
          sig_yy[ip] = invII[ip] != 0 ? (2000./std::sqrt(invII[ip]) + 2. * 50.) * D_yy[ip] : 0.0;

          double cc = 0.0;

		      ptcls.dprops["F_11"][ip] =  cc * sig_xx[ip] - 0.5 * data.rho * data.g *  (ptcls.dprops["hp"][ip]  ) ;
		      ptcls.dprops["F_12"][ip] =  cc * sig_xy[ip] ;
		      ptcls.dprops["F_21"][ip] =  cc * sig_xy[ip]  ;
		      ptcls.dprops["F_22"][ip] =  cc * sig_yy[ip] - 0.5 * data.rho * data.g *   (ptcls.dprops["hp"][ip]    )   ;
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


 my_timer.toc("step 8");

	        ptcls.p2g (Plotvars,std::vector<std::string>{"Mp","vpx","vpy","apx","apy"},
             std::vector<std::string>{"rho_v","vvx","vvy","avx","avy"},true);

 my_timer.tic("step 8b");
        ptcls.p2g (vars,std::vector<std::string>{"hp"}, std::vector<std::string>{"HV"});

	my_timer.toc ("step 8b");

	my_timer.tic ("save vts");
        filename = "nc_grid_";
        filename = filename + std::to_string (it);
      filename = filename + ".vts";
	 // grid.vtk_export(filename.c_str(), Plotvars);
	grid.vtk_export(filename.c_str(), vars);
        t +=dt;
	my_timer.toc ("save vts");

      }

    my_timer.print_report ();
    //  ptcls.print<particles_t::output_format::csv>(std::cout);

    return 0;
}
