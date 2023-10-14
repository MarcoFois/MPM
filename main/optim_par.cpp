#include "counter.h"
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/tbb/global_control.h>
#include <fstream>
#include <functional>
#include <map>
#include <cmath>
#include <iostream>
#include <particles.h>
#include <quadgrid_cpp.h>
#include <timer.h>
#include "mpm_data.h"


struct stress_tensor_t {

  const std::vector<double>& vpx;
  const std::vector<double>& vpy;
  const std::vector<double>& hp;
  const std::vector<double>& vpx_dx;
  const std::vector<double>& vpy_dy;
  const std::vector<double>& vpx_dy;
  const std::vector<double>& vpy_dx;
  std::vector<double>& F_11;
  std::vector<double>& F_12;
  std::vector<double>& F_21;
  std::vector<double>& F_22;
  const DATA& data;
  
  stress_tensor_t (particles_t &ptcls, DATA &data_) :
    vpx{ptcls.dprops.at ("vpx")},
    vpy{ptcls.dprops.at ("vpy")},
    hp{ptcls.dprops.at ("hp")},
    vpx_dx{ptcls.dprops.at ("vpx_dx")},
    vpy_dy{ptcls.dprops.at ("vpy_dy")},
    vpx_dy{ptcls.dprops.at ("vpx_dy")},
    vpy_dx{ptcls.dprops.at ("vpy_dx")},
    F_11{ptcls.dprops.at ("F_11")},
    F_12{ptcls.dprops.at ("F_12")},
    F_21{ptcls.dprops.at ("F_21")},
    F_22{ptcls.dprops.at ("F_22")},
    data{data_} { }

  
  void operator()(int ip) {

    double A{3./2.}, B{-114./32.}, C{65./32.};
    double nrm{0.0}, ALF{0.0}, Z1{0.0}, Z2{0.0}, ZZ{0.0};
    double s_xx{0.0}, s_xy{0.0}, s_yy{0.0};
    double D_xx{0.0}, D_yx{0.0}, D_zx{0.0}, D_xy{0.0}, D_yy{0.0}, D_zy{0.0}, D_xz{0.0}, D_yz{0.0}, D_zz{0.0};
    double invII{0.0}, sig_yy{0.0}, sig_xy{0.0}, sig_xx{0.0}, cc{1.0};
      
    nrm = std::sqrt(vpx[ip] * vpx[ip] +vpy[ip] * vpy[ip] );
	  
    ALF = hp[ip] > 1.e-3
      ? (6. * 50. * nrm)/((hp[ip]+0.001) * 2000.)
      : 0.0;
    B = -114./32. - ALF;
    Z1 = (-B + std::sqrt(B * B - 4. * A * C))/(2. * A);
    Z2 = (-B - std::sqrt(B * B - 4. * A * C))/(2. * A);
    ZZ = std::abs (Z1 - .5) <= .5 ? Z1 : Z2;
 
    s_xx = vpx_dx[ip];
    s_xy = 0.5 * (vpx_dy[ip] + vpy_dx[ip]);
    s_yy = vpy_dy[ip];

    D_xx = s_xx;
    D_yx = s_xy;
    D_zx = hp[ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ)) * (vpx[ip] / (hp[ip] + 0.001) ) : 0.0;

    D_xy = s_xy;
    D_yy = s_yy;
    D_zy = hp[ip] > 1.e-3 ? 0.5 * (3. / (2. + ZZ))  * (vpy[ip] / (hp[ip]+ 0.001) ) : 0.0;

    D_xz =  0.0;
    D_yz =  0.0;
    D_zz = - (vpx_dx[ip] + vpy_dy[ip] );

    invII = 0.5 * (D_xx * D_xx + D_yy * D_yy + D_zz * D_zz +
		   D_xz * D_xz + D_yz * D_yz + D_xy * D_xy);

    sig_xx = invII != 0 ? (2000./std::sqrt(invII) + 2. * 50.) * D_xx : 0.0;
    sig_xy = invII != 0 ? (2000./std::sqrt(invII) + 2. * 50.) * D_xy : 0.0;
    sig_yy = invII != 0 ? (2000./std::sqrt(invII) + 2. * 50.) * D_yy : 0.0;

    F_11[ip] =  cc * sig_xx - .5 * data.rho * data.g *  (hp[ip]);
    F_12[ip] =  cc * sig_xy;
    F_21[ip] =  cc * sig_xy;
    F_22[ip] =  cc * sig_yy - .5 * data.rho * data.g *  (hp[ip]);
	  
  }
  
};

using idx_t = quadgrid_t<std::vector<double>>::idx_t;
cdf::timer::timer_t my_timer{};
auto policy = std::execution::unseq;

int main ()
{

  //tbb::global_control(tbb::global_control::max_allowed_parallelism, 1);
 
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

    my_timer.tic ("g2p");
    ptcls.g2p (vars, {"Z"}, {"Zp"});
    my_timer.toc ("g2p");
    
    ptcls.dprops.at("dZxp").assign(ptcls.num_particles, 0.0);
    ptcls.dprops.at("dZyp").assign(ptcls.num_particles, 0.0);

    my_timer.tic ("g2p");
    ptcls.g2p (vars, {"dZdx","dZdy"}, {"dZxp","dZyp"});
    my_timer.toc ("g2p");
	
    for (idx_t ip = 0; ip < num_particles; ++ip) {
      ptcls.dprops["hpZ"][ip] = ptcls.dprops["hp"][ip] + ptcls.dprops["Zp"][ip];
    }

    for (idx_t ip = 0; ip<num_particles; ++ip)  {
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
    
    while (t < data.T) 
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
        cel = std::abs(max_vel);
	
        if (it > 0)
          dt = 0.05 * data.hx / (1e-2 + cel); 
        std::cout << "time = " << t << "  " << " dt = " <<  dt << std::endl;
	std::cout << "cel = " << cel << std::endl;
        my_timer.toc ("update dt");

	my_timer.tic ("save csv");
        std::string filename = "nc_particles_";
        filename = filename + std::to_string (it++);
        filename = filename + ".csv";
        if (it % 25 == 0)
        {
          std::ofstream OF (filename.c_str ());
          ptcls.print<particles_t::output_format::csv>(OF);
          OF.close ();
        }
	my_timer.toc ("save csv");

        //  (0)  CONNECTIVITY and BASIS FUNCTIONS
        my_timer.tic ("reorder");
        ptcls.init_particle_mesh ();

        // ordering.resize (ptcls.num_particles);
        // idx_t iordering = 0;
        // for (auto const & ii : ptcls.grd_to_ptcl) {
        //   for (auto const & jj : ii.second) {
        //     ordering[iordering++] = jj;
        //   }
        // }

        // ptcls.reorder (ordering);
	// iordering = 0;
        // for (auto & ii : ptcls.grd_to_ptcl) {
        //   for (auto & jj : ii.second) {
        //     jj = iordering++;
        //   }
        // }
	// ptcls.init_particle_mesh ();
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

	my_timer.tic ("p2g");	
        ptcls.p2g (vars, {"Mp","mom_px","mom_py"}, {"Mv","mom_vx","mom_vy"});
	my_timer.toc ("p2g");

	my_timer.tic ("g2pd");
        ptcls.g2pd (vars, {"Z"}, {"dZxp"}, {"dZyp"});
        my_timer.toc ("g2pd");
	

        // (2)  EXTERNAL FORCES ON VERTICES (P2G)
        my_timer.tic ("step 2a");
        dpl::transform (policy, ptcls.dprops["Mp"].begin (), ptcls.dprops["Mp"].end (), ptcls.dprops["dZxp"].begin (),ptcls.dprops["Fpx"].begin (),
			[&] (double mp, double gzx) { return - data.g * mp * gzx; });
        dpl::transform (policy, ptcls.dprops["Mp"].begin (), ptcls.dprops["Mp"].end (), ptcls.dprops["dZyp"].begin (),ptcls.dprops["Fpy"].begin (),
			[&] (double mp, double gzy) { return - data.g * mp * gzy; });

        dpl::transform (policy, ptcls.dprops["Ap"].begin (), ptcls.dprops["Ap"].end (), ptcls.dprops["Fb_x"].begin (),
			ptcls.dprops["Fric_px"].begin (), std::multiplies<double> ());
        dpl::transform (policy, ptcls.dprops["Ap"].begin (), ptcls.dprops["Ap"].end (), ptcls.dprops["Fb_y"].begin (),
			ptcls.dprops["Fric_py"].begin (), std::multiplies<double> ());
        my_timer.toc ("step 2a");

	my_timer.tic ("p2g");
        ptcls.p2g (vars, {"Fpx","Fpy","Fric_px","Fric_py"}, {"FPxv","FPyv","Fric_x","Fric_y"});
	my_timer.toc ("p2g");
	
	my_timer.tic ("step 2b");
	dpl::transform (policy, vars["FPxv"].begin (), vars["FPxv"].end (), vars["Fric_x"].begin (),
			vars["F_ext_vx"].begin (), std::plus<double> ());
	dpl::transform (policy, vars["FPyv"].begin (), vars["FPyv"].end (), vars["Fric_y"].begin (),
			vars["F_ext_vy"].begin (), std::plus<double> ());

        my_timer.toc ("step 2b");

        // (3) INTERNAL FORCES (p2gd) and MOMENTUM BALANCE
        my_timer.tic ("p2gd");
        ptcls.p2gd (vars, {"F_11","F_21"}, {"F_12","F_22"}, "Vp", {"F_int_vx","F_int_vy"});
        my_timer.toc ("p2gd");
	
        my_timer.tic ("step 3");

	dpl::transform (policy, vars["F_ext_vx"].begin (), vars["F_ext_vx"].end (), vars["F_int_vx"].begin (), vars["Ftot_vx"].begin (), std::minus<double> ());
	dpl::transform (policy, vars["F_ext_vy"].begin (), vars["F_ext_vy"].end (), vars["F_int_vy"].begin (), vars["Ftot_vy"].begin (), std::minus<double> ());

	dpl::transform (policy, vars["Ftot_vx"].begin (), vars["Ftot_vx"].end (), vars["Ftot_vx"].begin (), vars["mom_vx"].begin (), [=] (double x, double y) { return dt*x + y; });
	dpl::transform (policy, vars["Ftot_vy"].begin (), vars["Ftot_vy"].end (), vars["Ftot_vy"].begin (), vars["mom_vy"].begin (), [=] (double x, double y) { return dt*x + y; });

        my_timer.toc ("step 3");

        // (4)  COMPUTE NODAL ACCELERATIONS AND VELOCITIES
        my_timer.tic ("step 4");

	dpl::transform (policy, vars["Ftot_vx"].begin (), vars["Ftot_vx"].end (), vars["Mv"].begin (), vars["avx"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	dpl::transform (policy, vars["Ftot_vy"].begin (), vars["Ftot_vy"].end (), vars["Mv"].begin (), vars["avy"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	dpl::transform (policy, vars["mom_vx"].begin (), vars["mom_vx"].end (), vars["Mv"].begin (), vars["vvx"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	dpl::transform (policy, vars["mom_vy"].begin (), vars["mom_vy"].end (), vars["Mv"].begin (), vars["vvy"].begin (), [] (double x, double y) { return y > 1.e-6 ? x/y : 0.0; });
	
        my_timer.toc ("step 4");

        // (5) BOUNDARY CONDITIONS - TO DO
        my_timer.tic ("step 5");
        /* We assume the slide will never reach the boundary and do nothing here! */
        my_timer.toc ("step 5");

        // (6) RETURN TO POINTS (G2P) and UPDATE POS AND VEL ON PARTICLES
        my_timer.tic ("step 6");
        ptcls.dprops.at("vpx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("vpy").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("apy").assign(ptcls.num_particles, 0.0);
	my_timer.toc ("step 6");

	my_timer.tic ("g2p");
        ptcls.g2p (vars, {"vvx","vvy","avx","avy"}, {"vpx","vpy","apx","apy"});
	my_timer.toc ("g2p");
	// if (it > 0) {
	//   for (auto const & iiii : vars.at ("avx")) {
	//     std::cout << iiii << std::endl;
	//   }

	//   assert (false);
	// }
	my_timer.tic ("step 6b");
	
	dpl::transform (policy, ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (),  ptcls.dprops["apx"].begin (), ptcls.dprops["vpx"].begin (), [=] (double x, double y) { return x + dt * y; } );
	dpl::transform (policy, ptcls.dprops["vpy"].begin (), ptcls.dprops["vpy"].end (),  ptcls.dprops["apy"].begin (), ptcls.dprops["vpy"].begin (), [=] (double x, double y) { return x + dt * y; } );
       
	dpl::transform (policy, ptcls.x.begin (), ptcls.x.end (),  ptcls.dprops["vpx"].begin (), ptcls.x.begin (), [=] (double x, double y) { return x + dt * y; } );
	dpl::transform (policy, ptcls.y.begin (), ptcls.y.end (),  ptcls.dprops["vpy"].begin (), ptcls.y.begin (), [=] (double x, double y) { return x + dt * y; } );

        my_timer.toc ("step 6b");

        // (7) COMPUTE HEIGHT WITH STRAIN (divergence of velocities)
        my_timer.tic ("step 7a");
        ptcls.dprops.at("vpx_dx").assign(ptcls.num_particles, 0.0);
        ptcls.dprops.at("vpy_dy").assign(ptcls.num_particles, 0.0);
	my_timer.toc ("step 7a");

	my_timer.tic ("g2pd");        
        ptcls.g2pd (vars, {"vvx","vvy"}, {"vpx_dx","vpy_dx"}, {"vpx_dy","vpy_dy"});
	my_timer.toc ("g2pd");
 
        my_timer.tic ("step 7");
	dpl::transform (policy, ptcls.dprops["vpx_dx"].begin (), ptcls.dprops["vpx_dx"].end (), ptcls.dprops["vpy_dy"].begin (), norm_v.begin (), std::plus<double> ());
	dpl::transform (policy, ptcls.dprops["hp"].begin (), ptcls.dprops["hp"].end (), norm_v.begin (), ptcls.dprops["hp"].begin (), [=] (double x, double y) { return x / (1 + dt * y); } );
	dpl::transform (policy, ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (), ptcls.dprops["Mp"].begin (), ptcls.dprops["mom_px"].begin (), std::multiplies<double> () );
	dpl::transform (policy, ptcls.dprops["vpy"].begin (), ptcls.dprops["vpy"].end (), ptcls.dprops["Mp"].begin (), ptcls.dprops["mom_py"].begin (), std::multiplies<double> () );
	dpl::transform (policy, ptcls.dprops["Vp"].begin (), ptcls.dprops["Vp"].end (), norm_v.begin (), ptcls.dprops["Vp"].begin (), [=] (double x, double y) { return x / (1 + dt * y); } );
	dpl::transform (policy, ptcls.dprops["Vp"].begin (), ptcls.dprops["Vp"].end (), ptcls.dprops["hp"].begin (), ptcls.dprops["Ap"].begin (), std::divides<double> () );
        my_timer.toc ("step 7");

	// (7b) UPDATE FRICTIONS
	my_timer.tic ("step 7b");
	dpl::transform (policy, ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (), ptcls.dprops["vpy"].begin (), norm_v.begin (), [] (double x, double y) { return std::sqrt (x*x + y*y); });
	dpl::transform (policy, ptcls.dprops["vpx"].begin (), ptcls.dprops["vpx"].end (), ptcls.dprops["hp"].begin (), ptcls.dprops["Fb_x"].begin (),
			[&] (double x, double y) { return - data.rho * data.g * ( y * std::tan(fric_ang) + x * x  / data.xi) * x; });
	dpl::transform (policy, ptcls.dprops["Fb_x"].begin (), ptcls.dprops["Fb_x"].end (), norm_v.begin (), ptcls.dprops["Fb_x"].begin (), std::divides<double> ());
	dpl::transform (policy, ptcls.dprops["vpy"].begin (), ptcls.dprops["vpy"].end (), ptcls.dprops["hp"].begin (), ptcls.dprops["Fb_y"].begin (),
			[&] (double x, double y) { return - data.rho * data.g * ( y * std::tan(fric_ang) + x * x  / data.xi) * x; });
	dpl::transform (policy, ptcls.dprops["Fb_y"].begin (), ptcls.dprops["Fb_y"].end (), norm_v.begin (), ptcls.dprops["Fb_y"].begin (), std::divides<double> ());
	my_timer.toc ("step 7b");


	// (8) UPDATE PARTICLE STRESS (USL)
        my_timer.tic ("step 8");
	
	range<idx_t> rng(0, ptcls.num_particles);
	stress_tensor_t st(ptcls, data);
	dpl::for_each (policy, rng.begin (), rng.end (), st);
	
	ptcls.dprops.at("hpZ").assign(ptcls.num_particles, 0.0);
	ptcls.dprops.at("Zp").assign(ptcls.num_particles, 0.0);
	my_timer.toc ("step 8");

	my_timer.tic ("g2p");
	ptcls.g2p (vars, {"Z"}, {"Zp"});
	my_timer.toc ("g2p");

	my_timer.tic ("step 8b");
	dpl::transform (policy, ptcls.dprops["hp"].begin (), ptcls.dprops["hp"].end (), ptcls.dprops["Zp"].begin (), ptcls.dprops["hpZ"].begin (), std::plus<double> ());
        my_timer.toc ("step 8b");

	my_timer.tic ("p2g");
        ptcls.p2g (vars,std::vector<std::string>{"hp"}, std::vector<std::string>{"HV"});
	my_timer.toc ("p2g");


	// my_timer.tic ("save vts");
	// filename = "nc_grid_";
	// filename = filename + std::to_string (it);
	// filename = filename + ".vts";
	// if (it % 250 == 0) {
        //  grid.vtk_export(filename.c_str(), Plotvars);
        //  grid.vtk_export(filename.c_str(), vars);
	// }
        t +=dt;
	// my_timer.toc ("save vts");

        //my_timer.print_report ();
	if (it % 3000 == 0) break;
	
      }

    my_timer.print_report ();
    //  ptcls.print<particles_t::output_format::csv>(std::cout);

    return 0;
}
