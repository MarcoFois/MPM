#ifndef HAVE_MPM_DATA
#define HAVE_MPM_DATA

#include <json.hpp>

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
  double xi;
  std::vector<double> Z;
  std::vector<double> dZdx;
  std::vector<double> dZdy;

  DATA (const char* filename);

};

void
from_json (const nlohmann::json &j, DATA &d)
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
  j.at("xi").get_to(d.xi);
  j.at("Z").get_to(d.Z);
  j.at("dZdx").get_to(d.dZdx);
  j.at("dZdy").get_to(d.dZdy);
}

DATA::DATA (const char* filename) {
  nlohmann::json json;
  std::ifstream json_file (filename);
  json_file>>json;
  json_file.close ();
  json.get_to (*this);
};

#endif
