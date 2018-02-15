#include <vector>

namespace mech {
  
  /* thermodynamics */
  void getSpecificHeats_R(double& T, std::vector<double>& cp0_R);
  void getEnthalpies_RT(double& T, std::vector<double>& h0_RT);
  void getEntropies_R(double& T, std::vector<double>& s0_R);
  void getGibbsFunctions_RT(double& T, std::vector<double>& g0_RT);
  void getEquilibriumConstants(double& T, std::vector<double>& keq);
  void getTemperature(double& h, double& Told, std::vector<double>& y, double& T);
  
  /* rates coefficients */
  void getFalloffRates(double& T, std::vector<double>& C, std::vector<double>& kfwd);
  void getRateCoefficients(double& T, std::vector<double>& C,
			   std::vector<double>& kfwd, std::vector<double>& krev);
  
  /* net production rates */
  void getNetRatesOfProgress(double& T, std::vector<double>& C,
			     std::vector<double>& Rnet);
  void getNetProductionRates(double& p, double& T, std::vector<double>& y,
			     std::vector<double>& wdot);

}
