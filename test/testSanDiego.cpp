#include <fstream>
#include "sanDiego.h"
#include "cantera/IdealGasMix.h"

void testProductionRates(Cantera::IdealGasMix& gas) {

  int     kk = gas.nSpecies();
  double  p  = Cantera::OneAtm;
  double  T0 = 300;
  double  dT = 100;
  double* x       = new double[kk];
  double* wdot_ct = new double[kk];
  std::vector<double> y(kk);
  std::vector<double> mw(kk);
  std::vector<double> wdot(kk,0.0);
  std::ofstream out;

  out.setf(std::ios::scientific);
  out.precision(4);

  /* get properties from cantera */
  out.open("outputData/sanDiego/speciesSourceTerm.ct.dat");

  for(int k = 0; k < kk; ++k) { x[k] = 1.0 / kk; }
  
  for(int m = 0; m < 20; ++m) {

    double T = T0 + 2 * m * dT;
    gas.setState_TPX(T, p, x);
    gas.getNetProductionRates(wdot_ct);
    gas.getMolecularWeights(&mw[0]);
    double rho = gas.density();

    out << T << "\t";
    for(int k = 0; k < kk; ++k) { out << wdot_ct[k] << "\t"; }
    out << std::endl;

  }
  out.close();

  /* now get properties from source.cpp */
  out.open("outputData/sanDiego/speciesSourceTerm.mech.dat");
  
  for(int m = 0; m < 40; ++m) {

    double T = T0 + m * dT;

    /* use cantera to get species densities */
    gas.setState_TPX(T, p, x);
    gas.getMassFractions(&y[0]);

    /* now call mech */
    mech::getNetProductionRates(p, T, y, wdot);

    out << T << "\t";
    for(int k = 0; k < kk; ++k) { out << wdot[k] << "\t"; }
    out << std::endl;

  }
  out.close();

}

void testNetRatesOfProgress(Cantera::IdealGasMix& gas) {

  int     ii = gas.nReactions();
  int     kk = gas.nSpecies();
  double  p  = Cantera::OneAtm;
  double  T0 = 300;
  double  dT = 100;
  double* x       = new double[kk];
  double* Rnet_ct = new double[ii];
  std::vector<double> y(kk);
  std::vector<double> c(kk);
  std::vector<double> mw(kk);
  std::vector<double> Rnet(ii);
  std::ofstream out;

  out.setf(std::ios::scientific);
  out.precision(4);

  /* set composition */
  for(int k = 0; k < kk; ++k) { x[k] = 1.0 / kk; }
  
  /* get properties from cantera */
  out.open("outputData/sanDiego/netRates.ct.dat");
  
  for(int m = 0; m < 20; ++m) {

    double T = T0 + 2 * m * dT;
    gas.setState_TPX(T, p, x);
    gas.getNetRatesOfProgress(Rnet_ct);

    out << T << "\t";
    for(int i = 0; i < ii; ++i) { out << Rnet_ct[i] << "\t"; }
    out << std::endl;

  }
  out.close();

  /* now get properties from source.cpp */
  out.open("outputData/sanDiego/netRates.mech.dat");
  
  for(int m = 0; m < 40; ++m) {

    double T = T0 + m * dT;

    /* use cantera to get species concentrations */
    gas.setState_TPX(T, p, x);
    gas.getMassFractions(&y[0]);
    gas.getMolecularWeights(&mw[0]);
    double rho = gas.density();
    for(int k = 0; k < kk; ++k) { c[k] = rho * y[k] / mw[k]; }

    /* now call mech */
    mech::getNetRatesOfProgress(T, c, Rnet);

    out << T << "\t";
    for(int i = 0; i < ii; ++i) { out << Rnet[i] << "\t"; }
    out << std::endl;

  }
  out.close();
  
}

void testRateCoefficients(Cantera::IdealGasMix& gas) {

  int     kk = gas.nSpecies();
  int     nr = gas.nReactions();
  double  p  = Cantera::OneAtm;
  double  T0 = 300;
  double  dT = 100;
  double* x       = new double[kk];
  double* kfwd_ct = new double[nr];
  double* krev_ct = new double[nr];
  std::vector<double> y(kk);
  std::vector<double> c(kk);
  std::vector<double> mw(kk);
  std::vector<double> kfwd(nr);
  std::vector<double> krev(nr);
  std::ofstream outFwd;
  std::ofstream outRev;

  outFwd.setf(std::ios::scientific);
  outFwd.precision(4);

  outRev.setf(std::ios::scientific);
  outRev.precision(4);

  /* set composition */
  for(int k = 0; k < kk; ++k) { x[k] = 1.0 / kk; }

  /* get fwd properties from cantera */
  outFwd.open("outputData/sanDiego/fwdRateCoeffs.ct.dat");
  outRev.open("outputData/sanDiego/revRateCoeffs.ct.dat");
  
  for(int m = 0; m < 20; ++m) {

    double T = T0 + 2 * m * dT;
    gas.setState_TPX(T, p, x);
    gas.getFwdRateConstants(kfwd_ct);
    gas.getRevRateConstants(krev_ct);

    outFwd << T << "\t";
    outRev << T << "\t";
    for(int i = 0; i < nr; ++i) {
      outFwd << kfwd_ct[i] << "\t";
      outRev << krev_ct[i] << "\t";
    }
    outFwd << std::endl;
    outRev << std::endl;

  }
  outFwd.close();
  outRev.close();
    
  /* now get fwd properties from source.cpp */
  outFwd.open("outputData/sanDiego/fwdRateCoeffs.mech.dat");
  outRev.open("outputData/sanDiego/revRateCoeffs.mech.dat");

  for(int m = 0; m < 40; ++m) {

    double T = T0 + m * dT;

    /* use cantera to get species concentrations */
    gas.setState_TPX(T, p, x);
    gas.getMassFractions(&y[0]);
    gas.getMolecularWeights(&mw[0]);
    double rho = gas.density();
    for(int k = 0; k < kk; ++k) { c[k] = rho * y[k] / mw[k]; }
    
    mech::getRateCoefficients(T, c, kfwd, krev);

    outFwd << T << "\t";
    outRev << T << "\t";
    for(int i = 0; i < nr; ++i) {
      outFwd << kfwd[i] << "\t";
      outRev << krev[i] << "\t";
    }
    outFwd << std::endl;
    outRev << std::endl;

  }
  
  outFwd.close();
  outRev.close();

}

void testEquilibriumConstants(Cantera::IdealGasMix& gas) {
  
  int     nr = gas.nReactions();
  double  p  = Cantera::OneAtm;
  double  T0 = 300;
  double  dT = 100;
  double* keq_ct = new double[nr];
  std::vector<double> keq(nr);
  std::ofstream out;

  out.setf(std::ios::scientific);
  out.precision(4);

  /* get properties from cantera */
  out.open("outputData/sanDiego/equilConstants.ct.dat");
  
  out << "T ";
  for(int i = 0; i < nr; ++i) { out << "keq"+std::to_string(i+1)+" "; }
  out << std::endl;
  
  for(int m = 0; m < 20; ++m) {

    double T = T0 + 2 * m * dT;
    gas.setState_TP(T, p);
    gas.getEquilibriumConstants(keq_ct);

    out << T << "\t";
    for(int i = 0; i < nr; ++i) { out << std::max(1.0e-18, keq_ct[i]) << "\t"; }
    out << std::endl;

  }
  out.close();

  /* now get properties from source.cpp */
  out.open("outputData/sanDiego/equilConstants.mech.dat");
  
  out << "T ";
  for(int i = 0; i < nr; ++i) { out << "keq"+std::to_string(i+1)+" "; }
  out << std::endl;

  for(int m = 0; m < 40; ++m) {

    double T = T0 + m * dT;
    mech::getEquilibriumConstants(T, keq);

    out << T << "\t";
    for(int i = 0; i < nr; ++i) { out << std::max(1.0e-18, 1.0 / keq[i]) << "\t"; }
    out << std::endl;
    
  }
  out.close();
  
}

void testThermo(Cantera::IdealGasMix& gas) {

  std::cout.setf(std::ios::scientific);
  std::cout.precision(4);
  
  int    ns  = gas.nSpecies();
  double p   = Cantera::OneAtm;
  double T0  = 300;
  double dT  = 100;
  double* cp_ct = new double[ns];
  double* h_ct  = new double[ns];
  double* s_ct  = new double[ns];
  std::ofstream outc;
  std::ofstream outh;
  std::ofstream outs;

  outc.setf(std::ios::scientific);
  outc.precision(4);

  outh.setf(std::ios::scientific);
  outh.precision(4);

  outs.setf(std::ios::scientific);
  outs.precision(4);
  
  /* get properties from cantera */
  outc.open("outputData/sanDiego/thermo.specificHeat.ct.dat");
  outh.open("outputData/sanDiego/thermo.enthalpy.ct.dat");
  outs.open("outputData/sanDiego/thermo.entropy.ct.dat");
  
  outc << "T ";
  for(int k = 0; k < ns; ++k) { outc << "cp"+std::to_string(k+1)+" "; }
  outc << std::endl;
  
  outh << "T\t";
  for(int k = 0; k < ns; ++k) { outh << "h"+std::to_string(k+1)+"\t"; }
  outh << std::endl;
  
  outs << "T\t";
  for(int k = 0; k < ns; ++k) { outs << "s"+std::to_string(k+1)+"\t"; }
  outs << std::endl;
  
  for(int m = 0; m < 20; ++m) {

    double T = T0 + 2 * m * dT;
    gas.setState_TP(T, p);
    gas.getCp_R(cp_ct);
    gas.getEnthalpy_RT(h_ct);
    gas.getEntropy_R(s_ct);

    outc << T << "\t";
    outh << T << "\t";
    outs << T << "\t";
    for(int k = 0; k < ns; ++k) { outc << cp_ct[k] << "\t"; }
    for(int k = 0; k < ns; ++k) { outh << h_ct[k]  << "\t"; }
    for(int k = 0; k < ns; ++k) { outs << s_ct[k]  << "\t"; }
    outc << std::endl;
    outh << std::endl;
    outs << std::endl;

  }
  outc.close();
  outh.close();
  outs.close();

  /* now get properties from source.cpp */
  outc.open("outputData/sanDiego/thermo.specificHeat.mech.dat");
  outh.open("outputData/sanDiego/thermo.enthalpy.mech.dat");
  outs.open("outputData/sanDiego/thermo.entropy.mech.dat");
  
  for(int m = 0; m < 40; ++m) {

    std::vector<double> cp(ns);
    std::vector<double> h(ns);
    std::vector<double> s(ns);

    double T = T0 + m * dT;
    mech::getSpecificHeats_R(T, cp);
    mech::getEnthalpies_RT(T, h);
    mech::getEntropies_R(T, s);

    outc << T << "\t";
    outh << T << "\t";
    outs << T << "\t";
    for(int k = 0; k < ns; ++k) { outc << cp[k] << "\t"; }
    for(int k = 0; k < ns; ++k) { outh << h[k]  << "\t"; }
    for(int k = 0; k < ns; ++k) { outs << s[k]  << "\t"; }
    outc << std::endl;
    outh << std::endl;
    outs << std::endl;

  }
  outc.close();
  outh.close();
  outs.close();
  
}

int main() {

  Cantera::IdealGasMix gas("ctis/sanDiego.cti","gas");

  testThermo(gas);
  testEquilibriumConstants(gas);
  testRateCoefficients(gas);
  testNetRatesOfProgress(gas);
  testProductionRates(gas);
  
}
