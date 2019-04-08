#include "Prometheus.h"

void Prometheus::greetings() {

  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  std::cout << "Welcome to Prometheus..." << std::endl;
  std::cout << "Working with " << m_mech << " mech..." << std::endl;
  std::cout << "Number of species:   " << m_gas->nSpecies()   << std::endl;
  std::cout << "Number of elements:  " << m_gas->nElements()  << std::endl;
  std::cout << "Number of reactions: " << m_gas->nReactions() << std::endl;
  std::cout << "The implementation will ";
  if( m_ooriented == false ) { std::cout << "not "; }
  std::cout << "be object-oriented..." << std::endl;
  std::cout << "The implementation will ";
  if( m_templated == false ) { std::cout << "not "; }
  std::cout << "be templated..." << std::endl;
  
}

void Prometheus::writeMech() {

  std::string   filename;
  std::ofstream out;
  out.setf(std::ios::scientific);
  out.precision(6);

  greetings();

  /* thermo */
  if( m_ooriented == true ) {

    filename = m_mech+"Thermo.h";
    out.open( filename.c_str() );
    writeThermoHeader( out );
    
  } else {

    filename = m_mech+".h";
    out.open( filename.c_str() );
    writeDefs( out );
    
  }

  writeThermo(out);

  /* kinetics */
  if( m_ooriented == true ) {

    out.close();
    filename = m_mech+"Kinetics.h";
    out.open( filename.c_str() );
    writeKineticsHeader( out );
  }
  
  writeKinetics( out );
  
}

void Prometheus::writeDefs(std::ostream& out) {

  double* mw = new double[m_gas->nSpecies()];
  m_gas->getMolecularWeights(mw);

  out << "#include <cmath>"    << std::endl;
  out << "#include <vector>"   << std::endl;
  out << "#include <iostream>" << std::endl;
  out << std::endl;
  out << "namespace mech {" << std::endl;
  out << std::endl;
  out << "  int    mm = " << m_gas->nElements()  << ";" << std::endl;
  out << "  int    kk = " << m_gas->nSpecies()   << ";" << std::endl;
  out << "  int    ii = " << m_gas->nReactions() << ";" << std::endl;
  out << "  int    ff = " << nFalloff()          << ";" << std::endl;
  out << "  double OneAtm      = 1.01325e5;"     << std::endl;
  out << "  double OneThird    = 1.0 / 3.0;"     << std::endl;
  out << "  double GasConstant = 8314.4621;"     << std::endl;
  out << "  double BigNumber   = 1.0e300;"       << std::endl; 
  out << std::endl;
  out << "  std::vector<double> mw = { ";
  for(int k = 0; k < m_gas->nSpecies(); ++k) {
    if(k != m_gas->nSpecies()-1) {
      out << mw[k] << ", ";
    } else {
      out << mw[k];
    }
  }
  out << " };" << std::endl;
  out << std::endl;
  
}

void Prometheus::writeKineticsHeader(std::ostream& out) {

  int ii = m_gas->nReactions();
  int ff = nFalloff();

  out << "#ifndef MECH_THERMO" << std::endl;
  out << "#define MECH_THERMO" << std::endl;

  out << "include <cmath>" << std::endl;
  out << "include <vector>" << std::endl;
  out << "include \"gasThermo.h\"" << std::endl;
  out  << std::endl;

  out << "namespace mech {" << std::endl;
  out << std::endl;
  out << "  template <class ModelType>" << std::endl;
  out << "    class GasKinetics {" << std::endl;
  out << std::endl;
  out << "  public:" << std::endl;
  out << std::endl;
  out << "  GasKinetics() : "
      << "m_ii(" << ii << "), "
      << "m_ff(" << ff << ") {" << std::endl;
  out << std::endl;

  for(int i = 0; i < ii; ++i) {
    out << "      m_reactionString[" << i << "] = \"";
    out << m_gas->reactionString(i) << "\";"
	<< std::endl;    
  }
  
  out << std::endl;
  out << "    }" << std::endl;
  out << std::endl;
  
  out << "    int nReactions() { return m_ii; }"  << std::endl;
  out << std::endl;
  out << "    int nFalloff() { return m_ff; }" << std::endl;
  out << std::endl;
  
  
}

void Prometheus::writeThermoHeader(std::ostream& out) {

  int mm = m_gas->nElements();
  int kk = m_gas->nSpecies();
  std::vector<double> mw( kk, 0.0 );
  m_gas->getMolecularWeights( &mw[0] );
  

  out << "#ifndef MECH_THERMO" << std::endl;
  out << "#define MECH_THERMO" << std::endl;
  out << std::endl;
  
  out << "#include <cmath>" << std::endl;
  out << "#include <vector>" << std::endl;
  out << "#include <Eigen/Dense>" << std::endl;
  out << "#include <cppad/cppad.hpp>" << std::endl;
  out << "#include mechDefs.h" << std::endl;
  out << std::endl;

  out << "namespace mech {" << std::endl;
  out << std::endl;
  out << "  template <class ModelType>" << std::endl;
  out << "    class GasThermo {" << std::endl;
  out << std::endl;
  out << "  public:" << std::endl;
  out << std::endl;
  out << "  GasThermo() : "
      << "m_kk(" << kk << "), "
      << "m_mm(" << mm << ") {"
      << std::endl;
  out << std::endl;

  out << "      m_mw.resize( m_kk );" << std::endl;
  out << "      m_y0.resize( m_kk );" << std::endl;
  out << "      m_Emat.setZero( m_kk, m_mm );" << std::endl;
  out << std::endl;
  for(int k = 0; k < kk; ++k) {
    out << "      m_speciesIndex[\""
	<< m_gas->speciesName(k) << "\"] = "
	<< k << ";"
	<< std::endl;
  }
  out << std::endl;
  for(int k = 0; k < kk; ++k) {
    out << "      m_speciesName[" << k << "] = \""
	<< m_gas->speciesName(k) << "\";"
	<< std::endl;
  }
  out << std::endl;
  for(int k = 0; k < kk; ++k) {
    out << "      m_mw[" << k << "] = " << mw[k] << ";" << std::endl;
  }
  out << std::endl;
  for(int k = 0; k < kk; ++k) {
    for(int m = 0; m < mm; ++m) {
      out << "      m_Emat(" << k << "," << m << ") = "
	  << m_gas->nAtoms( k, m ) << ";"
	  << std::endl;
    }
    out << std::endl;
  }
  out << std::endl;
  out << "    }" << std::endl;
  out << std::endl;

  out << "    int nSpecies() { return m_kk; }" << std::endl;
  out << std::endl;
  out << "    int nElements() { return m_mm; }" << std::endl;
  out << std::endl;
  out << "    int nVariables() { return m_kk; }" << std::endl;
  out << std::endl;
  out << "    int nConserved() { return m_mm; }" << std::endl;
  out << std::endl;
  out << "    ModelType temperature() { return m_T; }" << std::endl;
  out << std::endl;
  out << "    ModelType pressure() { return m_p; }" << std::endl;
  out << std::endl;
  out << "    ModelType enthalpy() { return m_h; }" << std::endl;
  out << std::endl;
  out << "    std::vector<double> molecularWeights() { return m_mw; }" << std::endl;
  out << std::endl;
  out << "    std::vector<ModelType> initialMassFractions() { return m_y0; }" << std::endl;
  out << std::endl;
  out << "    void setTemperature(ModelType& T);" << std::endl;
  out << std::endl;
  out << "    void setPressrure(ModelType& p);" << std::endl;
  out << std::endl;
  out << "    void setEnthalpyMass(ModelType& T0, std::vector<ModelType>& x0);" << std::endl;
  out << std::endl;
  out << "    void moleToMassFractions(std::vector<ModelType>& x, std::vector<ModelType>& y);"
      << std::endl;
  out << std::endl;
  out << "    void setInitialState_TPPhi(std::string& fuel,"    << std::endl;
  out << "                               std::string& diluent," << std::endl;
  out << "                               ModelType& T0,"        << std::endl;
  out << "                               ModelType& p,"         << std::endl;
  out << "                               ModelType& phi);"      << std::endl;
  out << std::endl;
  out << "    void setInitialState_TPX(ModelType& T0,"           << std::endl;
  out << "                             ModelType& p,"            << std::endl;
  out << "                             std::vector<double>& x);" << std::endl;
  out << std::endl;
  out << "    template<class DataType>" << std::endl;
  out << "      void getDensity(DataType& T," << std::endl;
  out << "                      std::vector<DataType>& y," << std::endl;
  out << "                      DataType& rho);" << std::endl;
  out << std::endl;
  out << "    template<class DataType>" << std::endl;
  out << "      void getTemperature(ModelType& Told," << std::endl;
  out << "                          std::vector<DataType>& y," << std::endl;
  out << "                          DataType& T);" << std::endl;
  out << std::endl;
  out << "    template<class DataType>" << std::endl;
  out << "      void getSpecificHeats_R(DataType& T," << std::endl;
  out << "                              std::vector<DataType>& cp0_R);" << std::endl;
  out << std::endl;
  out << "    template<class DataType>" << std::endl;
  out << "      void getEnthalpies_RT(DataType& T," << std::endl;
  out << "                            std::vector<DataType>& h0_RT);" << std::endl;
  out << std::endl;
  out << "    template<class DataType>" << std::endl;
  out << "      void getEntropies_R(DataType& T," << std::endl;
  out << "                          std::vector<DataType>& s0_R);" << std::endl;
  out << std::endl;
  out << "    template<class DataType>" << std::endl;
  out << "      void getGibbsFunctions_RT(DataType& T," << std::endl;
  out << "                                std::vector<DataType>& g0_RT);" << std::endl;
  out << std::endl;
  out << "    template<class DataType>" << std::endl;
  out << "      void getEquilibriumConstants(DataType& T," << std::endl;
  out << "                                   std::vector<DataType>& keq);" << std::endl;
  out << std::endl;
  
  out << "  protected:" << std::endl;
  out << "   int m_kk;" << std::endl;
  out << "   int m_mm;" << std::endl;
  out << "   ModelType                 m_T;" << std::endl;
  out << "   ModelType                 m_p;" << std::endl;
  out << "   ModelType                 m_h;" << std::endl;
  out << "   std::vector<double>       m_mw;" << std::endl;
  out << "   std::vector<ModelType>    m_y0;" << std::endl;
  out << "   Eigen::MatrixXd           m_Emat;" << std::endl;
  out << "   std::map<std::string,int> m_speciesIndex;" << std::endl;
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;
  
  // out << "  template<class ModelType>" << std::endl;
  // out << "    void GasThermo<ModelType>::setTemperature(ModelType& T) { m_T = T; };"
  //     << std::endl;
  // out << std::endl;
  
  // out << "  template<class ModelType>" << std::endl;
  // out << "    void GasThermo<ModelType>::setPressure(ModelType& p) { m_p = p; };"
  //     << std::endl;
  
  // out << "  template<class ModelType>" << std::endl;
  // out << "    void GasThermo<ModelType>::setEnthalpyMass(ModelType& T0," << std::endl;
  // out << "					       std::vector<ModelType>& x0) {"
  //     << std::endl;
  // out << std::endl;  
  // out << "    ModelType              RT = GasConstant * T0;" << std::endl;
  // out << "    std::vector<ModelType> hk( m_kk, 0.0 );" << std::endl;
  // out << std::endl;
  // out << "    m_h = 0.0;" << std::endl;
  // out << std::endl;
  // out << "    moleToMassFractions( x0, m_y0 );" << std::endl;
  // out << std::endl;
  // out << "    getEnthalpies_RT( T0, hk );" << std::endl;
  // out << "    for(int k = 0; k < m_kk; ++k) {" << std::endl;
  // out << "      hk[k] = RT * hk[k];" << std::endl;
  // out << "      m_h  += hk[k] * m_y0[k] / m_mw[k];" << std::endl;
  // out << "    }" << std::endl;
  // out << std::endl;
  // out << "  };" << std::endl;
  // out << std::endl;

  // out << "  template<class ModelType>" << std::endl;
  // out << "    void GasThermo<ModelType>::moleToMassFractions(std::vector<ModelType>& x,"
  //     << std::endl;
  // out << "						   std::vector<ModelType>& y) {"
  //     << std::endl;
  // out << std::endl;
  // out << "   ModelType xk;  " << std::endl;
  // out << "   ModelType mmw = 0.0;" << std::endl;
  // out << std::endl;
  // out << "    for(int k = 0; k < m_kk; ++k) {" << std::endl;
  // out << "      xk = std::max( x[k], 0.0 );" << std::endl;
  // out << "      y[k]  = xk;" << std::endl;
  // out << "      mmw  += m_mw[k] * xk;" << std::endl;
  // out << "    }" << std::endl;
  // out << std::endl;
  // out << "  };" << std::endl;
  // out << std::endl;
  // out << "  template<class ModelType>" << std::endl;
  // out << "    void GasThermo<ModelType>::setInitialState_TPPhi(ModelType& T0," << std::endl;
  // out << "						     ModelType& p," << std::endl;
  // out << "						     ModelType& phi) {" << std::endl;
  // out << std::endl;
  
}

void Prometheus::writeKinetics(std::ostream& out) {

  std::cout << "\t ... Writing falloff kinetics ... " << std::endl;
  writeFalloffKinetics(out);
  std::cout << "\t ... Writing rate coefficients ... " << std::endl;
  writeRateCoefficients(out);
  std::cout << "\t ... Writing reation rates ... " << std::endl;
  writeNetRatesOfProgress(out);
  std::cout << "\t ... Writing species production rates ... " << std::endl;
  writeNetProductionRates(out);

}

void Prometheus::writeThermo(std::ostream& out) {

  std::cout << "\t ... Writing specific heats ... " << std::endl;
  writeSpeciesSpecificHeats(out);

  std::cout << "\t ... Writing enthalpies ... " << std::endl;
  writeSpeciesEnthalpies(out);

  std::cout << "\t ... Writing entropies ... " << std::endl;
  writeSpeciesEntropies(out);

  std::cout << "\t ... Writing gibbs functions ... " << std::endl;
  writeSpeciesGibbs(out);

  std::cout << "\t ... Writing equilibrium constants ... " << std::endl;
  writeEquilibriumConstants(out);

  std::cout << "\t ... Writing Newton method for temperature inversion ... " << std::endl;
  writeNewtonTemperature(out);
  
  
}

void Prometheus::writeFalloffKinetics(std::ostream& out) {

  std::vector<int>                          falloffIndices;
  std::vector<int>                          falloffType;
  std::vector<std::vector<double> >         falloffEfficiencies;
  std::vector<std::vector<double> >         falloffParams;
  std::shared_ptr<Cantera::Reaction>        reaction;
  std::shared_ptr<Cantera::FalloffReaction> work;

  if( m_templated == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getFalloffRates(" << m_baseType << "& T, std::vector<" << m_baseType << ">& C,";
  out << " std::vector<" << m_baseType << ">& kfwd) {" << std::endl;
  out << std::endl;
  out << "    int TROE = 110;" << std::endl;
  out << "    " << m_baseType << " lpr;" << std::endl;
  out << "    " << m_baseType << " cc;"  << std::endl;
  out << "    " << m_baseType << " nn;"  << std::endl;
  out << "    " << m_baseType << " f1;"  << std::endl;
  out << "    " << m_baseType << " logT = std::log(T);" << std::endl;
  out << "    " << m_baseType << " invT = 1.0/T;"       << std::endl;
  out << "    std::vector<" << m_baseType << "> khi(ff,0.0);"  << std::endl;
  out << "    std::vector<" << m_baseType << "> klo(ff,0.0);"  << std::endl;
  out << "    std::vector<" << m_baseType << "> pr(ff,0.0);"   << std::endl;
  out << "    std::vector<" << m_baseType << "> work(ff,0.0);" << std::endl;
  out << "    std::vector<int> falloffType(ff,100);" << std::endl;
  out << std::endl;

  /* identify falloff reactions */
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    reaction = m_gas->reaction(i);
    if(reaction->reaction_type == Cantera::FALLOFF_RXN) { falloffIndices.push_back(i); }

  }

  /* Parameterization */
  for(int i = 0; i < falloffIndices.size(); ++i) {

    reaction = m_gas->reaction( falloffIndices[i] );
    work     = std::dynamic_pointer_cast<Cantera::FalloffReaction>(reaction);

    /* Arrhenius expressions */
    double Af  = work->high_rate.preExponentialFactor();
    double bf  = work->high_rate.temperatureExponent();
    double Taf = (-1) * work->high_rate.activationEnergy_R();

    double A0  = work->low_rate.preExponentialFactor();
    double b0  = work->low_rate.temperatureExponent();
    double Ta0 = (-1) * work->low_rate.activationEnergy_R();

    out << "    " << fmt("khi",i) << " = std::exp(" << std::log(Af);
    if(bf  != 0.0) { writeCoeffTimesVar(true, bf, " * logT", out); }
    if(Taf != 0.0) { writeCoeffTimesVar(true, Taf," * invT", out); }
    out << ");" << std::endl;

    out << "    " << fmt("klo",i) << " = std::exp(" << std::log(A0);
    if(b0  != 0.0) { writeCoeffTimesVar(true, b0, " * logT", out); }
    if(Ta0 != 0.0) { writeCoeffTimesVar(true, Ta0," * invT", out); }
    out << ");" << std::endl;
    out << std::endl;

    /* Third-body efficiencies */
    std::vector<double> efficiencies;
    for(int k = 0; k < m_gas->nSpecies(); ++k) {
      const std::string name = m_gas->speciesName(k);
      double eff = work->third_body.efficiency(name);
      efficiencies.push_back(eff);
    }
    falloffEfficiencies.push_back(efficiencies);

    /* falloff function parameters */
    std::vector<double> params(4);
    work->falloff->getParameters(&params[0]);
    falloffParams.push_back(params);

    /* store falloff function type */
    falloffType.push_back( work->falloff->getType() );
    
  }

  /* Third-body efficiencies */
  out.precision(2);
  for(int i = 0; i < falloffIndices.size(); ++i) {

    out << "    " << fmt("pr",i) << " = ";

    const std::string Ck   = fmt("C",0);
    double eff = falloffEfficiencies[i][0];
    if(eff != 0.0) { writeCoeffTimesVar(false, eff, " * "+Ck, out); }
    
    for(int k = 1; k < m_gas->nSpecies(); ++k) {
      const std::string Ck   = fmt("C",k);
      double eff = falloffEfficiencies[i][k];
      if(eff != 0.0) { out << " + "; writeCoeffTimesVar(false, eff, " * "+Ck, out); }
    }

    out << "; " << std::endl;

  }
  out << std::endl;
  out << "    for(int i = 0; i < ff; ++i) { pr[i] *= (klo[i] / khi[i]); }" << std::endl;
  out << std::endl;

  /* falloff function */
  for(int i = 0; i < falloffIndices.size(); ++i) {

    out << "    " << fmt("falloffType",i) << " = " << falloffType[i] << ";" << std::endl;;
    
  }
  out << std::endl;
  
  out.precision(6);
  for(int i = 0; i < falloffIndices.size(); ++i) {

    if( falloffType[i] == Cantera::SIMPLE_FALLOFF ) {

      out << "    " << fmt("work",i) << " = 1.0;" << std::endl;

    } else if( falloffType[i] == Cantera::TROE_FALLOFF ) {

      std::vector<double> params = falloffParams[i];
      params[1] = (-1) * 1.0 / params[1];
      params[2] = (-1) * 1.0 / params[2];
      params[3] = (-1) * params[3];
    
      out << "    " << fmt("work",i) << " = ";
      out << "(1.0 - " << params[0] << ") * std::exp(";
      writeCoeffTimesVar(true,params[1]," * T",out);
      out << ")";
      out << " + " << params[0] << " * std::exp(";
      writeCoeffTimesVar(true,params[2]," * T",out);
      out << ")";
      if(params[3] != 0.0) {
	out << " + std::exp(";
	writeCoeffTimesVar(true,params[3]," * invT",out);
	out<< ");" << std::endl;
      } else {
	out << ";" << std::endl;
      }

    }
    
  }
  out << std::endl;
  
  /* falloff calculator */
  out << "    for(int i = 0; i < ff; ++i) {" << std::endl;
  out << "      lpr =  std::log10(pr[i]);"                     << std::endl;
  out << "      if(falloffType[i] == TROE) {"                      << std::endl;
  out << "        cc      = -0.40 - 0.67 * std::log10(work[i]);"   << std::endl;
  out << "        nn      =  0.75 - 1.27 * std::log10(work[i]);"   << std::endl;
  out << "        f1      =  (lpr + cc)/(nn - 0.14 * (lpr + cc));" << std::endl;
  out << "        work[i] =  std::log10(work[i])/(1 + f1 * f1);"   << std::endl;
  out << "        work[i] =  std::pow(10.0, work[i]);"             << std::endl;
  out << "      }"                                                 << std::endl;
  out << "      work[i] =  (pr[i] * work[i])/(1 + pr[i]);"         << std::endl;
  out << "    }" << std::endl;
  out << std::endl;
  
  /* fwd rate coefficients */
  for(int i = 0; i < falloffIndices.size(); ++i) {
    int rxn = falloffIndices[i];
    out << "    "
	<< fmt("kfwd",rxn) << " = "
	<< fmt("khi",i)    << " * "
	<< fmt("work",i)   << ";"
	<< std::endl;
  }
  
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;
  

}

void Prometheus::writeRateCoefficients(std::ostream& out) {

  std::vector<int>                   tbrIndices;
  std::vector<std::vector<double> >  tbrEfficiencies;
  std::shared_ptr<Cantera::Reaction> reaction;

  if( m_templated == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getRateCoefficients(" << m_baseType
      << "& T, std::vector<" << m_baseType << ">& C,";
  out << " std::vector<" << m_baseType << ">& kfwd, std::vector<"
      << m_baseType << ">& krev) {" << std::endl;
  out << std::endl;
  out << "    " << m_baseType << " logT = std::log(T);"         << std::endl;
  out << "    " << m_baseType << " invT = 1.0/T;"               << std::endl;
  out << "    std::vector<" << m_baseType << "> keq(ii,0.0);"   << std::endl;
  out << std::endl;
  out << "    getEquilibriumConstants(T, keq);" << std::endl;
  out << "    for(int i = 0; i < ii; ++i) { keq[i] = exp( keq[i] ); }" << std::endl;
  out << "    for(int i = 0; i < ii; ++i) { if( keq[i] > BigNumber ) { keq[i] = BigNumber; } }" << std::endl;
  out << std::endl;

  /* Arrhenius kinetics */
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    reaction = m_gas->reaction(i);

    switch(reaction->reaction_type) {

    case Cantera::ELEMENTARY_RXN:

      writeArrheniusKinetics(i, reaction, out);
      break;
      
    case Cantera::THREE_BODY_RXN:

      tbrIndices.push_back(i);
      getThirdBodyEfficiencies(reaction, tbrEfficiencies);
      writeArrheniusKinetics(i, reaction, out);
      break;
      
    }

  }
  out << std::endl;

  /* Three-body reactions */
  out.precision(3);
  for(int i = 0; i < tbrIndices.size(); ++i) {

    out << "    " << fmt("kfwd",tbrIndices[i]) << " *= ( ";
    
    const std::string Ck = fmt("C",0);
    double eff = tbrEfficiencies[i][0];
    if(eff != 0.0) { writeCoeffTimesVar(false, eff, " * "+Ck, out); }
    
    for(int k = 1; k < m_gas->nSpecies(); ++k) {
      const std::string Ck = fmt("C",k);
      double eff = tbrEfficiencies[i][k];
      if(eff != 0.0) { writeCoeffTimesVar(true, eff, " * "+Ck, out); }
    }
    
    out << " ); " << std::endl;
    
  }
  out << std::endl;

  out << "    getFalloffRates(T, C, kfwd);" << std::endl;
  out << std::endl;
  out << "    for(int i = 0; i < ii; ++i) { krev[i] = kfwd[i] * keq[i]; }" << std::endl;
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;

}

void Prometheus::writeArrheniusKinetics(int& rxn, std::shared_ptr<Cantera::Reaction>& reaction,
				     std::ostream& out) {

  out.precision(6);
  
  std::shared_ptr<Cantera::ElementaryReaction> work;
  work = std::dynamic_pointer_cast<Cantera::ElementaryReaction>(reaction);

  double A  = work->rate.preExponentialFactor();
  double b  = work->rate.temperatureExponent();
  double Ta = (-1) * work->rate.activationEnergy_R();

  out << "    " << fmt("kfwd",rxn) << " = std::exp(" << std::log(A);
  if(b  != 0.0) { writeCoeffTimesVar(true, b, " * logT", out); }
  if(Ta != 0.0) { writeCoeffTimesVar(true, Ta," * invT", out); }
  out << ");" << std::endl;

}

void Prometheus::writeNetRatesOfProgress(std::ostream& out) {

  std::vector<int> reacIndices;
  std::vector<int> prodIndices;
  std::shared_ptr<Cantera::Reaction> reaction;

  if( m_templated == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getNetRatesOfProgress(" << m_baseType << "& T, std::vector<"
      << m_baseType << ">& C, std::vector<" << m_baseType << ">& Rnet) {";
  out << std::endl;
  out << std::endl;
  out << "    std::vector<" << m_baseType << "> kfwd(ii,0.0);"  << std::endl;
  out << "    std::vector<" << m_baseType << "> krev(ii,0.0);"  << std::endl;
  out << "    std::vector<" << m_baseType << "> Rfwd(ii,0.0);"  << std::endl;
  out << "    std::vector<" << m_baseType << "> Rrev(ii,0.0);"  << std::endl;
  out << std::endl;
  out << "    getRateCoefficients(T, C, kfwd, krev);"       << std::endl;
  out << std::endl;

  for(int i = 0; i < m_gas->nReactions(); ++i) {

    reaction = m_gas->reaction(i);
    getParticipatingSpecies(reaction, reacIndices, prodIndices);
    std::string fwdSpecies = multiply("C",reacIndices);
    std::string revSpecies = multiply("C",prodIndices);

    out << "    " << fmt("Rfwd",i)
	<< " = "  << fmt("kfwd",i)
	<< " * "  << fwdSpecies
	<< ";"
	<< std::endl;
    if(reaction->reversible == true) {
      out << "    " << fmt("Rrev",i)
	  << " = "  << fmt("krev",i)
	  << " * "  << revSpecies
	  << ";"
	  << std::endl;
    }
    out << std::endl;
    
  }

  out << "    for(int i = 0; i < ii; ++i) { Rnet[i] = Rfwd[i] - Rrev[i]; }" << std::endl;
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;
  
  
}

void Prometheus::writeNetProductionRates(std::ostream& out) {

  std::vector<int> reacIndices;
  std::vector<int> prodIndices;
  std::shared_ptr<Cantera::Reaction> reaction;

  for(int k = 0; k < m_gas->nSpecies(); ++k) {
    m_rhs.insert( std::pair<int, std::vector<std::string> >(k, std::vector<std::string>()) );
  }

  if( m_templated == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getNetProductionRates(double& p, " << m_baseType << "& T, std::vector<"
      << m_baseType << ">& y, std::vector<" << m_baseType << ">& omega) {";
  out << std::endl;
  out << std::endl;
  out << "    " << m_baseType << " W;"   << std::endl;
  out << "    " << m_baseType << " rho;" << std::endl;
  out << "    std::vector<" << m_baseType << "> C(kk,0.0);"     << std::endl;
  out << "    std::vector<" << m_baseType << "> Rnet(ii,0.0);"  << std::endl;
  out << std::endl;

  out << "    W   = 0.0;" << std::endl;
  out << "    for(int k = 0; k < kk; ++k) { W += y[k] / mw[k]; }" << std::endl;
  out << "    W   = 1.0 / W;" << std::endl;
  out << "    rho = (p * W) / (GasConstant * T);" << std::endl;
  out << "    for(int k = 0; k < kk; ++k) { C[k] = rho * y[k] / mw[k]; }" << std::endl;
  out << std::endl;
  out << "    getNetRatesOfProgress(T, C, Rnet);"       << std::endl;
  out << std::endl;

  for(int i = 0; i < m_gas->nReactions(); ++i) {

    reaction = m_gas->reaction(i);
    getParticipatingSpecies(reaction, reacIndices, prodIndices);

    for(int k = 0; k < reacIndices.size(); ++k) {
      std::string rhs = " - " + fmt("Rnet",i);
      m_rhs[reacIndices[k]].push_back(rhs);
    }

    for(int k = 0; k < prodIndices.size(); ++k) {
      std::string rhs = " + " + fmt("Rnet",i); 
      m_rhs[prodIndices[k]].push_back(rhs);
    }
    
  }

  for(std::map<int, std::vector<std::string> >::iterator iter = m_rhs.begin();
      iter != m_rhs.end();
      ++iter)
    {
      std::string rhs = wrapString(iter->second);
      if(rhs != "") {
	out << "    " << fmt("omega",iter->first)
	    << " = "  << rhs
	    << ";"    << std::endl;
      }
    }
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;
  out << "}" << std::endl;
  
  
}

void Prometheus::writeSpeciesSpecificHeats(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);
  //std::vector<double> c9(34);

  /* write function name */
  writeFunctionName( out, "GasThermo", "getSpecificHeats_R" );

  /* arguments */
  out << "(";
  out << m_dataType << "& T, std::vector<"
      << m_dataType << ">& cp0_R";
  out << ") {" << std::endl;

  /* declarations */
  out << std::endl;
  out << "    " << m_dataType << " tt0 = T;"       << std::endl;
  out << "    " << m_dataType << " tt1 = T * tt0;" << std::endl;
  out << "    " << m_dataType << " tt2 = T * tt1;" << std::endl;
  out << "    " << m_dataType << " tt3 = T * tt2;" << std::endl;
  out << "    " << m_dataType << " tt4 = 1.0 / T;" << std::endl;
  out << "    " << m_dataType << " tt5 = tt4 / T;" << std::endl;
  out << std::endl;

  for(int k = 0; k < m_gas->nSpecies(); ++k) {

    type = m_gas->thermo().speciesThermo().reportType(k);

    if(type == NASA) { 

      m_gas->thermo().speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);
      
      out << "    if(tt0 > " << c[0] << ") {" << std::endl;
      out << "      " << fmt("cp0_R", k) << " = ";
      writeCoeffTimesVar(false, c[1], "", out);
      writeCoeffTimesVar(true,  c[2], " * tt0", out);
      writeCoeffTimesVar(true,  c[3], " * tt1", out);
      writeCoeffTimesVar(true,  c[4], " * tt2", out);
      writeCoeffTimesVar(true,  c[5], " * tt3;", out);
      out << std::endl;
      out << "    } else {" << std::endl;
      out << "      " << fmt("cp0_R", k) << " = ";
      writeCoeffTimesVar(false, c[8], "", out);
      writeCoeffTimesVar(true,  c[9],  " * tt0",  out);
      writeCoeffTimesVar(true,  c[10], " * tt1", out);
      writeCoeffTimesVar(true,  c[11], " * tt2", out);
      writeCoeffTimesVar(true,  c[12], " * tt3;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    } else if(type == NASA9) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      /* write out single zone */
      out << "    if(tt0 > " << c9[1] << " && tt0 < " << c9[2] << ") {" << std::endl;
      out << "      " << fmt("cp0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[5],  "", out);
      writeCoeffTimesVar(true,  c9[6],  " * tt0", out);
      writeCoeffTimesVar(true,  c9[7],  " * tt1", out);
      writeCoeffTimesVar(true,  c9[8],  " * tt2", out);
      writeCoeffTimesVar(true,  c9[9],  " * tt3", out);
      writeCoeffTimesVar(true,  c9[3],  " * tt5", out);
      writeCoeffTimesVar(true,  c9[4],  " * tt4;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];

      /* write out zone 1 */
      out << "    if(tt0 > " << c9[1] << " && tt0 < " << c9[2] << ") {" << std::endl;
      out << "      " << fmt("cp0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[5],  "", out);
      writeCoeffTimesVar(true,  c9[6],  " * tt0", out);
      writeCoeffTimesVar(true,  c9[7],  " * tt1", out);
      writeCoeffTimesVar(true,  c9[8],  " * tt2", out);
      writeCoeffTimesVar(true,  c9[9],  " * tt3", out);
      writeCoeffTimesVar(true,  c9[3],  " * tt5", out);
      writeCoeffTimesVar(true,  c9[4],  " * tt4;", out);
      out << std::endl;
      
      /* write out other zones */
      for(int z = 1; z < nzones; ++z) { 

	int minTempIdx = 11 * z + 1;
	int maxTempIdx = 11 * z + 2;	

	out << "    } else if(tt0 > " << c9[ minTempIdx ]
	    << " && tt0 < " << c9[ maxTempIdx ]
	    <<  ") {" << std::endl;
	
	out << "      " << fmt("cp0_R", k) << " = ";
	writeCoeffTimesVar(false, c9[ maxTempIdx + 3 ], "", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 4 ], " * tt0",  out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 5 ], " * tt1", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 6 ], " * tt2", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 7 ], " * tt3", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 1 ], " * tt5", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 2 ], " * tt4;", out);
	out << std::endl;

      }

      out << "    };" << std::endl;
      out << std::endl;
      
    }

  }

  out << "  };" << std::endl;
  out << std::endl;

}

void Prometheus::writeSpeciesEnthalpies(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  writeFunctionName( out, "GasThermo", "getEnthalpies_RT" );

  /* arguments */
  out << "(";
  out << m_dataType << "& T, std::vector<"
      << m_dataType << ">& h0_RT";
  out << ") {" << std::endl;

  /* declarations */
  out << "    " << m_dataType << " tt0 = T;"                   << std::endl;
  out << "    " << m_dataType << " tt1 = T * tt0;"             << std::endl;
  out << "    " << m_dataType << " tt2 = T * tt1;"             << std::endl;
  out << "    " << m_dataType << " tt3 = T * tt2;"             << std::endl;
  out << "    " << m_dataType << " tt4 = 1.0 / T;"             << std::endl;
  out << "    " << m_dataType << " tt5 = tt4 / T;"             << std::endl;
  out << "    " << m_dataType << " tt6 = std::log(tt0) * tt4;" << std::endl;
  out << std::endl;

  /* write out the code */

  for(int k = 0; k < m_gas->nSpecies(); ++k) {

    type = m_gas->thermo().speciesThermo().reportType(k);
    
    if(type == NASA) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);
    
      out << "    if(tt0 > " << c[0] << ") {" << std::endl;
      out << "      " << fmt("h0_RT", k) << " = ";
      writeCoeffTimesVar(false, c[1], "", out);
      writeCoeffTimesVar(true,  c[2], " * 0.50 * tt0", out);
      writeCoeffTimesVar(true,  c[3], " * OneThird * tt1", out);
      writeCoeffTimesVar(true,  c[4], " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c[5], " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c[6], " * tt4;", out);
      out << std::endl;
      out << "    } else {" << std::endl;
      out << "      " << fmt("h0_RT", k) << " = ";
      writeCoeffTimesVar(false, c[8], "", out);
      writeCoeffTimesVar(true,  c[9], " * 0.50 * tt0",  out);
      writeCoeffTimesVar(true,  c[10], " * OneThird * tt1", out);
      writeCoeffTimesVar(true,  c[11], " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c[12], " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c[13], " * tt4;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    } else if(type == NASA9) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3] *= -1.0;

      /* write single zone */
      out << "    if(tt0 > " << c9[1] << " && tt0 < " << c9[2] << ") {" << std::endl;
      out << "      " << fmt("h0_RT", k) << " = ";
      writeCoeffTimesVar(false, c9[5],  "", out);
      writeCoeffTimesVar(true,  c9[6],  " * 0.50 * tt0", out);
      writeCoeffTimesVar(true,  c9[7],  " * OneThird * tt1", out);
      writeCoeffTimesVar(true,  c9[8],  " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c9[9],  " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c9[10], " * tt4", out);
      writeCoeffTimesVar(true,  c9[3],  " * tt5", out);
      writeCoeffTimesVar(true,  c9[4],  " * tt6;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];
      
      c9[3] *= -1.0;      

      /* write zone 1 */
      out << "    if(tt0 > " << c9[1] << " && tt0 < " << c9[2] << ") {" << std::endl;
      out << "      " << fmt("h0_RT", k) << " = ";
      writeCoeffTimesVar(false, c9[5],  "", out);
      writeCoeffTimesVar(true,  c9[6],  " * 0.50 * tt0", out);
      writeCoeffTimesVar(true,  c9[7],  " * OneThird * tt1", out);
      writeCoeffTimesVar(true,  c9[8],  " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c9[9],  " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c9[10], " * tt4", out);
      writeCoeffTimesVar(true,  c9[3],  " * tt5", out);
      writeCoeffTimesVar(true,  c9[4],  " * tt6;", out);
      out << std::endl;

      /* write out other zones */
      for(int z = 1; z < nzones; ++z) {

	int minTempIdx = 11 * z + 1;
	int maxTempIdx = 11 * z + 2;

	c9[ maxTempIdx + 1 ] *= -1.0;

	out << "    } else if(tt0 > " << c9[ minTempIdx ]
	    << " && tt0 < " << c9[ maxTempIdx ]
	    <<  ") {" << std::endl;
      
	out << "      " << fmt("h0_RT", k) << " = ";
	writeCoeffTimesVar(false, c9[ maxTempIdx + 3 ], "", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 4 ], " * 0.50 * tt0",  out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 5 ], " * OneThird * tt1", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 6 ], " * 0.25 * tt2", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 7 ], " * 0.20 * tt3", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 8 ], " * tt4", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 1 ], " * tt5", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 2 ], " * tt6;", out);
	out << std::endl;
	
      }      
      
      out << "    };" << std::endl;
      out << std::endl;

    }

  }

  out << "  };" << std::endl;
  out << std::endl;
  
}

void Prometheus::writeSpeciesEntropies(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  writeFunctionName( out, "GasThermo", "getEntropies_R" );

  /* arguments */
  out << "(";
  out << m_dataType << "& T, std::vector<"
      << m_dataType << ">& s0_R";
  out << ") {" << std::endl;

  /* declarations */
  out << "    " << m_dataType << " tt0 = T;"           << std::endl;
  out << "    " << m_dataType << " tt1 = T * tt0;"     << std::endl;
  out << "    " << m_dataType << " tt2 = T * tt1;"     << std::endl;
  out << "    " << m_dataType << " tt3 = T * tt2;"     << std::endl;
  out << "    " << m_dataType << " tt4 = 1.0 / T;"     << std::endl;
  out << "    " << m_dataType << " tt5 = tt4 / T;"     << std::endl; 
  out << "    " << m_dataType << " tt6 = std::log(T);" << std::endl;
  out << std::endl;

  for(int k = 0; k < m_gas->nSpecies(); ++k) {

    type = m_gas->thermo().speciesThermo().reportType(k);

    if(type == NASA) {
    
      m_gas->thermo().speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);

      out << "    if(tt0 > " << c[0] << ") {" << std::endl;
      out << "      " << fmt("s0_R", k) << " = ";
      writeCoeffTimesVar(false, c[1], " * tt6", out);
      writeCoeffTimesVar(true,  c[2], " * tt0", out);
      writeCoeffTimesVar(true,  c[3], " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c[4], " * OneThird * tt2", out);
      writeCoeffTimesVar(true,  c[5], " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c[7], ";", out);
      out << std::endl;
      out << "    } else {" << std::endl;
      out << "      " << fmt("s0_R", k) << " = ";
      writeCoeffTimesVar(false, c[8], " * tt6", out);
      writeCoeffTimesVar(true,  c[9], " * tt0",  out);
      writeCoeffTimesVar(true,  c[10], " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c[11], " * OneThird * tt2", out);
      writeCoeffTimesVar(true,  c[12], " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c[14], ";", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    } else if(type == NASA9) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3] *= -1.0;
      c9[4] *= -1.0;

      /* write single zone */      
      out << "    if(tt0 > " << c9[1] << " && tt0 < " << c9[2] << ") {" << std::endl;
      out << "      " << fmt("s0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[11], "", out);
      writeCoeffTimesVar(true,  c9[5],  " * tt6", out);
      writeCoeffTimesVar(true,  c9[6],  " * tt0", out);
      writeCoeffTimesVar(true,  c9[7],  " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c9[8],  " * OneThird * tt2", out);
      writeCoeffTimesVar(true,  c9[9],  " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c9[3],  " * 0.50 * tt5", out);
      writeCoeffTimesVar(true,  c9[4],  " * tt4;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];
      
      c9[3] *= -1.0;
      c9[4] *= -1.0;

      /* write out zone 1 */      
      out << "    if(tt0 > " << c9[1] << " && tt0 < " << c9[2] << ") {" << std::endl;
      out << "      " << fmt("s0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[11], "", out);
      writeCoeffTimesVar(true,  c9[5],  " * tt6", out);
      writeCoeffTimesVar(true,  c9[6],  " * tt0", out);
      writeCoeffTimesVar(true,  c9[7],  " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c9[8],  " * OneThird * tt2", out);
      writeCoeffTimesVar(true,  c9[9],  " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c9[3],  " * 0.50 * tt5", out);
      writeCoeffTimesVar(true,  c9[4],  " * tt4;", out);
      out << std::endl;

      /* write out other zones */
      for(int z = 1; z < nzones; ++z) {

	int minTempIdx = 11 * z + 1;
	int maxTempIdx = 11 * z + 2;

	c9[ maxTempIdx + 1 ] *= -1.0;
	c9[ maxTempIdx + 2 ] *= -1.0;
	
	out << "    } else if(tt0 > " << c9[ minTempIdx ]
	    << " && tt0 < " << c9[ maxTempIdx ]
	    <<  ") {" << std::endl;

	out << "      " << fmt("s0_R", k) << " = ";
	writeCoeffTimesVar(false, c9[ maxTempIdx + 9 ], "", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 3 ], " * tt6",  out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 4 ], " * tt0", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 5 ], " * 0.50 * tt1", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 6 ], " * OneThird * tt2", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 7 ], " * 0.25 * tt3", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 1 ], " * 0.50 * tt5", out);
	writeCoeffTimesVar(true,  c9[ maxTempIdx + 2 ], " * tt4;", out);
	out << std::endl;    

      }
      
      out << "    };" << std::endl;
      out << std::endl;

    }

  }

  out << "  };" << std::endl;
  out << std::endl;
  
}

void Prometheus::writeSpeciesGibbs(std::ostream& out) {

  /* write function name */
  writeFunctionName( out, "GasThermo", "getGibbsFunctions_RT" );

  /* arguments */
  out << "("
      << m_dataType << "& T, std::vector<"
      << m_dataType << ">& g0_RT"
      << ") {" << std::endl;
  out << std::endl;

  /* declarations */
  out << "    std::vector<" << m_dataType << "> h0_RT(kk, 0.0);" << std::endl;
  out << "    std::vector<" << m_dataType << "> s0_R(kk, 0.0);"  << std::endl;
  out << std::endl;

  /* write out the code */
  out << "    getEnthalpies_RT(T, h0_RT);" << std::endl;
  out << "    getEntropies_R(T, s0_R);"    << std::endl;
  out << "    for(int k = 0; k < kk; ++k) { g0_RT[k] = h0_RT[k] - s0_R[k]; }" << std::endl;
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;
  
}


void Prometheus::writeEquilibriumConstants(std::ostream& out) {

  std::shared_ptr<Cantera::Reaction> reaction;
  std::vector<int>         reacIndices;
  std::vector<int>         prodIndices;

  /* write function name */
  writeFunctionName( out, "GasThermo", "getEquilibriumConstants" );
  
  /* arguments*/
  out << "(";
  if( m_ooriented == true ) { out << m_baseType; } else { out << m_dataType; }
  out << "& T, std::vector<"
      << m_dataType << ">& keq";
  out << ") {" << std::endl;
  out << std::endl;

  /* declarations*/
  out << "    ";
  if( m_ooriented == true ) { out << m_baseType; } else { out << m_dataType; }
  out << " p0 = OneAtm;"          << std::endl;
  out << "    " << m_dataType << " RT = GasConstant * T;" << std::endl;
  out << "    " << m_dataType << " C0 = p0 / RT;"         << std::endl;
  out << "    std::vector<" << m_dataType << "> g0_RT(kk, 0.0);"       << std::endl;
  out << std::endl;

  /* write out the code */
  out << "    getGibbsFunctions_RT( T, g0_RT );"    << std::endl;
  //out << "    for(int k = 0; k < kk; ++k) { g0_RT[k] = exp( g0_RT[k] ); }" << std::endl;
  out << std::endl;
  
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    reaction = m_gas->reaction(i);
    
    if(reaction->reversible == true) {

      getParticipatingSpecies(reaction, reacIndices, prodIndices);
      
      int dn = prodIndices.size() - reacIndices.size();
      std::string revSpecies = increment("g0_RT", prodIndices); 
      std::string fwdSpecies = increment("g0_RT", reacIndices);
 
      out << "    " << fmt("keq",i) << " = ";
      if(dn < 0) { out << " C0 + "; }
      out << "( " << revSpecies << " ) - "
	  << "( " << fwdSpecies;
      if(dn > 0) { for(int n = 0; n < dn; ++n) { out << " + C0 "; } }
      out << " );" << std::endl;
      
    }

  }

  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;

}

void Prometheus::writeNewtonTemperature(std::ostream& out) {

  /* write function name */
  writeFunctionName( out, "GasThermo", "getTemperature" );
  
  /* arguments */
  out << "(";
  if( m_ooriented == true ) { out << m_baseType; } else { out << m_dataType; }
  out << "& h, ";
  if( m_ooriented == true ) { out << m_baseType; } else { out << m_dataType; }
  out <<"& Told, std::vector<"
      << m_dataType << ">& y, "
      << m_dataType << "& T";
  out << ") {" << std::endl;
  out << std::endl;
  
  /* declarations */
  out << "    double tol   = 1.0e-06;" << std::endl;
  out << "    int    niter = 500;"     << std::endl;
  out << "    " << m_dataType << " RT;"              << std::endl;
  out << "    " << m_dataType << " To;"              << std::endl;
  out << "    " << m_dataType << " Tp;"              << std::endl;
  out << "    " << m_dataType << " dT = 1.0;"        << std::endl;
  out << "    " << m_dataType << " fy = h;"          << std::endl;
  out << "    " << m_dataType << " jy = 0.0;"        << std::endl;
  out << "    std::vector<" << m_dataType << "> hk(kk, 0.0);"     << std::endl;
  out << "    std::vector<" << m_dataType << "> cpk(kk, 0.0);"    << std::endl;
  out << std::endl;

  /* write out the code */
  out << "    To = Told;" << std::endl;
  out << "    Tp = Told;" << std::endl;
  out << std::endl;
  out << "    for(int it = 0; it < niter; ++it) { " << std::endl;
  out << std::endl;
  out << "      RT = GasConstant * To;"       << std::endl;
  out << "      getSpecificHeats_R(To, cpk);" << std::endl;
  out << "      getEnthalpies_RT(To, hk);"    << std::endl;
  out << "      for(int k = 0; k < kk; ++k) { " << "hk[k]  = RT * hk[k] / mw[k]; }"
      << std::endl;
  out << "      for(int k = 0; k < kk; ++k) { " << "cpk[k] = GasConstant * cpk[k] / mw[k]; }"
      << std::endl;
  out << "      for(int k = 0; k < kk; ++k) { fy -=  hk[k] * y[k]; }" << std::endl;
  out << "      for(int k = 0; k < kk; ++k) { jy -= cpk[k] * y[k]; }" << std::endl;
  out << "      dT  = -fy / jy;" << std::endl;
  out << "      Tp  =  To + dT;" << std::endl;
  out << "      To  =  Tp;"      << std::endl;
  out << std::endl;
  out << "      if( (std::fabs(dT) < tol)) {" << std::endl;
  out << "	T = Tp;" << std::endl;
  out << "	break; " << std::endl;
  out << "      }" << std::endl;
  out << std::endl;
  out << "      fy = h;"   << std::endl;
  out << "      jy = 0.0;" << std::endl;
  out << std::endl;
  out << "    }" << std::endl;
  out << std::endl;
  out << "    T = Tp;" << std::endl;
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;

}

void Prometheus::getParticipatingSpecies(std::shared_ptr<Cantera::Reaction>& reaction,
				      std::vector<int>& reacIndices,
				      std::vector<int>& prodIndices) {

  reacIndices.clear();
  prodIndices.clear(); 
  
  Cantera::Composition reactantsMap = reaction->reactants;
  Cantera::Composition productsMap  = reaction->products;

  for(Cantera::Composition::const_iterator iter = reactantsMap.begin();
      iter != reactantsMap.end();
      ++iter) {
    reacIndices.push_back( m_gas->speciesIndex(iter->first) );
    if(iter->second > 1.0) {
      int nu = (int)iter->second;      
      for(int k = 0; k < nu-1; ++k) { reacIndices.push_back( m_gas->speciesIndex(iter->first) ); }
    }
  }

  // for(int k = 0; k < reacIndices.size(); ++k) { std::cout << m_gas->speciesName( reacIndices[k] ) << "\t"; }
  // std::cout << std::endl;
  
  for(Cantera::Composition::const_iterator iter = productsMap.begin();
      iter != productsMap.end();
      ++iter) {
    prodIndices.push_back( m_gas->speciesIndex(iter->first) );
    if(iter->second > 1.0) {
      int nu = (int)iter->second;     
      for(int k = 0; k < nu-1; ++k) { prodIndices.push_back( m_gas->speciesIndex(iter->first) ); }
    }
  }
  
}

void Prometheus::getThirdBodyEfficiencies(std::shared_ptr<Cantera::Reaction>& reaction,
				       std::vector<std::vector<double> >& tbrEfficiencies) {

  std::vector<double> efficiencies;
  std::shared_ptr<Cantera::ThreeBodyReaction> work;
  work = std::dynamic_pointer_cast<Cantera::ThreeBodyReaction>(reaction);
  
  for(int k = 0; k < m_gas->nSpecies(); ++k) {

    const std::string name = m_gas->speciesName(k);
    double eff = work->third_body.efficiency(name);
    efficiencies.push_back(eff);
    
  }

  tbrEfficiencies.push_back(efficiencies);

}

int Prometheus::nFalloff() {

  int nfall = 0;
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    if(m_gas->reactionType(i) == Cantera::FALLOFF_RXN) { nfall += 1; }
    
  }

  return(nfall);
  
};

void Prometheus::writeFunctionName(std::ostream& out,
				const std::string& className,
				const std::string& funName) {

  if( m_templated == true && m_ooriented == true ) {
    
    out << "  template <class " << m_baseType << ">" << std::endl;
    out << "  template <class " << m_dataType << ">" << std::endl;
    
  } else if( m_templated == true && m_ooriented == false ) {
    
    out << "  template <class " << m_dataType << ">" << std::endl;
    
  }

  if( m_ooriented == true ) {
    out << "  void " << className << "<" << m_baseType << ">::" << funName;
  } else {
    out << "  void " << funName;
  }
  
}
