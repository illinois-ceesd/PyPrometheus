#include "chemGen.h"

void chemGen::writeMech(std::ostream& out) {

  out.setf(std::ios::scientific);
  out.precision(6);

  writeDefs(out);
  writeThermo(out);
  writeKinetics(out);
  
}

void chemGen::writeDefs(std::ostream& out) {

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

void chemGen::writeKinetics(std::ostream& out) {

  writeFalloffKinetics(out);
  writeRateCoefficients(out);
  writeNetRatesOfProgress(out);
  writeNetProductionRates(out);

}

void chemGen::writeThermo(std::ostream& out) {

  writeSpeciesSpecificHeats(out);
  writeSpeciesEnthalpies(out);
  writeSpeciesEntropies(out);
  writeSpeciesGibbs(out);
  writeEquilibriumConstants(out);
  writeNewtonTemperature(out);
  
  
}

void chemGen::writeFalloffKinetics(std::ostream& out) {

  std::vector<int>                          falloffIndices;
  std::vector<int>                          falloffType;
  std::vector<std::vector<double> >         falloffEfficiencies;
  std::vector<std::vector<double> >         falloffParams;
  std::shared_ptr<Cantera::Reaction>        reaction;
  std::shared_ptr<Cantera::FalloffReaction> work;

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
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

void chemGen::writeRateCoefficients(std::ostream& out) {

  std::vector<int>                   tbrIndices;
  std::vector<std::vector<double> >  tbrEfficiencies;
  std::shared_ptr<Cantera::Reaction> reaction;

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
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

void chemGen::writeArrheniusKinetics(int& rxn, std::shared_ptr<Cantera::Reaction>& reaction,
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

void chemGen::writeNetRatesOfProgress(std::ostream& out) {

  std::vector<int> reacIndices;
  std::vector<int> prodIndices;
  std::shared_ptr<Cantera::Reaction> reaction;

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
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

void chemGen::writeNetProductionRates(std::ostream& out) {

  std::vector<int> reacIndices;
  std::vector<int> prodIndices;
  std::shared_ptr<Cantera::Reaction> reaction;

  for(int k = 0; k < m_gas->nSpecies(); ++k) {
    m_rhs.insert( std::pair<int, std::vector<std::string> >(k, std::vector<std::string>()) );
  }

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
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

void chemGen::writeSpeciesSpecificHeats(std::ostream& out) {

  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(34);

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getSpecificHeats_R(" << m_baseType << "& T, std::vector<"
      << m_baseType << ">& cp0_R) {" << std::endl;
  out << std::endl;
  out << "    " << m_baseType << " tt0 = T;"       << std::endl;
  out << "    " << m_baseType << " tt1 = T * tt0;" << std::endl;
  out << "    " << m_baseType << " tt2 = T * tt1;" << std::endl;
  out << "    " << m_baseType << " tt3 = T * tt2;" << std::endl;
  out << "    " << m_baseType << " tt4 = 1.0 / T;"             << std::endl;
  out << "    " << m_baseType << " tt5 = tt4 / T;"             << std::endl;
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

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

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
      out << "    } else if(tt0 > " << c9[12] << " && tt0 < " << c9[13] <<  ") {" << std::endl;
      out << "      " << fmt("cp0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[16], "", out);
      writeCoeffTimesVar(true,  c9[17], " * tt0",  out);
      writeCoeffTimesVar(true,  c9[18], " * tt1", out);
      writeCoeffTimesVar(true,  c9[19], " * tt2", out);
      writeCoeffTimesVar(true,  c9[20], " * tt3", out);
      writeCoeffTimesVar(true,  c9[14], " * tt5", out);
      writeCoeffTimesVar(true,  c9[15], " * tt4;", out);
      out << std::endl;
      out << "    } else if(tt0 > " << c9[23] << " && tt0 < " << c9[24] << ") {" << std::endl;
      out << "      " << fmt("cp0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[27], "", out);
      writeCoeffTimesVar(true,  c9[28], " * tt0",  out);
      writeCoeffTimesVar(true,  c9[29], " * tt1", out);
      writeCoeffTimesVar(true,  c9[30], " * tt2", out);
      writeCoeffTimesVar(true,  c9[31], " * tt3", out);
      writeCoeffTimesVar(true,  c9[25], " * tt5", out);
      writeCoeffTimesVar(true,  c9[26], " * tt4;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;
      
    }

  }

  out << "  };" << std::endl;
  out << std::endl;

}

void chemGen::writeSpeciesEnthalpies(std::ostream& out) {

  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(34);

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getEnthalpies_RT(" << m_baseType << "& T, std::vector<"
      << m_baseType << ">& h0_RT) {" << std::endl;
  out << std::endl;
  out << "    " << m_baseType << " tt0 = T;"                   << std::endl;
  out << "    " << m_baseType << " tt1 = T * tt0;"             << std::endl;
  out << "    " << m_baseType << " tt2 = T * tt1;"             << std::endl;
  out << "    " << m_baseType << " tt3 = T * tt2;"             << std::endl;
  out << "    " << m_baseType << " tt4 = 1.0 / T;"             << std::endl;
  out << "    " << m_baseType << " tt5 = tt4 / T;"             << std::endl;
  out << "    " << m_baseType << " tt6 = std::log(tt0) * tt4;" << std::endl;
  out << std::endl;

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

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3]  = (-1) * c9[3];
      c9[14] = (-1) * c9[14];
      c9[25] = (-1) * c9[25];
      
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
      out << "    } else if(tt0 > " << c9[12] << " && tt0 < " << c9[13] <<  ") {" << std::endl;
      out << "      " << fmt("h0_RT", k) << " = ";
      writeCoeffTimesVar(false, c9[16], "", out);
      writeCoeffTimesVar(true,  c9[17], " * 0.50 * tt0",  out);
      writeCoeffTimesVar(true,  c9[18], " * OneThird * tt1", out);
      writeCoeffTimesVar(true,  c9[19], " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c9[20], " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c9[21], " * tt4", out);
      writeCoeffTimesVar(true,  c9[14], " * tt5", out);
      writeCoeffTimesVar(true,  c9[15], " * tt6;", out);
      out << std::endl;
      out << "    } else if(tt0 > " << c9[23] << " && tt0 < " << c9[24] << ") {" << std::endl;
      out << "      " << fmt("h0_RT", k) << " = ";
      writeCoeffTimesVar(false, c9[27], "", out);
      writeCoeffTimesVar(true,  c9[28], " * 0.50 * tt0",  out);
      writeCoeffTimesVar(true,  c9[29], " * OneThird * tt1", out);
      writeCoeffTimesVar(true,  c9[30], " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c9[31], " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c9[32], " * tt4", out);
      writeCoeffTimesVar(true,  c9[25], " * tt5", out);
      writeCoeffTimesVar(true,  c9[26], " * tt6;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    }

  }

  out << "  };" << std::endl;
  out << std::endl;
  
}

void chemGen::writeSpeciesEntropies(std::ostream& out) {

  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(34);

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getEntropies_R(" << m_baseType << "& T, std::vector<"
      << m_baseType << ">& s0_R) {" << std::endl;
  out << std::endl;
  out << "    " << m_baseType << " tt0 = T;"           << std::endl;
  out << "    " << m_baseType << " tt1 = T * tt0;"     << std::endl;
  out << "    " << m_baseType << " tt2 = T * tt1;"     << std::endl;
  out << "    " << m_baseType << " tt3 = T * tt2;"     << std::endl;
  out << "    " << m_baseType << " tt4 = 1.0 / T;"     << std::endl;
  out << "    " << m_baseType << " tt5 = tt4 / T;"     << std::endl; 
  out << "    " << m_baseType << " tt6 = std::log(T);" << std::endl;
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

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3]  = (-1) * c9[3];
      c9[4]  = (-1) * c9[4];
      c9[14] = (-1) * c9[14];
      c9[15] = (-1) * c9[15];
      c9[25] = (-1) * c9[25];
      c9[26] = (-1) * c9[26];

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
      out << "    } else if(tt0 > " << c9[12] << " && tt0 < " << c9[13] <<  ") {" << std::endl;
      out << "      " << fmt("s0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[22], "", out);
      writeCoeffTimesVar(true,  c9[16], " * tt6",  out);
      writeCoeffTimesVar(true,  c9[17], " * tt0", out);
      writeCoeffTimesVar(true,  c9[18], " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c9[19], " * OneThird * tt2", out);
      writeCoeffTimesVar(true,  c9[20], " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c9[14], " * 0.50 * tt5", out);
      writeCoeffTimesVar(true,  c9[15], " * tt4;", out);
      out << std::endl;
      out << "    } else if(tt0 > " << c9[23] << " && tt0 < " << c9[24] << ") {" << std::endl;
      out << "      " << fmt("s0_R", k) << " = ";
      writeCoeffTimesVar(false, c9[33], "", out);
      writeCoeffTimesVar(true,  c9[27], " * tt6",  out);
      writeCoeffTimesVar(true,  c9[28], " * tt0", out);
      writeCoeffTimesVar(true,  c9[29], " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c9[30], " * OneThird * tt2", out);
      writeCoeffTimesVar(true,  c9[31], " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c9[25], " * 0.50 * tt5", out);
      writeCoeffTimesVar(true,  c9[26], " * tt4;", out);
      out << std::endl;
      out << "    };" << std::endl;
      out << std::endl;

    }

  }

  out << "  };" << std::endl;
  out << std::endl;
  
}

void chemGen::writeSpeciesGibbs(std::ostream& out) {

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; } 
  out << "  void getGibbsFunctions_RT(" << m_baseType << "& T, std::vector<"
      << m_baseType << ">& g0_RT) {" << std::endl;
  out << std::endl;
  out << "    std::vector<" << m_baseType << "> h0_RT(kk, 0.0);" << std::endl;
  out << "    std::vector<" << m_baseType << "> s0_R(kk, 0.0);"  << std::endl;
  out << std::endl;
  out << "    getEnthalpies_RT(T, h0_RT);" << std::endl;
  out << "    getEntropies_R(T, s0_R);"    << std::endl;
  out << "    for(int k = 0; k < kk; ++k) { g0_RT[k] = h0_RT[k] - s0_R[k]; }" << std::endl;
  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;
  
}


void chemGen::writeEquilibriumConstants(std::ostream& out) {

  std::shared_ptr<Cantera::Reaction> reaction;
  std::vector<int>         reacIndices;
  std::vector<int>         prodIndices;

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getEquilibriumConstants(" << m_baseType << "& T, std::vector<"
      << m_baseType << ">& keq) {" << std::endl;
  out << std::endl;
  out << "    double            p0 = OneAtm;"          << std::endl;
  out << "    " << m_baseType << "              RT = GasConstant * T;" << std::endl;
  out << "    " << m_baseType << "              C0 = p0 / RT;"         << std::endl;
  out << "    std::vector<" << m_baseType << "> g0_RT(kk, 0.0);"       << std::endl;
  out << std::endl;
  out << "    getGibbsFunctions_RT(T, g0_RT);"    << std::endl;
  out << "    for(int k = 0; k < kk; ++k) { g0_RT[k] = exp(g0_RT[k]); }" << std::endl;
  out << std::endl;
  
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    reaction = m_gas->reaction(i);
    
    if(reaction->reversible == true) {
      
      getParticipatingSpecies(reaction, reacIndices, prodIndices);
      
      int dn = prodIndices.size() - reacIndices.size();
      std::string revSpecies = multiply("g0_RT", prodIndices); 
      std::string fwdSpecies = multiply("g0_RT", reacIndices);
 
      out << "    " << fmt("keq",i) << " = ";
      if(dn < 0) { out << " C0 * "; }
      out << "( " << revSpecies << " ) / "
	  << "( " << fwdSpecies;
      if(dn > 0) { out << " * C0 "; }
      out << " );" << std::endl;
      
    }

  }

  out << std::endl;
  out << "  };" << std::endl;
  out << std::endl;

}

void chemGen::writeNewtonTemperature(std::ostream& out) {

  if( m_temp == true ) { out << "  template <class Type>" << std::endl; }
  out << "  void getTemperature(double& h, double& Told, std::vector<"
      << m_baseType << ">& y, " << m_baseType << "& T) {" << std::endl;
  out << std::endl;
  out << "    double            tol   = 1.0e-06;" << std::endl;
  out << "    int               niter = 500;"     << std::endl;
  out << "    " << m_baseType << "              RT;"              << std::endl;
  out << "    " << m_baseType << "              To;"              << std::endl;
  out << "    " << m_baseType << "              Tp;"              << std::endl;
  out << "    " << m_baseType << "              dT = 1.0;"        << std::endl;
  out << "    " << m_baseType << "              fy = h;"          << std::endl;
  out << "    " << m_baseType << "              jy = 0.0;"        << std::endl;
  out << "    std::vector<" << m_baseType << "> hk(kk, 0.0);"     << std::endl;
  out << "    std::vector<" << m_baseType << "> cpk(kk, 0.0);"    << std::endl;
  out << std::endl;
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
  out << "	return;" << std::endl;
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

void chemGen::getParticipatingSpecies(std::shared_ptr<Cantera::Reaction>& reaction,
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
    if(iter->second == 2.0) { reacIndices.push_back( m_gas->speciesIndex(iter->first) ); }
  }

  for(Cantera::Composition::const_iterator iter = productsMap.begin();
      iter != productsMap.end();
      ++iter) {
    prodIndices.push_back( m_gas->speciesIndex(iter->first) );
    if(iter->second == 2.0) { prodIndices.push_back( m_gas->speciesIndex(iter->first) ); }
  }
  
}

void chemGen::getThirdBodyEfficiencies(std::shared_ptr<Cantera::Reaction>& reaction,
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

int chemGen::nFalloff() {

  int nfall = 0;
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    if(m_gas->reactionType(i) == Cantera::FALLOFF_RXN) { nfall += 1; }
    
  }

  return(nfall);
  
};
