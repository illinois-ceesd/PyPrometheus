#include "PrometheusPython.h"

void PyPrometheus::greetings() {

  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  std::cout << "Welcome to PyPrometheus..." << std::endl;
  std::cout << "Working with " << m_mech << " mech..." << std::endl;
  std::cout << "Number of species:   " << m_gas->nSpecies()   << std::endl;
  std::cout << "Number of elements:  " << m_gas->nElements()  << std::endl;
  std::cout << "Number of reactions: " << m_gas->nReactions() << std::endl;
  
}

void PyPrometheus::writeMech() {

  std::string   filename;
  std::ofstream out;
  out.setf(std::ios::scientific);
  out.precision(6);

  greetings();

  /* defs */
  filename = m_mech + ".py";
  out.open( filename.c_str() );
  writeDefs( out );
  
  /* thermo */
  writeThermo( out );

  /* kinetics */  
  //writeKinetics( out );
  
}

void PyPrometheus::writeDefs(std::ostream& out) {

  double* mw = new double[m_gas->nSpecies()];
  m_gas->getMolecularWeights(mw);

  out << "import numpy as np"    << std::endl;
  out << std::endl;
  out << "class prometheusChemistry():" << std::endl;
  out << std::endl;
  out << "    def __init__(self):" << std::endl;
  out << std::endl;
  out << "        self.mm = " << m_gas->nElements()  << ";" << std::endl;
  out << "        self.kk = " << m_gas->nSpecies()   << ";" << std::endl;
  out << "        self.ii = " << m_gas->nReactions() << ";" << std::endl;
  out << "        self.one_atm   = 1.01325e5;"     << std::endl;
  out << "        self.one_third = 1.0 / 3.0;"     << std::endl;
  out << "        self.gas_constant = 8314.4621;"     << std::endl;
  out << "        self.big_number   = 1.0e300;"       << std::endl;
  out << std::endl;
  out << "        self.wts = np.array( [ ";
  for(int k = 0; k < m_gas->nSpecies(); ++k) {
    if(k != m_gas->nSpecies()-1) {
      out << mw[k] << ", ";
    } else {
      out << mw[k];
    }
  }
  out << " ] )" << std::endl;
  out << std::endl;
  out << "        return" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::writeKinetics(std::ostream& out) {

  std::cout << "\t ... Writing rate coefficients ... " << std::endl;
  writeRateCoefficients(out);
  std::cout << "\t ... Writing reation rates ... " << std::endl;
  writeNetRatesOfProgress(out);
  std::cout << "\t ... Writing species production rates ... " << std::endl;
  writeNetProductionRates(out);

}

void PyPrometheus::writeThermo(std::ostream& out) {

  std::cout << "\t ... Writing mixture specific heat ... " << std::endl;
  writeMixtureSpecificHeat(out);
  
  std::cout << "\t ... Writing species specific heats ... " << std::endl;
  writeSpeciesSpecificHeats(out);

  std::cout << "\t ... Writing species enthalpies ... " << std::endl;
  writeSpeciesEnthalpies(out);

  std::cout << "\t ... Writing species entropies ... " << std::endl;
  writeSpeciesEntropies(out);

  std::cout << "\t ... Writing species gibbs functions ... " << std::endl;
  writeSpeciesGibbs(out);

  std::cout << "\t ... Writing equilibrium constants ... " << std::endl;
  writeEquilibriumConstants(out);
    
}

void PyPrometheus::writeRateCoefficients(std::ostream& out) {

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

void PyPrometheus::writeArrheniusKinetics(int& rxn, std::shared_ptr<Cantera::Reaction>& reaction,
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

void PyPrometheus::writeNetRatesOfProgress(std::ostream& out) {

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

void PyPrometheus::writeNetProductionRates(std::ostream& out) {

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

void PyPrometheus::writeMixtureSpecificHeat(std::ostream& out) {

  /* write function name */
  writeFunctionName( out,  "get_mixture_specific_heat_mass" );

  /* arguments */
  out << "(self,T,Y):" << std::endl;

  /* allocate buffers */
  out << std::endl;
  out << "        nz,ny,nx = T.shape" << std::endl;
  out << "        cp = np.zeros( [ nz,ny,nx ] )" << std::endl;
  out << std::endl;

  /* write out the code */
  out << "        " << "cp0_R = self.get_species_specific_heats_R( T )" << std::endl;
  for(int k = 0; k < m_gas->nSpecies(); ++k) {
    out << "        cp += "
	<< fmt( "Y", k )     << " * "
	<< fmt( "cp0_R", k ) << " / "
	<< fmt( "self.wts", k )
	<< std::endl;
  }
  out << "        cp *= self.gas_constant" << std::endl;
  out << std::endl;
  out << "        return cp" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::writeSpeciesSpecificHeats(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  writeFunctionName( out,  "get_species_specific_heats_R" );

  /* arguments */
  out << "(self,T):" << std::endl;

  /* declarations */
  out << std::endl;
  out << "        " << "tt0 = T"                   << std::endl;
  out << "        " << "tt1 = T * tt0"             << std::endl;
  out << "        " << "tt2 = T * tt1"             << std::endl;
  out << "        " << "tt3 = T * tt2"             << std::endl;
  out << "        " << "tt4 = np.power( T, -1.0 )" << std::endl;
  out << "        " << "tt5 = tt4 * tt4"           << std::endl;
  out << std::endl;

  /* allocate buffers */
  out << "        nz,ny,nx = T.shape" << std::endl;
  out << "        cp0_R = np.zeros( [ self.kk,nz,ny,nx ] )" << std::endl;
  out << std::endl;

  for(int k = 0; k < m_gas->nSpecies(); ++k) {

    type = m_gas->thermo().speciesThermo().reportType(k);

    if(type == NASA) { 

      m_gas->thermo().speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);

      out << "        cp_high  = ";      
      writeCoeffTimesVar(false, c[1], "", out);
      writeCoeffTimesVar(true,  c[2], " * tt0", out);
      writeCoeffTimesVar(true,  c[3], " * tt1", out);
      writeCoeffTimesVar(true,  c[4], " * tt2", out);
      writeCoeffTimesVar(true,  c[5], " * tt3", out);
      out << std::endl;

      out << "        cp_low   = ";
      writeCoeffTimesVar(false, c[8], "", out);
      writeCoeffTimesVar(true,  c[9],  " * tt0", out);
      writeCoeffTimesVar(true,  c[10], " * tt1", out);
      writeCoeffTimesVar(true,  c[11], " * tt2", out);
      writeCoeffTimesVar(true,  c[12], " * tt3", out);
      out << std::endl;

      out << "        " << fmt( "cp0_R", k )
	  << " = np.where( tt0 < " << c[0]
	  << ", cp_low, cp_high )"
	  << std::endl;
      out << std::endl;

    } else if(type == NASA9) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      /* write out single zone */
      out << "        " << fmt( "cp0_R", k ) << " = ";
      writeCoeffTimesVar(false, c9[5],  "", out);
      writeCoeffTimesVar(true,  c9[6],  " * tt0", out);
      writeCoeffTimesVar(true,  c9[7],  " * tt1", out);
      writeCoeffTimesVar(true,  c9[8],  " * tt2", out);
      writeCoeffTimesVar(true,  c9[9],  " * tt3", out);
      writeCoeffTimesVar(true,  c9[3],  " * tt5", out);
      writeCoeffTimesVar(true,  c9[4],  " * tt4", out);
      out << std::endl;

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];

      /* nothing here yet */
      
    }

  }

  out << "        return cp0_R" << std::endl;
  out << std::endl;

}

void PyPrometheus::writeSpeciesEnthalpies(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  writeFunctionName( out, "get_species_enthalpies_RT" );

  /* arguments */
  out << "(self,T):" << std::endl;

  /* declarations */
  out << std::endl;
  out << "        " << "tt0 = T"                   << std::endl;
  out << "        " << "tt1 = T * tt0"             << std::endl;
  out << "        " << "tt2 = T * tt1"             << std::endl;
  out << "        " << "tt3 = T * tt2"             << std::endl;
  out << "        " << "tt4 = np.power( T, -1.0 )" << std::endl;
  out << "        " << "tt5 = tt4 * tt4"           << std::endl;
  out << "        " << "tt6 = np.log(tt0) * tt4"   << std::endl;
  out << std::endl;

  /* allocate buffers */
  out << "        nz,ny,nx = T.shape" << std::endl;
  out << "        h0_RT = np.zeros( [ self.kk,nz,ny,nx ] )" << std::endl;
  out << std::endl;
  
  /* write out the code */

  for(int k = 0; k < m_gas->nSpecies(); ++k) {

    type = m_gas->thermo().speciesThermo().reportType(k);
    
    if(type == NASA) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);

      out << "        h_high  = ";
      writeCoeffTimesVar(false, c[1], "", out);
      writeCoeffTimesVar(true,  c[2], " * 0.50 * tt0", out);
      writeCoeffTimesVar(true,  c[3], " * self.one_third * tt1", out);
      writeCoeffTimesVar(true,  c[4], " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c[5], " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c[6], " * tt4", out);
      out << std::endl;
      
      out << "        h_low   = ";
      writeCoeffTimesVar(false, c[8], "", out);
      writeCoeffTimesVar(true,  c[9], " * 0.50 * tt0",  out);
      writeCoeffTimesVar(true,  c[10], " * self.one_third * tt1", out);
      writeCoeffTimesVar(true,  c[11], " * 0.25 * tt2", out);
      writeCoeffTimesVar(true,  c[12], " * 0.20 * tt3", out);
      writeCoeffTimesVar(true,  c[13], " * tt4", out);
      out << std::endl;
      
      out << "        " << fmt( "h_RT", k )
	  << " = np.where( tt0 < " << c[0]
	  << ", h_low, h_high )"
	  << std::endl;
      out << std::endl;      

    } else if(type == NASA9) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3] *= -1.0;

      /* nothing here yet */      

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];
      
      c9[3] *= -1.0;      

      /* hoting here yet*/
	
      out << std::endl;

    }

  }

  out << "        return h0_RT" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::writeSpeciesEntropies(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  writeFunctionName( out, "get_species_entropies_R" );

  /* arguments */
  out << "(self,T):" << std::endl;;  

  /* declarations */
  out << std::endl;
  out << "        " << "tt0 = T"                   << std::endl;
  out << "        " << "tt1 = T * tt0"             << std::endl;
  out << "        " << "tt2 = T * tt1"             << std::endl;
  out << "        " << "tt3 = T * tt2"             << std::endl;
  out << "        " << "tt4 = np.power( T, -1.0 )" << std::endl;
  out << "        " << "tt5 = tt4 * tt4"           << std::endl;
  out << "        " << "tt6 = np.log(tt0)"         << std::endl;
  out << std::endl;
  
  /* allocate buffers */
  out << "        nz,ny,nx = T.shape" << std::endl;
  out << "        s0_R = np.zeros( [ self.kk,nz,ny,nx ] )" << std::endl;
  out << std::endl;

  /* write out the code */
  
  for(int k = 0; k < m_gas->nSpecies(); ++k) {

    type = m_gas->thermo().speciesThermo().reportType(k);

    if(type == NASA) {
    
      m_gas->thermo().speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);

      out << "        s_high  = ";
      writeCoeffTimesVar(false, c[1], " * tt6", out);
      writeCoeffTimesVar(true,  c[2], " * tt0", out);
      writeCoeffTimesVar(true,  c[3], " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c[4], " * self.one_third * tt2", out);
      writeCoeffTimesVar(true,  c[5], " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c[7], "", out);
      out << std::endl;
      
      out << "        s_low   = ";
      writeCoeffTimesVar(false, c[8], " * tt6", out);
      writeCoeffTimesVar(true,  c[9], " * tt0",  out);
      writeCoeffTimesVar(true,  c[10], " * 0.50 * tt1", out);
      writeCoeffTimesVar(true,  c[11], " * self.one_third * tt2", out);
      writeCoeffTimesVar(true,  c[12], " * 0.25 * tt3", out);
      writeCoeffTimesVar(true,  c[14], "", out);
      out << std::endl;
      
      out << "        " << fmt( "s0_R", k )
	  << " = np.where( tt0 < " << c[0]
	  << ", s_low, s_high )"
	  << std::endl;
      out << std::endl;

    } else if(type == NASA9) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3] *= -1.0;
      c9[4] *= -1.0;

      /* nothing here yet */            

    } else if(type == NASA9MULTITEMP) {

      m_gas->thermo().speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];
      
      c9[3] *= -1.0;
      c9[4] *= -1.0;

      /* nothing here yet */            

    }

  }

  out << "        return s0_R" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::writeSpeciesGibbs(std::ostream& out) {

  /* write function name */
  writeFunctionName( out, "get_species_gibbs_RT" );

  /* arguments */
  out << "(self,T):" << std::endl;

  /* write out the code */
  out << std::endl;
  out << "        " << "h0_RT = self.get_enthalpies_RT( T )" << std::endl;
  out << "        " << "s0_R  = self.get_entropies_R( T )" << std::endl;
  out << "        " << "g0_RT = h0_RT - s0_R" << std::endl;
  out << std::endl;
  out << "        return g0_RT" << std::endl;
  out << std::endl;
  
}


void PyPrometheus::writeEquilibriumConstants(std::ostream& out) {

  std::shared_ptr<Cantera::Reaction> reaction;
  std::vector<int>         reacIndices;
  std::vector<int>         prodIndices;

  /* write function name */
  writeFunctionName( out, "get_equilibrium_constants" );
  
  /* arguments*/
  out << "(self,T):" << std::endl;

  /* declarations*/
  out << std::endl;
  out << "        " << "RT = self.gas_constant * T" << std::endl;
  out << "        " << "C0 = self.one_atm * np.power( RT, -1.0 )" << std::endl;
  out << std::endl;

  /* allocate buffers */
  out << "        nz,ny,nx = T.shape" << std::endl;
  out << "        keq = np.zeros( [ self.ii,nz,ny,nx ] )" << std::endl;
  out << std::endl;

  /* write out the code */
  out << "        " << "g0_RT = self.get_gibbs_RT( T )" << std::endl;
  out << std::endl;
  
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    reaction = m_gas->reaction(i);
    
    if(reaction->reversible == true) {

      getParticipatingSpecies(reaction, reacIndices, prodIndices);
      
      int dn = prodIndices.size() - reacIndices.size();
      std::string revSpecies = increment("g0_RT", prodIndices); 
      std::string fwdSpecies = increment("g0_RT", reacIndices);
 
      out << "        " << fmt("keq",i) << " = ";
      if(dn < 0) { out << " C0 + "; }
      out << "( " << revSpecies << " ) - "
	  << "( " << fwdSpecies;
      if(dn > 0) { for(int n = 0; n < dn; ++n) { out << " + C0 "; } }
      out << " )" << std::endl;
      
    }

  }

  out << std::endl;
  out << "        return keq" << std::endl;
  out << std::endl;

}

void PyPrometheus::getParticipatingSpecies(std::shared_ptr<Cantera::Reaction>& reaction,
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

void PyPrometheus::getThirdBodyEfficiencies(std::shared_ptr<Cantera::Reaction>& reaction,
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

int PyPrometheus::nFalloff() {

  int nfall = 0;
  for(int i = 0; i < m_gas->nReactions(); ++i) {

    if(m_gas->reactionType(i) == Cantera::FALLOFF_RXN) { nfall += 1; }
    
  }

  return(nfall);
  
};

void PyPrometheus::writeFunctionName(std::ostream& out,
				   const std::string& funName) {

  out << "    def " << funName;
  
  
}
