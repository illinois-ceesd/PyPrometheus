#include "PyPrometheus.h"

void PyPrometheus::Greetings() {

  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  
  std::string prometheus = "Prometheus: ";
  std::cout << prometheus << "Welcome..." << std::endl;
  std::cout << prometheus << "Generating Python thermochemistry source code..." << std::endl;
  std::cout << prometheus << "Working with " << m_mech  << " mech..." << std::endl;
  std::cout << prometheus << "Number of species:   "    << m_thermo->nSpecies()     << std::endl;
  std::cout << prometheus << "Number of elements:  "    << m_thermo->nElements()    << std::endl;
  std::cout << prometheus << "Number of reactions: "    << m_kinetics->nReactions() << std::endl;
  
}

void PyPrometheus::WriteMech() {

  std::string   filename;
  std::ofstream out;
  out.setf(std::ios::scientific);
  out.precision(6);

  Greetings();

  /* defs */
  filename = m_mech + ".py";
  out.open( filename.c_str() );
  WriteDefinitions( out );

  ///
  /// State
  ///
  WriteState( out );
  
  ///
  /// Thermo
  ///
  WriteThermo( out );

  ///
  /// Kinetics
  WriteKinetics( out );

  
}

void PyPrometheus::WriteDefinitions(std::ostream& out) {

  double* wts = new double[m_thermo->nSpecies()];
  m_thermo->getMolecularWeights( wts );

  out << "import numpy as np"    << std::endl;
  out << std::endl;
  out << "class prometheus_thermochemistry_kernel():" << std::endl;
  out << std::endl;
  out << "    def __init__(self):" << std::endl;
  out << std::endl;
  out << "        self.model_name    = \"" << m_mech << "\"" << std::endl;
  out << "        self.num_elements  = " << m_thermo->nElements()    << std::endl;
  out << "        self.num_species   = " << m_thermo->nSpecies()     << std::endl;
  out << "        self.num_reactions = " << m_kinetics->nReactions() << std::endl;
  out << "        self.num_falloff   = " << nFalloff() << std::endl;
  out << "        self.one_atm       = 1.01325e5;"     << std::endl;
  out << "        self.one_third     = 1.0 / 3.0;"     << std::endl;
  out << "        self.gas_constant  = 8314.4621;"     << std::endl;
  out << "        self.big_number    = 1.0e300;"       << std::endl;
  out << std::endl;
  out << "        self.wts = np.array( [ ";
  for(int k = 0; k < m_thermo->nSpecies(); ++k) {
    if(k != m_thermo->nSpecies()-1) {
      out << wts[k] << ", ";
    } else {
      out << wts[k];
    }
  }
  out << " ] )" << std::endl;
  out << "        self.iwts = np.reciprocal( self.wts )" << std::endl;
  out << std::endl;
  out << "        return" << std::endl;
  out << std::endl;
  
}

///
/// State
///
void PyPrometheus::WriteDensity(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_density" );
  WriteFunctionArguments( out, { "p", "T", "Y" } );

  out << "        " << "mmw = self.get_mix_molecular_weight( Y )" << std::endl;
  out << "        " << "RT  = self.gas_constant * T" << std::endl;
  out << "        " << "rho = p * mmw / RT" << std::endl;
  out << "        " << "return rho" << std::endl;
  out << std::endl;
    
}

void PyPrometheus::WritePressure(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_pressure" );
  WriteFunctionArguments( out, { "rho", "T", "Y" } );

  out << "        " << "mmw = self.get_mix_molecular_weight( Y )" << std::endl;
  out << "        " << "RT  = self.gas_constant * T" << std::endl;
  out << "        " << "p   = rho * RT / mmw" << std::endl;
  out << "        " << "return p" << std::endl;
  out << std::endl;
    
}

void PyPrometheus::WriteMixMolecularWeight(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_mix_molecular_weight" );
  WriteFunctionArguments( out, { "Y" } );

  out << "        " << "mmw = np.reciprocal( np.dot( self.iwts, Y ) )" << std::endl;
  out << "        " << "return mmw" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteConcentrations(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_concentrations" );
  WriteFunctionArguments( out, { "rho", "Y" } );

  out << "        " << "C = self.iwts * rho * Y" << std::endl;
  out << "        " << "return C" << std::endl;
  out << std::endl;
  
}

///
/// Thermo
///
void PyPrometheus::WriteMixtureSpecificHeatConstantPressure(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_mixture_specific_heat_cp_mass" );
  WriteFunctionArguments( out, { "T", "Y" } );

  /* write out the code */
  out << "        " << "cp0_R = self.get_species_specific_heats_R( T )" << std::endl;
  out << "        " << "cp  = 0.0" << std::endl;
  for(int k = 0; k < m_thermo->nSpecies(); ++k) {
    out << "        cp += "
	<< WriteAccess( "Y", k )     << " * "
	<< WriteAccess( "cp0_R", k ) << " * "
	<< WriteAccess( "self.iwts", k )
	<< std::endl;
  }
  out << "        cp *= self.gas_constant" << std::endl;
  out << std::endl;
  out << "        return cp" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteMixtureSpecificHeatConstantVolume(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_mixture_specific_heat_cv_mass" );
  WriteFunctionArguments( out, { "T", "Y" } );

  /* write out the code */
  out << "        " << "cp0_R = self.get_species_specific_heats_R( T ) - 1.0" << std::endl;
  out << "        " << "cp  = 0.0" << std::endl;
  for(int k = 0; k < m_thermo->nSpecies(); ++k) {
    out << "        cp += "
	<< WriteAccess( "Y", k )     << " * "
	<< WriteAccess( "cp0_R", k ) << " * "
	<< WriteAccess( "self.iwts", k )
	<< std::endl;
  }
  out << "        cp *= self.gas_constant" << std::endl;
  out << std::endl;
  out << "        return cp" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteMixtureEnthalpy(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_mixture_enthalpy_mass" );
  WriteFunctionArguments( out, { "T", "Y" } );

  /* write out the code */
  out << "        " << "h0_RT = self.get_species_enthalpies_RT( T )" << std::endl;
  out << "        " << "h  = 0.0" << std::endl;
  for(int k = 0; k < m_thermo->nSpecies(); ++k) {
    out << "        h += "
	<< WriteAccess( "Y", k )     << " * "
	<< WriteAccess( "h0_RT", k ) << " * "
	<< WriteAccess( "self.iwts", k )
	<< std::endl;
  }
  out << "        h *= self.gas_constant * T" << std::endl;
  out << std::endl;
  out << "        return h" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteMixtureInternalEnergy(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out,  "get_mixture_internal_energy_mass" );
  WriteFunctionArguments( out, { "T", "Y" } );

  /* write out the code */
  out << "        " << "e0_RT = self.get_species_enthalpies_RT( T ) - 1.0" << std::endl;
  out << "        " << "e  = 0.0" << std::endl;
  for(int k = 0; k < m_thermo->nSpecies(); ++k) {
    out << "        e += "
	<< WriteAccess( "Y", k )     << " * "
	<< WriteAccess( "e0_RT", k ) << " * "
	<< WriteAccess( "self.iwts", k )
	<< std::endl;
  }
  out << "        e *= self.gas_constant * T" << std::endl;
  out << std::endl;
  out << "        return e" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteSpeciesSpecificHeats(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  WriteFunctionName( out,  "get_species_specific_heats_R" );
  WriteFunctionArguments( out, { "T" } );

  /* declarations */
  out << "        " << "tt0 = T"                   << std::endl;
  out << "        " << "tt1 = T * tt0"             << std::endl;
  out << "        " << "tt2 = T * tt1"             << std::endl;
  out << "        " << "tt3 = T * tt2"             << std::endl;
  out << "        " << "tt4 = np.power( T, -1.0 )" << std::endl;
  out << "        " << "tt5 = tt4 * tt4"           << std::endl;
  out << std::endl;

  /* allocate buffers */
  out << "        cp0_R = np.zeros( self.num_species )" << std::endl;
  out << std::endl;

  for(int k = 0; k < m_thermo->nSpecies(); ++k) {

    type = m_thermo->speciesThermo().reportType(k);

    if(type == NASA) { 

      m_thermo->speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);

      out << "        cp_high  = ";      
      WriteCoeffTimesVar(false, c[1], "", out);
      WriteCoeffTimesVar(true,  c[2], " * tt0", out);
      WriteCoeffTimesVar(true,  c[3], " * tt1", out);
      WriteCoeffTimesVar(true,  c[4], " * tt2", out);
      WriteCoeffTimesVar(true,  c[5], " * tt3", out);
      out << std::endl;

      out << "        cp_low   = ";
      WriteCoeffTimesVar(false, c[8], "", out);
      WriteCoeffTimesVar(true,  c[9],  " * tt0", out);
      WriteCoeffTimesVar(true,  c[10], " * tt1", out);
      WriteCoeffTimesVar(true,  c[11], " * tt2", out);
      WriteCoeffTimesVar(true,  c[12], " * tt3", out);
      out << std::endl;

      out << "        " << WriteAccess( "cp0_R", k )
	  << " = np.where( tt0 < " << c[0]
	  << ", cp_low, cp_high )"
	  << std::endl;
      out << std::endl;

    } else if(type == NASA9) {

      m_thermo->speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      /* write out single zone */
      out << "        " << WriteAccess( "cp0_R", k ) << " = ";
      WriteCoeffTimesVar(false, c9[5],  "", out);
      WriteCoeffTimesVar(true,  c9[6],  " * tt0", out);
      WriteCoeffTimesVar(true,  c9[7],  " * tt1", out);
      WriteCoeffTimesVar(true,  c9[8],  " * tt2", out);
      WriteCoeffTimesVar(true,  c9[9],  " * tt3", out);
      WriteCoeffTimesVar(true,  c9[3],  " * tt5", out);
      WriteCoeffTimesVar(true,  c9[4],  " * tt4", out);
      out << std::endl;

    } else if(type == NASA9MULTITEMP) {

      m_thermo->speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];

      /* nothing here yet */
      
    }

  }

  out << "        return cp0_R" << std::endl;
  out << std::endl;

}

void PyPrometheus::WriteSpeciesEnthalpies(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  WriteFunctionName( out, "get_species_enthalpies_RT" );
  WriteFunctionArguments( out, { "T" } );

  /* declarations */
  out << "        " << "tt0 = T"                   << std::endl;
  out << "        " << "tt1 = T * tt0"             << std::endl;
  out << "        " << "tt2 = T * tt1"             << std::endl;
  out << "        " << "tt3 = T * tt2"             << std::endl;
  out << "        " << "tt4 = np.power( T, -1.0 )" << std::endl;
  out << "        " << "tt5 = tt4 * tt4"           << std::endl;
  out << "        " << "tt6 = np.log(tt0) * tt4"   << std::endl;
  out << std::endl;

  /* allocate buffers */
  out << "        h0_RT = np.zeros( self.num_species )" << std::endl;
  out << std::endl;
  
  /* write out the code */

  for(int k = 0; k < m_thermo->nSpecies(); ++k) {

    type = m_thermo->speciesThermo().reportType(k);
    
    if(type == NASA) {

      m_thermo->speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);

      out << "        h_high  = ";
      WriteCoeffTimesVar(false, c[1], "", out);
      WriteCoeffTimesVar(true,  c[2], " * 0.50 * tt0", out);
      WriteCoeffTimesVar(true,  c[3], " * self.one_third * tt1", out);
      WriteCoeffTimesVar(true,  c[4], " * 0.25 * tt2", out);
      WriteCoeffTimesVar(true,  c[5], " * 0.20 * tt3", out);
      WriteCoeffTimesVar(true,  c[6], " * tt4", out);
      out << std::endl;
      
      out << "        h_low   = ";
      WriteCoeffTimesVar(false, c[8], "", out);
      WriteCoeffTimesVar(true,  c[9], " * 0.50 * tt0",  out);
      WriteCoeffTimesVar(true,  c[10], " * self.one_third * tt1", out);
      WriteCoeffTimesVar(true,  c[11], " * 0.25 * tt2", out);
      WriteCoeffTimesVar(true,  c[12], " * 0.20 * tt3", out);
      WriteCoeffTimesVar(true,  c[13], " * tt4", out);
      out << std::endl;
      
      out << "        " << WriteAccess( "h0_RT", k )
	  << " = np.where( tt0 < " << c[0]
	  << ", h_low, h_high )"
	  << std::endl;
      out << std::endl;      

    } else if(type == NASA9) {

      m_thermo->speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3] *= -1.0;

      /* nothing here yet */      

    } else if(type == NASA9MULTITEMP) {

      m_thermo->speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];
      
      c9[3] *= -1.0;      

      /* hoting here yet*/
	
      out << std::endl;

    }

  }

  out << "        return h0_RT" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteSpeciesEntropies(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);

  /* write function name */
  WriteFunctionName( out, "get_species_entropies_R" );
  WriteFunctionArguments( out, { "T" } );

  /* declarations */
  out << "        " << "tt0 = T"                   << std::endl;
  out << "        " << "tt1 = T * tt0"             << std::endl;
  out << "        " << "tt2 = T * tt1"             << std::endl;
  out << "        " << "tt3 = T * tt2"             << std::endl;
  out << "        " << "tt4 = np.power( T, -1.0 )" << std::endl;
  out << "        " << "tt5 = tt4 * tt4"           << std::endl;
  out << "        " << "tt6 = np.log(tt0)"         << std::endl;
  out << std::endl;
  
  /* allocate buffers */
  out << "        s0_R = np.zeros( self.num_species )" << std::endl;
  out << std::endl;

  /* write out the code */
  
  for(int k = 0; k < m_thermo->nSpecies(); ++k) {

    type = m_thermo->speciesThermo().reportType(k);

    if(type == NASA) {
    
      m_thermo->speciesThermo().reportParams(k, type, &c[0], minTemp, maxTemp, refPres);

      out << "        s_high  = ";
      WriteCoeffTimesVar(false, c[1], " * tt6", out);
      WriteCoeffTimesVar(true,  c[2], " * tt0", out);
      WriteCoeffTimesVar(true,  c[3], " * 0.50 * tt1", out);
      WriteCoeffTimesVar(true,  c[4], " * self.one_third * tt2", out);
      WriteCoeffTimesVar(true,  c[5], " * 0.25 * tt3", out);
      WriteCoeffTimesVar(true,  c[7], "", out);
      out << std::endl;
      
      out << "        s_low   = ";
      WriteCoeffTimesVar(false, c[8], " * tt6", out);
      WriteCoeffTimesVar(true,  c[9], " * tt0",  out);
      WriteCoeffTimesVar(true,  c[10], " * 0.50 * tt1", out);
      WriteCoeffTimesVar(true,  c[11], " * self.one_third * tt2", out);
      WriteCoeffTimesVar(true,  c[12], " * 0.25 * tt3", out);
      WriteCoeffTimesVar(true,  c[14], "", out);
      out << std::endl;
      
      out << "        " << WriteAccess( "s0_R", k )
	  << " = np.where( tt0 < " << c[0]
	  << ", s_low, s_high )"
	  << std::endl;
      out << std::endl;

    } else if(type == NASA9) {

      m_thermo->speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      c9[3] *= -1.0;
      c9[4] *= -1.0;

      /* nothing here yet */            

    } else if(type == NASA9MULTITEMP) {

      m_thermo->speciesThermo().reportParams(k, type, &c9[0], minTemp, maxTemp, refPres);

      int nzones = c9[0];
      
      c9[3] *= -1.0;
      c9[4] *= -1.0;

      /* nothing here yet */            

    }

  }

  out << "        return s0_R" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteSpeciesGibbs(std::ostream& out) {

  /* write function name */
  WriteFunctionName( out, "get_species_gibbs_RT" );
  WriteFunctionArguments( out, { "T" } );

  /* write out the code */
  out << "        " << "h0_RT = self.get_species_enthalpies_RT( T )" << std::endl;
  out << "        " << "s0_R  = self.get_species_entropies_R( T )" << std::endl;
  out << "        " << "g0_RT = h0_RT - s0_R" << std::endl;
  out << std::endl;
  out << "        return g0_RT" << std::endl;
  out << std::endl;
  
}


void PyPrometheus::WriteEquilibriumConstants(std::ostream& out) {

  std::shared_ptr<Cantera::Reaction> reaction;
  std::vector<int>         reacIndices;
  std::vector<int>         prodIndices;

  /* write function name */
  WriteFunctionName( out, "get_equilibrium_constants" );
  WriteFunctionArguments( out, { "T" } );

  /* declarations */
  out << "        " << "RT = self.gas_constant * T" << std::endl;
  out << "        " << "C0 = np.log( self.one_atm / RT )" << std::endl;
  out << std::endl;

  /* allocate buffers */
  out << "        k_eq = np.zeros( self.num_reactions )" << std::endl;
  out << std::endl;

  /* write out the code */
  out << "        " << "g0_RT = self.get_species_gibbs_RT( T )" << std::endl;
  out << std::endl;
  
  for(int i = 0; i < m_kinetics->nReactions(); ++i) {

    reaction = m_kinetics->reaction(i);
    
    if(reaction->reversible == true) {

      GetParticipatingSpecies(reaction, reacIndices, prodIndices);
      
      int dn = prodIndices.size() - reacIndices.size();
      std::string revSpecies = Increment("g0_RT", prodIndices); 
      std::string fwdSpecies = Increment("g0_RT", reacIndices);
 
      out << "        " << WriteAccess("k_eq",i) << " = ";
      if(dn < 0) { out << " C0 + "; }
      out << "( " << revSpecies << " ) - "
	  << "( " << fwdSpecies;
      if(dn > 0) { for(int n = 0; n < dn; ++n) { out << " + C0 "; } }
      out << " )" << std::endl;
      
    }

  }

  out << std::endl;
  out << "        return k_eq" << std::endl;
  out << std::endl;

}

void PyPrometheus::WriteNewtonTemperature(std::ostream& out) {

  /* write function name */
  std::vector<std::string> arguments;
  arguments = { "H_or_E", "T_guess", "Y", "do_energy = False" };
  
  WriteFunctionName( out, "get_temperature" );
  WriteFunctionArguments( out, arguments );
  
  /* declarations */  
  out << "        if do_energy == False:" << std::endl;
  out << "            pv_fun = self.get_mixture_specific_heat_cp_mass" << std::endl;
  out << "            he_fun = self.get_mixture_enthalpy_mass" << std::endl;
  out << "        else:" << std::endl;
  out << "            pv_fun = self.get_mixture_specific_heat_cv_mass" << std::endl;
  out << "            he_fun = self.get_mixture_internal_energy_mass" << std::endl;
  out << std::endl;
  out << "        num_iter = 500" << std::endl;
  out << "        tol = 1.0e-6"   << std::endl;
  out << "        T_i = T_guess" << std::endl;
  out << "        dT = 1.0"    << std::endl;
  out << "        F  = H_or_E" << std::endl;
  out << "        J  = 0.0"    << std::endl;
  out << std::endl;
  
  out << "        for iter in range( 0, num_iter ):"   << std::endl;
  out << "            F    -= he_fun( T_i, Y )" << std::endl;
  out << "            J    -= pv_fun( T_i, Y )" << std::endl;
  out << "            dT    = - F / J"           << std::endl;
  out << "            T_i  += dT"                << std::endl;
  out << "            if np.abs( dT ) < tol:"    << std::endl;
  out << "                T = T_i"               << std::endl;
  out << "                break"                 << std::endl;
  out << "            F = H_or_E" << std::endl;
  out << "            J = 0.0"    << std::endl;
  out << std::endl;
  out << "        T = T_i" << std::endl;
  out << std::endl;
  out << "        return T" << std::endl;
  out << std::endl;
}

///
/// Kinetics
///
void PyPrometheus::WriteFalloffKinetics(std::ostream& out) {

  std::vector<int>                          falloffIndices;
  std::vector<int>                          falloffType;
  std::vector<std::vector<double> >         falloffEfficiencies;
  std::vector<std::vector<double> >         falloffParams;
  std::shared_ptr<Cantera::Reaction>        reaction;
  std::shared_ptr<Cantera::FalloffReaction> work;

  WriteFunctionName( out, "get_falloff_rates" );
  WriteFunctionArguments( out, { "T", "C", "k_fwd" } );

  /* declarations */
  out << "        TROE = 110" << std::endl;
  out << "        log_T = np.log( T )" << std::endl;
  out << "        inv_T = 1.0 / T"     << std::endl;
  out << "        k_hi = np.zeros( self.num_falloff )" << std::endl;
  out << "        k_lo = np.zeros( self.num_falloff )" << std::endl;
  out << "        pr   = np.zeros( self.num_falloff )" << std::endl;
  out << "        work = np.zeros( self.num_falloff )" << std::endl;
  out << "        falloff_type = 100 * np.ones( self.num_falloff, dtype = int )" << std::endl;
  out << std::endl;
  
  /* identify falloff reactions */
  for(int i = 0; i < m_kinetics->nReactions(); ++i) {

    reaction = m_kinetics->reaction(i);
    if(reaction->reaction_type == Cantera::FALLOFF_RXN) { falloffIndices.push_back(i); }

  }

  /* Parameterization */
  for(int i = 0; i < falloffIndices.size(); ++i) {

    reaction = m_kinetics->reaction( falloffIndices[i] );
    work     = std::dynamic_pointer_cast<Cantera::FalloffReaction>(reaction);

    /* Arrhenius expressions */
    double Af  = work->high_rate.preExponentialFactor();
    double bf  = work->high_rate.temperatureExponent();
    double Taf = (-1) * work->high_rate.activationEnergy_R();

    double A0  = work->low_rate.preExponentialFactor();
    double b0  = work->low_rate.temperatureExponent();
    double Ta0 = (-1) * work->low_rate.activationEnergy_R();

    out << "        " << WriteAccess("k_hi",i) << " = np.exp(" << std::log(Af);
    if(bf  != 0.0) { WriteCoeffTimesVar(true, bf, " * log_T", out); }
    if(Taf != 0.0) { WriteCoeffTimesVar(true, Taf," * inv_T", out); }
    out << ")" << std::endl;

    out << "        " << WriteAccess("k_lo",i) << " = np.exp(" << std::log(A0);
    if(b0  != 0.0) { WriteCoeffTimesVar(true, b0, " * log_T", out); }
    if(Ta0 != 0.0) { WriteCoeffTimesVar(true, Ta0," * inv_T", out); }
    out << ")" << std::endl;
    out << std::endl;

    /* Third-body efficiencies */
    std::vector<double> efficiencies;
    for(int k = 0; k < m_thermo->nSpecies(); ++k) {
      const std::string name = m_thermo->speciesName(k);
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

    out << "        " << WriteAccess("pr",i) << " = ";

    const std::string Ck   = WriteAccess("C",0);
    double eff = falloffEfficiencies[i][0];
    if(eff != 0.0) { WriteCoeffTimesVar(false, eff, " * "+Ck, out); }
    
    for(int k = 1; k < m_thermo->nSpecies(); ++k) {
      const std::string Ck   = WriteAccess("C",k);
      double eff = falloffEfficiencies[i][k];
      if(eff != 0.0) { out << " + "; WriteCoeffTimesVar(false, eff, " * "+Ck, out); }
    }

    out << " " << std::endl;

  }
  out << std::endl;
  out << "        for i_falloff in range( 0, self.num_falloff ):" << std::endl;
  out << "            pr[i_falloff] *= ( k_lo[i_falloff] / k_hi[i_falloff] )" << std::endl;
  out << std::endl;

  /* falloff function */
  for(int i = 0; i < falloffIndices.size(); ++i) {
    out << "        " << WriteAccess("falloff_type",i)
	<< " = " << falloffType[i] << std::endl;    
  }
  out << std::endl;

  /* falloff rate */
  out.precision(6);
  for(int i = 0; i < falloffIndices.size(); ++i) {

    if( falloffType[i] == Cantera::SIMPLE_FALLOFF ) {

      out << "        " << WriteAccess("work",i) << " = 1.0;" << std::endl;

    } else if( falloffType[i] == Cantera::TROE_FALLOFF ) {

      std::vector<double> params = falloffParams[i];
      params[1] = (-1) * 1.0 / params[1];
      params[2] = (-1) * 1.0 / params[2];
      params[3] = (-1) * params[3];
    
      out << "        " << WriteAccess("work",i) << " = ";
      out << "(1.0 - " << params[0] << ") * np.exp(";
      WriteCoeffTimesVar(true,params[1]," * T",out);
      out << ")";
      out << " + " << params[0] << " * np.exp(";
      WriteCoeffTimesVar(true,params[2]," * T",out);
      out << ")";
      if(params[3] != 0.0) {
	out << " + np.exp(";
	WriteCoeffTimesVar(true,params[3]," * inv_T",out);
	out<< ")" << std::endl;
      } else {
	out << std::endl;
      }

    }
    
  }
  out << std::endl;

  /* falloff calculation */
  out << "        for i_falloff in range( 0, self.num_falloff ):" << std::endl;
  out << "            lpr = np.log10(pr[i_falloff])"       << std::endl;
  out << "            if falloff_type[i_falloff] == TROE:" << std::endl;
  out << "                cc = -0.40 - 0.67 * np.log10(work[i_falloff])" << std::endl;
  out << "                nn =  0.75 - 1.27 * np.log10(work[i_falloff])" << std::endl;
  out << "                f1 =  (lpr + cc)/(nn - 0.14 * (lpr + cc))"     << std::endl;
  out << "                work[i_falloff] = np.log10(work[i_falloff])/(1 + f1 * f1)"
      << std::endl;
  out << "                work[i_falloff] = 10.0 ** work[i_falloff]"
      << std::endl;
  out << "            work[i_falloff] = (pr[i_falloff] * work[i_falloff])/(1 + pr[i_falloff])"
      << std::endl;
  out << std::endl;

  /* fwd rate coefficients */
  for(int i = 0; i < falloffIndices.size(); ++i) {
    int rxn = falloffIndices[i];
    out << "        "
	<< WriteAccess("k_fwd",rxn) << " = "
	<< WriteAccess("k_hi",i)    << " * "
	<< WriteAccess("work",i)    << std::endl;
  }
  
  out << std::endl;
  out << "        return" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteRateCoefficients(std::ostream& out) {

  std::vector<int>                   tbrIndices;
  std::vector<std::vector<double> >  tbrEfficiencies;
  std::shared_ptr<Cantera::Reaction> reaction;

  WriteFunctionName( out, "get_rate_coefficients" );
  WriteFunctionArguments( out, { "T", "C" } );

  /* declarations */
  out << "        log_T = np.log( T )" << std::endl;
  out << "        inv_T = 1.0 / T"     << std::endl;
  out << "        k_eq  = self.get_equilibrium_constants( T )" << std::endl;
  out << "        k_fwd = np.zeros( self.num_reactions )"   << std::endl;
  out << "        k_rev = np.zeros( self.num_reactions )"   << std::endl;
  out << std::endl;
  
  /* Arrhenius kinetics */
  for(int i = 0; i < m_kinetics->nReactions(); ++i) {

    reaction = m_kinetics->reaction(i);

    switch(reaction->reaction_type) {

    case Cantera::ELEMENTARY_RXN:

      WriteArrheniusKinetics(i, reaction, out);
      break;
      
    case Cantera::THREE_BODY_RXN:

      tbrIndices.push_back(i);
      GetThirdBodyEfficiencies(reaction, tbrEfficiencies);
      WriteArrheniusKinetics(i, reaction, out);
      break;
      
    }

  }
  out << std::endl;

  /* Three-body reactions */
  out.precision(3);
  for(int i = 0; i < tbrIndices.size(); ++i) {

    out << "        " << WriteAccess("k_fwd",tbrIndices[i]) << " *= ( ";
    
    const std::string Ck = WriteAccess("C",0);
    double eff = tbrEfficiencies[i][0];
    if(eff != 0.0) { WriteCoeffTimesVar(false, eff, " * "+Ck, out); }
    
    for(int k = 1; k < m_thermo->nSpecies(); ++k) {
      const std::string Ck = WriteAccess("C",k);
      double eff = tbrEfficiencies[i][k];
      if(eff != 0.0) { WriteCoeffTimesVar(true, eff, " * "+Ck, out); }
    }
    
    out << " ) " << std::endl;
    
  }
  out << std::endl;

  if( nFalloff() > 0 ) {
    out << "        self.get_falloff_rates(T, C, k_fwd)" << std::endl;
  }
  out << std::endl;
  out << "        for i in range( 0, self.num_reactions ):" << std::endl;
  out << "            if k_eq[i] > self.big_number:" << std::endl;
  out << "                k_eq[i] = self.big_number" << std::endl;
  out << "            k_rev[i] = k_fwd[i] * np.exp( k_eq[i] )" << std::endl;
  out << std::endl;
  out << "        return k_fwd, k_rev" << std::endl;
  out << std::endl;
  
}

void PyPrometheus::WriteArrheniusKinetics(int& rxn,
					  std::shared_ptr<Cantera::Reaction>& reaction,
					  std::ostream& out) {

  out.precision(6);
  
  std::shared_ptr<Cantera::ElementaryReaction> work;
  work = std::dynamic_pointer_cast<Cantera::ElementaryReaction>(reaction);

  double A  = work->rate.preExponentialFactor();
  double b  = work->rate.temperatureExponent();
  double Ta = (-1) * work->rate.activationEnergy_R();

  out << "        " << WriteAccess("k_fwd",rxn) << " = np.exp(" << std::log(A);
  if(b  != 0.0) { WriteCoeffTimesVar(true, b, " * log_T", out); }
  if(Ta != 0.0) { WriteCoeffTimesVar(true, Ta," * inv_T", out); }
  out << ")" << std::endl;
  
}

void PyPrometheus::WriteNetRatesOfProgress(std::ostream& out) {

  std::vector<int> reacIndices;
  std::vector<int> prodIndices;
  std::shared_ptr<Cantera::Reaction> reaction;

  WriteFunctionName( out, "get_net_rates_of_progress" );
  WriteFunctionArguments( out, { "T", "C" } );

  /* declarations */
  out << "        R_fwd = np.zeros( self.num_reactions )" << std::endl;
  out << "        R_rev = np.zeros( self.num_reactions )" << std::endl;
  out << "        R_net = np.zeros( self.num_reactions )" << std::endl;
  out << "        k_fwd, k_rev = self.get_rate_coefficients( T, C )" << std::endl;
  out << std::endl;

  for(int i = 0; i < m_kinetics->nReactions(); ++i) {

    reaction = m_kinetics->reaction(i);
    GetParticipatingSpecies(reaction, reacIndices, prodIndices);
    std::string fwdSpecies = Multiply("C",reacIndices);
    std::string revSpecies = Multiply("C",prodIndices);

    out << "        " << WriteAccess("R_fwd",i)
	<< " = "  << WriteAccess("k_fwd",i)
	<< " * "  << fwdSpecies
	<< std::endl;
    if(reaction->reversible == true) {
      out << "        " << WriteAccess("R_rev",i)
	  << " = "  << WriteAccess("k_rev",i)
	  << " * "  << revSpecies
	  << std::endl;
    }
    out << std::endl;
    
  }

  out << "        for i in range( 0, self.num_reactions ):" << std::endl;
  out << "            R_net[i] = R_fwd[i] - R_rev[i]" << std::endl;
  out << std::endl;
  out << "        return R_net" << std::endl;
  out << std::endl;    
  
}

void PyPrometheus::WriteNetProductionRates(std::ostream& out) {

  std::vector<int> reacIndices;
  std::vector<int> prodIndices;
  std::shared_ptr<Cantera::Reaction> reaction;

  for(int k = 0; k < m_thermo->nSpecies(); ++k) {
    m_rhs.insert( std::pair<int, std::vector<std::string> >(k, std::vector<std::string>()) );
  }

  WriteFunctionName( out, "get_net_production_rates" );
  WriteFunctionArguments( out, { "rho", "T", "Y" } );

  /* declarations */
  //out << "        rho = self.get_density( p, T, Y )"        << std::endl;
  out << "        C   = self.get_concentrations( rho, Y )"  << std::endl;
  out << "        R_net = self.get_net_rates_of_progress( T, C )" << std::endl;
  out << "        omega = np.zeros( self.num_species )" << std::endl;
  out << std::endl;
  
  for(int i = 0; i < m_kinetics->nReactions(); ++i) {

    reaction = m_kinetics->reaction(i);
    GetParticipatingSpecies(reaction, reacIndices, prodIndices);

    for(int k = 0; k < reacIndices.size(); ++k) {
      std::string rhs = " - " + WriteAccess("R_net",i);
      m_rhs[reacIndices[k]].push_back(rhs);
    }

    for(int k = 0; k < prodIndices.size(); ++k) {
      std::string rhs = " + " + WriteAccess("R_net",i); 
      m_rhs[prodIndices[k]].push_back(rhs);
    }
    
  }

  for(std::map<int, std::vector<std::string> >::iterator iter = m_rhs.begin();
      iter != m_rhs.end();
      ++iter)
    {
      std::string rhs = WrapString(iter->second);
      if(rhs != "") {
	out << "        " << WriteAccess("omega",iter->first)
	    << " = "  << rhs
	    << std::endl;
      }
    }
  out << std::endl;
  out << "        return omega" << std::endl;
  out << std::endl;
  
  
}

void PyPrometheus::WriteFunctionName(std::ostream& out,
				     const std::string& funName) {

  out << "    def " << funName;
  
  
}

void PyPrometheus::WriteFunctionArguments(std::ostream& out,
					  const std::vector<std::string>& args) {

  int numArgs = args.size();
  
  out << "(self, ";
  for(int iArg = 0; iArg < numArgs; ++iArg) {
    if( iArg < numArgs - 1 ) {
      out << args[ iArg ] << ", ";
    } else {
      out << args[ iArg ] << "):";
    }    
  }
  out << "\n" << std::endl;
  
}
