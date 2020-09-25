#ifndef PROMETHEUS
#define PROMETHEUS

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "Utils.h"
#include "Config.h"

class Prometheus
{
 public:
    
  Prometheus() { };

  virtual void Greetings() { ErrorMessage( "Greetings" ); };  
  virtual void WriteMech() { ErrorMessage( "WriteMech"); };

  ///
  /// Definitions and class headers
  ///
  virtual void WriteDefinitions(std::ostream& out) { ErrorMessage( "WriteDefinitions" ); };
  virtual void WriteKineticsHeader(std::ostream& out) { ErrorMessage( "WriteKineticsHeader" ); };
  virtual void WriteThermoHeader(std::ostream& out) { ErrorMessage( "WriteThermoHeader" ); };

  ///
  /// Write drivers
  ///
  void WriteState(std::ostream& out) {

    std::string prometheus = "Prometheus:State: ";

    std::cout << prometheus << "Writing density ... " << std::endl;
    WriteDensity( out );

    std::cout << prometheus << "Writing pressure ... " << std::endl;
    WritePressure( out );

    std::cout << prometheus << "Writing mixture molecular weight ... " << std::endl;
    WriteMixMolecularWeight( out );
    
    std::cout << prometheus << "Writing concentrations ... " << std::endl;
    WriteConcentrations( out );
    
  }

  void WriteThermo(std::ostream& out) {

    std::string prometheus = "Prometheus:Thermo: ";

    std::cout << prometheus << "Writing mixture specific heat (constant pressure) ... " << std::endl;
    WriteMixtureSpecificHeatConstantPressure( out );

    std::cout << prometheus << "Writing mixture specific heat (constant volume) ... " << std::endl;
    WriteMixtureSpecificHeatConstantVolume( out );

    std::cout << prometheus << "Writing mixture enthalpy ... " << std::endl;
    WriteMixtureEnthalpy( out );

    std::cout << prometheus << "Writing mixture internal energy ... " << std::endl;
    WriteMixtureInternalEnergy( out );
    
    std::cout << prometheus << "Writing specific heats ... " << std::endl;
    WriteSpeciesSpecificHeats( out );
  
    std::cout << prometheus << "Writing enthalpies ... " << std::endl;
    WriteSpeciesEnthalpies( out );

    std::cout << prometheus << "Writing entropies ... " << std::endl;
    WriteSpeciesEntropies( out );

    std::cout << prometheus << "Writing Gibbs functions ... " << std::endl;
    WriteSpeciesGibbs( out );

    std::cout << prometheus << "Writing equilibrium constants ... " << std::endl;
    WriteEquilibriumConstants( out );

    std::cout << prometheus << "Writing Newton method for temperature inversion ... " << std::endl;
    WriteNewtonTemperature( out );

  };

  
  void WriteKinetics(std::ostream& out) {

    std::string prometheus = "Prometheus:Kinetics: ";

    if(nFalloff() > 0) {
      std::cout << prometheus << "Writing falloff kinetics ... " << std::endl;
      WriteFalloffKinetics( out );
    }

    std::cout << prometheus << "Writing rate coefficients ... " << std::endl;
    WriteRateCoefficients( out );

    std::cout << prometheus << "Writing reaction rates ... " << std::endl;
    WriteNetRatesOfProgress( out );

    std::cout << prometheus << "Writing species production rates ... " << std::endl;
    WriteNetProductionRates( out );
      
  };

  
  void WriteTransport(std::ostream& out) { ErrorMessage( "WriteTransport" ); };

  ///
  /// State
  ///
  virtual void WriteDensity(std::ostream& out) {
    ErrorMessage( "WriteDensity" );
  };
  virtual void WritePressure(std::ostream& out) {
    ErrorMessage( "WritePressure" );
  };
  virtual void WriteMixMolecularWeight(std::ostream& out) {
    ErrorMessage( "WriteMixMolecularWeight" );
  };
  virtual void WriteConcentrations(std::ostream& out) {
    ErrorMessage( "WriteConcetrations" );
  };
  
  ///
  /// Thermo
  ///
  virtual void WriteMixtureSpecificHeatConstantPressure(std::ostream& out) {
    ErrorMessage( "WriteMixtureSpecificHeatConstantPressure" );
  };
  virtual void WriteMixtureSpecificHeatConstantVolume(std::ostream& out) {
    ErrorMessage( "WriteMixtureSpecificHeatConstantVolume" );
  };
  virtual void WriteMixtureEnthalpy(std::ostream& out) {
    ErrorMessage( "WriteMixtureEnthalpy" );
  };
  virtual void WriteMixtureInternalEnergy(std::ostream& out) {
    ErrorMessage( "WriteMixtureInternalEnergy" );
  };
  virtual void WriteSpeciesSpecificHeats(std::ostream& out) {
    ErrorMessage( "WriteSpeciesSpecificHeats" );
  };
  virtual void WriteSpeciesEnthalpies(std::ostream& out) {
    ErrorMessage( "WriteSpeciesEnthalpies" );
  };
  virtual void WriteSpeciesEntropies(std::ostream& out) {
    ErrorMessage( "WriteSpeciesEntropies" );
  };
  virtual void WriteSpeciesGibbs(std::ostream& out) {
    ErrorMessage( "WriteSpeciesGibbs" );
  };
  virtual void WriteEquilibriumConstants(std::ostream& out) {
    ErrorMessage( "WriteEquilibriumConstants" );
  };
  virtual void WriteNewtonTemperature(std::ostream& out) {
    ErrorMessage( "WriteNewtonTemperature" );
  };
  
  ///
  /// Kinetics
  ///
  virtual void WriteFalloffKinetics(std::ostream& out) {
    ErrorMessage( "WriteFalloffKinetics" );
  };
  virtual void WriteRateCoefficients(std::ostream& out) {
    ErrorMessage( "WriteRateCoefficients" );
  };
  virtual void WriteArrheniusKinetics(int& rxn,
				      std::shared_ptr<Cantera::Reaction>& reaction,
				      std::ostream& out) {
    ErrorMessage( "WriteArrheniusKinetics" );
  };
  virtual void WriteNetRatesOfProgress(std::ostream& out) {
    ErrorMessage( "WriteNetRatesOfProgress" );
  };
  virtual void WriteNetProductionRates(std::ostream& out) {
    ErrorMessage( "WriteNetRatesOfProgress" );
  };

  void GetThirdBodyEfficiencies(std::shared_ptr<Cantera::Reaction>& reaction,
				std::vector<std::vector<double> >& tbrEfficiencies) {

    std::vector<double> efficiencies;
    std::shared_ptr<Cantera::ThreeBodyReaction> work;
    work = std::dynamic_pointer_cast<Cantera::ThreeBodyReaction>(reaction);
  
    for(int k = 0; k < m_thermo->nSpecies(); ++k) {

      const std::string name = m_thermo->speciesName(k);
      double eff = work->third_body.efficiency(name);
      efficiencies.push_back(eff);
    
    }

    tbrEfficiencies.push_back(efficiencies);
    
  };
  
  void GetParticipatingSpecies(std::shared_ptr<Cantera::Reaction>& reaction,
			       std::vector<int>& reacIndices,
			       std::vector<int>& prodIndices) {

    reacIndices.clear();
    prodIndices.clear(); 
  
    Cantera::Composition reactantsMap = reaction->reactants;
    Cantera::Composition productsMap  = reaction->products;

    for(Cantera::Composition::const_iterator iter = reactantsMap.begin();
	iter != reactantsMap.end();
	++iter) {
      reacIndices.push_back( m_thermo->speciesIndex(iter->first) );
      if(iter->second > 1.0) {
	int nu = (int)iter->second;      
	for(int k = 0; k < nu-1; ++k) {
	  reacIndices.push_back( m_thermo->speciesIndex(iter->first) );
	}
      }
    }

    for(Cantera::Composition::const_iterator iter = productsMap.begin();
	iter != productsMap.end();
	++iter) {
      prodIndices.push_back( m_thermo->speciesIndex(iter->first) );
      if(iter->second > 1.0) {
	int nu = (int)iter->second;     
	for(int k = 0; k < nu-1; ++k) {
	  prodIndices.push_back( m_thermo->speciesIndex(iter->first) );
	}
      }
    }
    
  }

  int nFalloff() {
    int nfall = 0;
    for(int i = 0; i < m_kinetics->nReactions(); ++i) {
      if( m_kinetics->reactionType(i) == Cantera::FALLOFF_RXN ) { nfall += 1; }
    }
    return(nfall);
  }

  ///
  /// Error message
  ///
  void ErrorMessage(const std::string& method) {
    std::cout << "Prometheus::" + method + " called from base class..." << std::endl;
  }
  
 protected:
  std::string m_mech;
  std::string m_thermoXY;
  std::map<int, std::vector<std::string> > m_rhs;
  std::shared_ptr<Cantera::ThermoPhase> m_thermo;
  std::shared_ptr<Cantera::Kinetics>    m_kinetics;
  std::shared_ptr<Cantera::Transport>   m_transport;
    
};

#endif
