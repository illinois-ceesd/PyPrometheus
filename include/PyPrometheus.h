#ifndef PYPROM
#define PYPROM

#include "Prometheus.h"

class PyPrometheus : public Prometheus
{
 public:
    
  PyPrometheus(Config& config,
	       std::shared_ptr<Cantera::Solution> solution)
    {
	/* save model name */      
	m_mech = config.mech();
	/* save fixed thermo variables */
	m_thermoXY = config.fixedThermoVars();
	/* pointer to managers (ugh) */
	m_thermo   = solution->thermo();
	m_kinetics = solution->kinetics();
	if(solution->transport() != NULL){
	  m_transport = solution->transport();
	}
	
      };

  void Greetings();
  void WriteMech();

  void WriteDefinitions(std::ostream& out);
  //void WriteKineticsHeader(std::ostream& out);
  //void WriteThermoHead(std::ostream& out);

  void WriteDensity(std::ostream& out);
  void WritePressure(std::ostream& out);
  void WriteMixMolecularWeight(std::ostream& out);
  void WriteConcentrations(std::ostream& out);

  void WriteMixtureGasConstant(std::ostream& out);
  void WriteMixtureSpecificHeatConstantPressure(std::ostream& out);
  void WriteMixtureSpecificHeatConstantVolume(std::ostream& out);
  void WriteMixtureEnthalpy(std::ostream& out);
  void WriteMixtureInternalEnergy(std::ostream& out);
  void WriteSpeciesSpecificHeats(std::ostream& out);
  void WriteSpeciesEnthalpies(std::ostream& out);
  void WriteSpeciesEntropies(std::ostream& out);
  void WriteSpeciesGibbs(std::ostream& out);
  void WriteEquilibriumConstants(std::ostream& out);
  void WriteNewtonTemperature(std::ostream& out);

  void WriteFalloffKinetics(std::ostream& out);
  void WriteRateCoefficients(std::ostream& out);
  void WriteArrheniusKinetics(int& rxn,
			      std::shared_ptr<Cantera::Reaction>& reaction,
			      std::ostream& out);
  void WriteNetRatesOfProgress(std::ostream& out);
  void WriteNetProductionRates(std::ostream& out);

  void WriteFunctionName(std::ostream& out,
			 const std::string& funName);
  void WriteFunctionArguments(std::ostream& out,
			      const std::vector<std::string>& args);
  
};

#endif
