#include "cantera/IdealGasMix.h"
#include "cantera/kinetics.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "utils.h"

class wolfChemGen
{
 public:
    
  wolfChemGen(std::string& mech, Cantera::IdealGasMix& gas) {
    m_gas  = &gas;
    m_mech = mech;
  };

  void greetings();
  
  void writeMech();

  void writeKinetics(std::ostream& out);
  
  void writeThermo(std::ostream& out);

  void writeFalloffKinetics(std::ostream& out);

  void writeRateCoefficients(std::ostream& out);

  void writeArrheniusKinetics(int& rxn,
			      std::shared_ptr<Cantera::Reaction>& reaction,
			      std::ostream& out);

  void writeNetRatesOfProgress(std::ostream& out);
  
  void writeNetProductionRates(std::ostream& out);

  void writeSpeciesSpecificHeats(std::ostream& out);

  void writeSpeciesEnthalpies(std::ostream& out);

  void writeSpeciesEntropies(std::ostream& out);

  void writeSpeciesGibbs(std::ostream& out);

  void writeEquilibriumConstants(std::ostream& out);

  void getThirdBodyEfficiencies(std::shared_ptr<Cantera::Reaction>& reaction,
				std::vector<std::vector<double> >& tbrEfficiencies);
  
  void getParticipatingSpecies(std::shared_ptr<Cantera::Reaction>& reaction,
			       std::vector<int>& reacIndices,
			       std::vector<int>& prodIndices);

  int nFalloff();
  
 private:
  std::string m_mech;
  std::map<int, std::vector<std::string> > m_rhs;
  Cantera::IdealGasMix* m_gas;
    
};
