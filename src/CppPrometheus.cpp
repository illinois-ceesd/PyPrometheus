#include "CppPrometheus.h"

void CppPrometheus::Greetings() {

  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  std::string prometheus = "Prometheus: ";
  std::cout << prometheus << "Welcome..." << std::endl;
  std::cout << prometheus << "Generating C++ thermochemistry source code..." << std::endl;
  std::cout << prometheus << "Working with " << m_mech  << " mech..." << std::endl;
  std::cout << prometheus << "Number of species:   "    << m_thermo->nSpecies()     << std::endl;
  std::cout << prometheus << "Number of elements:  "    << m_thermo->nElements()    << std::endl;
  std::cout << prometheus << "Number of reactions: "    << m_kinetics->nReactions() << std::endl;
  
  std::cout << prometheus << "The implementation will ";
  if( m_ooriented == false ) { std::cout << "not "; } 
  std::cout << "be object-oriented..." << std::endl;
  
  std::cout << prometheus << "The implementation will ";
  if( m_templated == false ) { std::cout << "not "; }
  std::cout << "be templated..." << std::endl;
  
}

void CppPrometheus::WriteMech() {

  std::string   filename;
  std::ofstream out;
  out.setf( std::ios::scientific );
  out.precision(6);

  Greetings();

  filename = m_mech + ".H";
  out.open( filename.c_str() );
  WriteDefinitions( out );
  
  ///
  /// Thermo
  ///
  WriteThermo( out );
  
  ///
  /// Kinetics
  ///
  WriteKinetics( out );

  ///
  /// Transport
  ///
  WriteTransport( out );
  
}

///
/// Definitions
///
void CppPrometheus::WriteDefinitions(std::ostream& out) {

  
  out << "... Nothing to see here yet ... " << std::endl;
  
}

///
/// State
///
void CppPrometheus::WriteDensity(std::ostream& out) {
    
}

void CppPrometheus::WritePressure(std::ostream& out) {
    
}

void CppPrometheus::WriteMixMolecularWeight(std::ostream& out) {
    
}

void CppPrometheus::WriteConcentrations(std::ostream& out) {

}

///
/// Thermo
///
void CppPrometheus::WriteMixtureSpecificHeatConstantPressure(std::ostream& out) {
  
}

void CppPrometheus::WriteMixtureSpecificHeatConstantVolume(std::ostream& out) {
  
}

void CppPrometheus::WriteMixtureEnthalpy(std::ostream& out) {
  
}

void CppPrometheus::WriteMixtureInternalEnergy(std::ostream& out) {
  
}

void CppPrometheus::WriteSpeciesSpecificHeats(std::ostream& out) {

  int    ncoeffGuess = 500;
  int    type;
  double minTemp;
  double maxTemp;
  double refPres;
  std::vector<double>  c(15);
  std::vector<double> c9(ncoeffGuess);
  //std::vector<double> c9(34);

  WriteFunctionName( out, "GasThermo", "" );
  
}

void CppPrometheus::WriteSpeciesEnthalpies(std::ostream& out) {
  
}

void CppPrometheus::WriteSpeciesEntropies(std::ostream& out) {
  
}

void CppPrometheus::WriteSpeciesGibbs(std::ostream& out) {
  
}

void CppPrometheus::WriteEquilibriumConstants(std::ostream& out) {
  
}

void CppPrometheus::WriteNewtonTemperature(std::ostream& out) {
  
}

void CppPrometheus::WriteFalloffKinetics(std::ostream& out) {

}

void CppPrometheus::WriteRateCoefficients(std::ostream& out) {

}

void CppPrometheus::WriteArrheniusKinetics(int& rxn,
					  std::shared_ptr<Cantera::Reaction>& reaction,
					  std::ostream& out) {

}

void CppPrometheus::WriteNetRatesOfProgress(std::ostream& out) {

}

void CppPrometheus::WriteNetProductionRates(std::ostream& out) {

}

void CppPrometheus::WriteFunctionName(std::ostream& out,
				   const std::string& className,
				   const std::string& funName) {
  
}
