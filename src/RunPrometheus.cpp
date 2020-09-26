#include <sstream>
#include <fstream>
#include <iostream>
#include "CppPrometheus.h"
#include "PyPrometheus.h"

using namespace std;

Prometheus* PrometheusFactory(Config& config,
			      std::shared_ptr<Cantera::Solution> solution) {
  
  if( config.language() == "C++" ) {
    return new CppPrometheus( config, solution );    
  } else if( config.language() == "Python" ) {
    return new PyPrometheus( config, solution );    
  }
  
};

int main(int argc, char** argv) {

  /* iostream format */
  cout.precision(4);
  cout.setf(ios::scientific);

  /* load options */
  Config config;
  if( argc > 1 ) {
    if( config.LoadConfig( argv[1] ) ) {
      std::cout << "Prometheus: Config failed, cannot continue ... " << std::endl;
      return(1);
    }
  }

  std::string mech  = config.mech();
  std::string ctif  = "../inputs/" + mech + ".cti";
  std::shared_ptr<Cantera::Solution> solution = Cantera::newSolution( ctif, "gas", "None" );

  /* run solver */
  Prometheus* prometheus = PrometheusFactory( config, solution );
  prometheus->WriteMech();
  
  return(0);
  
};
