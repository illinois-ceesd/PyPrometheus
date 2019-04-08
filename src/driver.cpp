#include <sstream>
#include <fstream>
#include <iostream>
#include "opts.h"
#include "Prometheus.h"
#include "cantera/IdealGasMix.h"

using namespace std;

int main(int argc, char** argv) {

  /* iostream format */
  cout.precision(4);
  cout.setf(ios::scientific);

  /* load options */
  opts genOpts;
  if( argc > 1 ) { genOpts.loadOpts( argv[1] ); }
  
  /* ideal gas mixture */
  bool ooriented    = genOpts.ooriented();
  bool templated    = genOpts.templated();
  std::string  mech = genOpts.mech();
  std::string  ctif = "ctis/" + mech + ".xml";
  Cantera::IdealGasMix gas( ctif, "gas" );

  /* run solver */
  Prometheus gen( ooriented, templated, mech, gas );
  gen.writeMech();
  
  return(0);
}
