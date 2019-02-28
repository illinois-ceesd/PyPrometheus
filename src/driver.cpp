#include <sstream>
#include <fstream>
#include <iostream>
#include "opts.h"
#include "chemGen.h"
#include "cantera/IdealGasMix.h"

using namespace std;

int main(int argc, char** argv) {

  /* iostream format */
  cout.precision(4);
  cout.setf(ios::scientific);

  /* load options */
  opts cgopts;
  if( argc > 1 ) { cgopts.loadOpts( argv[1] ); }
  
  /* ideal gas mixture */
  bool ooriented    = cgopts.ooriented();
  bool templated    = cgopts.templated();
  std::string  mech = cgopts.mech();
  std::string  ctif = "ctis/" + mech + ".cti";
  Cantera::IdealGasMix gas( ctif, "gas" );

  /* run solver */
  chemGen gen( ooriented, templated, mech, gas );
  gen.writeMech();
  
  return(0);
}
