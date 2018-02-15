#include <sstream>
#include <fstream>
#include <iostream>
#include "chemGen.h"
#include "cantera/IdealGasMix.h"

using namespace std;

int main() {

  /* iostream format */
  cout.precision(4);
  cout.setf(ios::scientific);

  /* ideal gas mixture */
  bool templated = false;
  std::string  mech = "gri30";
  std::string  ctif = "ctis/"+mech+".cti";
  std::string  outf = mech+".cpp";
  std::ofstream out(outf.c_str());
  Cantera::IdealGasMix gas(ctif,"gas");

  /* run solver */
  chemGen gen(templated, gas);
  gen.writeMech(out);
  
  return(0);
}
