#include "Utils.h"

std::string WriteAccess(const std::string& var, int i) {

  std::string out = var + "[" + std::to_string(i) + "]";
  return(out);

}

void WriteCoeffTimesVar(const bool& withSign, double& coeff,
			const std::string& var, std::ostream& out) {

  if(withSign == false) {
    out << coeff << var;
  } else {
    if(coeff > 0.0) {
      out << " + " << coeff << var;
    } else {
      out << " - " << std::fabs(coeff) << var;
    }
  }

}

std::string Increment(const std::string& var, std::vector<int>& indices) {

  std::string out;
  for(int i = 0; i < indices.size(); ++i) {
    if(i != indices.size()-1) {
      out += WriteAccess(var, indices[i]) + " + ";
    } else {
      out += WriteAccess(var, indices[i]);
    }
  }
  return(out);

}

std::string Multiply(const std::string& var, std::vector<int>& indices) {

  std::string out;
  for(int i = 0; i < indices.size(); ++i) {
    if(i != indices.size()-1) {
      out += WriteAccess(var, indices[i]) + " * ";
    } else {
      out += WriteAccess(var, indices[i]);
    }
  }
  return(out);
  
}

std::string WrapString(std::vector<std::string>& vstr) {

  std::string out;
  for(int s = 0; s < vstr.size(); ++s) { out += vstr[s]; }
  return(out);
  
}
