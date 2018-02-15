#include "utils.h"

std::string fmt(const std::string& var, int i) {

  std::string out = var + "[" + std::to_string(i) + "]";
  return(out);

}

void writeCoeffTimesVar(const bool& withSign, double& coeff,
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

std::string increment(const std::string& var, std::vector<int>& indices) {

  std::string out;
  for(int i = 0; i < indices.size(); ++i) {
    if(i != indices.size()-1) {
      out += fmt(var, indices[i]) + " + ";
    } else {
      out += fmt(var, indices[i]);
    }
  }
  return(out);

}

std::string multiply(const std::string& var, std::vector<int>& indices) {

  std::string out;
  for(int i = 0; i < indices.size(); ++i) {
    if(i != indices.size()-1) {
      out += fmt(var, indices[i]) + " * ";
    } else {
      out += fmt(var, indices[i]);
    }
  }
  return(out);
  
}

std::string wrapString(std::vector<std::string>& vstr) {

  std::string out;
  for(int s = 0; s < vstr.size(); ++s) { out += vstr[s]; }
  return(out);
  
}
