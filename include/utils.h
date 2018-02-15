#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

std::string fmt(const std::string& var, int i);

void writeCoeffTimesVar(const bool& withSign, double& coeff,
			const std::string& var, std::ostream& out);

std::string increment(const std::string& var, std::vector<int>& indices);

std::string multiply(const std::string& var, std::vector<int>& indices);

std::string wrapString(std::vector<std::string>& vstr);
