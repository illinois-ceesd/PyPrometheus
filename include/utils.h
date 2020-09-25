#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

std::string WriteAccess(const std::string& var, int i);

void WriteCoeffTimesVar(const bool& withSign, double& coeff,
			const std::string& var, std::ostream& out);

std::string Increment(const std::string& var, std::vector<int>& indices);

std::string Multiply(const std::string& var, std::vector<int>& indices);

std::string WrapString(std::vector<std::string>& vstr);
