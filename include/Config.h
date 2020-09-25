#include <iostream>
#include <fstream>
#include <string>

class Config {

public:
 Config() :
  m_filename(""),
  m_language("C++"),
  m_thermoXY("HP"),
  m_mech("sanDiego"),
  m_class(true),
  m_temp(true) { };
  
  int LoadConfig(const char* inputFile);

  int SetConfigOption(std::string& optName, std::string& optVal);

  int SetMech(std::string& optVal);

  int SetLanguage(std::string& optVal);

  int SetFixedThermoVars(std::string& optVal);
  
  int SetParadigm(std::string& optVal);
  
  int SetBaseType(std::string& optVal);
  
  std::string mech() { return m_mech; }

  std::string language() { return m_language; }

  std::string fixedThermoVars() { return m_thermoXY; } 
  
  bool ooriented() { return m_class; }
  
  bool templated() { return m_temp; }

  void ErrorMessage(const std::string& optName,
		    const std::string& optVal) {
    
    std::cout << "Prometheus:Config: Sorry, "
	      << optName << " cannot take value "
	      << optVal  << " ... " << std::endl;
    
  }
  
protected:
  bool m_class;
  bool m_temp;  
  std::string m_language;
  std::string m_thermoXY;
  std::string m_filename;
  std::string m_mech;

};
