#include <iostream>
#include <fstream>
#include <string>

class opts {

public:
 opts() :
  m_filename(""),
  m_mech("sanDiego"),
  m_class(true),
  m_temp(true) { };
  
  void loadOpts(const char* inputFile);

  void setOpt(std::string& optName, std::string& optVal);

  void setMech(std::string& optVal);

  void setParadigm(std::string& optVal);
  
  void setBaseType(std::string& optVal);

  std::string mech() { return m_mech; }

  bool ooriented() { return m_class; }
  
  bool templated() { return m_temp; }
  
protected:
  bool m_class;
  bool m_temp;  
  std::string m_filename;
  std::string m_mech;

};
