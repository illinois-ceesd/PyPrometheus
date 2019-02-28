#include "opts.h"

void opts::loadOpts(const char* inputFile) {

  std::ifstream inp;
  std::string   filename = inputFile;
  
  m_filename = filename;
  std::cout << std::endl;
  std::cout << " ... Input file provided ..." << std::endl;
  std::cout << " ... Loading options from " << m_filename << " ... " << std::endl;
  inp.open(m_filename.c_str());
  
  while(!inp.eof()) {

    std::string optName;
    std::string optVal;
    std::string f;
    inp >> optName >> optVal;

    if(optName.length() > 0) {
      
      f = optName.at(0);
      if(f != "#" && f.length() > 0) {
      
	setOpt(optName, optVal);
	
      }

    }
    
  }

  inp.close();

  std::cout << " ... Options loaded successfully ...\n" << std::endl;

}

void opts::setOpt(std::string& optName, std::string& optVal) {

  if(optName == "mech") {

    setMech(optVal);

  } else if(optName == "ooriented") {

    setParadigm(optVal);
    
  } else if(optName == "templated") {

    setBaseType(optVal);
        
  } else {
    
    std::cout << "\t Option " << optName << " not implemented..." << std::endl;
    
  }
  
}

void opts::setMech(std::string& optVal) {

  m_mech = optVal;

}

void opts::setParadigm(std::string& optVal) {

  if( optVal == "yes" ) {
    m_class = true;
  } else if( optVal == "no" ) {
    m_class = false;
  } else {
    std::cout << "Sorry, " << optVal << " is not known..." << std::endl;
  }

}

void opts::setBaseType(std::string& optVal) {

  if( optVal == "yes" ) {
    m_temp = true;
  } else if( optVal == "no" ) {
    m_temp = false;
  } else {
    std::cout << "Sorry, " << optVal << " is not known..." << std::endl;
  }
  
}
