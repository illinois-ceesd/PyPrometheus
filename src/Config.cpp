#include "Config.h"

int Config::LoadConfig(const char* inputFile) {

  int iret = 0;
  std::ifstream inp;
  std::string   filename = inputFile;
  
  m_filename = filename;
  std::cout << std::endl;
  std::cout << "Prometheus:Config: Input file provided ..." << std::endl;
  std::cout << "Prometheus:Config: Loading config options from "
	    << m_filename << " ... " << std::endl;
  inp.open(m_filename.c_str());
  
  while(!inp.eof()) {

    std::string optName;
    std::string optVal;
    std::string f;
    inp >> optName >> optVal;

    if(optName.length() > 0) {
      
      f = optName.at(0);
      if(f != "#" && f.length() > 0) {
      
	iret = SetConfigOption(optName, optVal);	
	
      }

    }

    if( iret ) { break; }
    
  }

  inp.close();

  if( iret == 0 ) {
    std::cout << "Prometheus:Config: Config options loaded successfully ...\n" << std::endl;
  } else {
    std::cout << "Prometheus:Config: Config options could not be loaded ...\n" << std::endl;
  }

  return (iret);  

}

int Config::SetConfigOption(std::string& optName, std::string& optVal) {

  int iret = 0;
  
  if( optName == "mech" ) {

    iret = SetMech( optVal );

  } else if( optName == "language" ) {

    iret = SetLanguage( optVal );

  } else if( optName == "thermoXY" ) {

    iret = SetFixedThermoVars( optVal );
    
  } else if( optName == "ooriented" ) {

    iret = SetParadigm( optVal );
    
  } else if( optName == "templated" ) {

    iret = SetBaseType( optVal );
        
  } else {
    
    std::cout << "Prometheus:Config: Option " << optName << " not implemented..." << std::endl;
    iret = 1;
    
  }

  return(iret);
  
}

int Config::SetMech(std::string& optVal) {

  m_mech = optVal;
  return(0);

}

int Config::SetLanguage(std::string& optVal) {

  m_language = optVal;
  if( m_language != "C++" && m_language != "Python" ) {
    ErrorMessage( "language", m_language );
    return(1);
  }
  return(0);
  
}

int Config::SetFixedThermoVars(std::string& optVal) {

  m_thermoXY = optVal;
  if( m_thermoXY != "HP" && m_thermoXY != "UV" ) {
    ErrorMessage( "thermoXY", m_thermoXY );
    return(1);
  }
  return(0);
  
}

int Config::SetParadigm(std::string& optVal) {
  
  if( optVal == "yes" ) {
    m_class = true;
  } else if( optVal == "no" ) {
    m_class = false;
  } else {
    std::cout << "Sorry, " << optVal << " is not known..." << std::endl;
    ErrorMessage( "ooriented", optVal );
    return(1);
  }
  return(0);

}

int Config::SetBaseType(std::string& optVal) {

  if( optVal == "yes" ) {
    m_temp = true;
  } else if( optVal == "no" ) {
    m_temp = false;
  } else {
    ErrorMessage( "templated", optVal );
    return(1);
    std::cout << "Sorry, " << optVal << " is not known..." << std::endl;
  }
  return(0);
  
}
