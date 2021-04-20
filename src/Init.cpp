#include "../include/init.h"


using namespace std;







MasterClass_analysis::MasterClass_analysis(string Path) {  //default constructor


  ifstream ReadInput(Path);

  //Read Header
  string Header;
  getline(ReadInput,Header);
  if(Header != "# Input file for the analysis of QCD+QED observables") crash("Input file does not have the expected header");

  string ReadLine;
  do {getline(ReadInput,ReadLine);} while(ReadLine=="" || ReadInput.eof());
  if(ReadLine != "#Analysis_Mode") crash("Error in input file. Expected: #Analysis_Mode, obtained: "+ReadLine+"|");
  ReadInput >> Analysis_Mode;


  
  if(Analysis_Mode=="Meson_masses") {
    do {getline(ReadInput, ReadLine);} while(ReadLine=="" || ReadInput.eof());
    if(ReadLine != "#Init_description_meson_masses") crash("Error in input file. Expected: #Init_description_meson_masses, obtained: "+ReadLine);
    ReadInput>>ReadLine;
    if(ReadLine != "#Meson_to_Analyze:") crash("Error in input file. Expected: #Meson_to_Analyze, obtained: "+ReadLine);
    ReadInput >> Meson_to_analyze;
    ReadInput >> ReadLine;
    if(ReadLine != "#CURRENT_TYPE:") crash ("Error in input file. Expected: #CURRENT_TYPE, obtained: "+ReadLine);
    ReadInput >> CURRENT_TYPE;
    ReadInput >> ReadLine;
    if(ReadLine != "#Include_disconnected?:") crash("Error in input file. Expected: #Include_disconnected?, obtained: "+ReadLine);
    ReadInput >> IncludeDisconnected;
  }

   
  else if(Analysis_Mode=="Meson_masses_twisted") {
    do {getline(ReadInput, ReadLine);} while(ReadLine=="" || ReadInput.eof());
    if(ReadLine != "#Init_description_meson_masses") crash("Error in input file. Expected: #Init_description_meson_masses, obtained: "+ReadLine);
    ReadInput>>ReadLine;
    if(ReadLine != "#Meson_to_Analyze:") crash("Error in input file. Expected: #Meson_to_Analyze, obtained: "+ReadLine);
    ReadInput >> Meson_to_analyze;
    ReadInput >> ReadLine;
    if(ReadLine != "#CURRENT_TYPE:") crash ("Error in input file. Expected: #CURRENT_TYPE, obtained: "+ReadLine);
    ReadInput >> CURRENT_TYPE;
    ReadInput >> ReadLine;
    if(ReadLine != "#Include_disconnected?:") crash("Error in input file. Expected: #Include_disconnected?, obtained: "+ReadLine);
    ReadInput >> IncludeDisconnected;
  }
    
  else if(Analysis_Mode=="Form_factors") {
    //do nothing  
  }

  else if(Analysis_Mode == "Axion_l7") {
    //do nothing
  }

  else if(Analysis_Mode == "Neutral_pi_TM") {
    //do nothing

  }

  else if(Analysis_Mode =="Eta_TM") {
    //do nothing

  }

  else if(Analysis_Mode=="VMD") {
    //do nothing
  }
  
  
  else crash("Analysis_Mode: "+Analysis_Mode+" not found");
  
  ReadInput.close();
  Analysis_manager();
 
  return;

}


void MasterClass_analysis::Analysis_manager() {

  if(Analysis_Mode=="Form_factors") Compute_form_factors();

  if(Analysis_Mode=="Meson_masses") {
    if(Meson_to_analyze =="PI") Pion_mass_analysis(this->CURRENT_TYPE, this->IncludeDisconnected);
  }

  if(Analysis_Mode=="Meson_masses_twisted") {
    if(Meson_to_analyze =="PI") Pion_mass_analysis_twisted(this->CURRENT_TYPE, this->IncludeDisconnected);
  }

  if(Analysis_Mode == "Axion_l7") {
    Axion_l7_analysis();
  }

  if(Analysis_Mode == "Neutral_pi_TM") {
    Neutral_pi_TM();
  }

  if(Analysis_Mode == "Eta_TM") {
    Eta_TM();
  }

  if(Analysis_Mode == "VMD") {
    Regge();
  }
  return;


}
