#include "output.h"
#include <set>
#include <numeric>
#include "mafalda.h"
#include "cell.h"
using namespace std;

// ministats
vector<ministat> Bcell_counts;  // number of Bcells
ministat Plasma_counts;

// Affinities
// This is for recording (average) affinities,number of divisions and mutations of B cells versus time, for plasma
// cells there is a different file of records
ministat Bcell_affinity;
ministat Bcell_BLIMP1;
ministat Bcell_IRF4;
ministat Bcell_BCL6;
ministat Bcell_nMut;
ministat Recycling_CC_nDiv;

output::output(string _parfname) : parfname(_parfname) { initialize_fileds(); }
output::~output() {}

string output::currentDateTime() {
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%y-%m-%d_%H-%M-%S", &tstruct);
  return buf;
}



void output::initialize_fileds() {
  // 0-9 -> CC states (7) + CB (2)
  for (int i = 0; i < 8; i++) {
    Bcell_counts.push_back(ministat());  // number of Bcells
  }
}

void output::clear_fileds() {
  for (int i = 0; i < 8; i++) {
    Bcell_counts.clear();  //.at(i).clear_ministat();   //number of Bcells
  }
  for (int i = 0; i < 8; i++) {
    Bcell_counts.push_back(ministat());  // number of Bcells
  }

    // Affinities
    Bcell_affinity.clear_ministat();
    Bcell_BLIMP1.clear_ministat();
    Bcell_BCL6.clear_ministat();
    Bcell_IRF4.clear_ministat();
    Bcell_nMut.clear_ministat();
    Recycling_CC_nDiv.clear_ministat();
}

// Take fields from simulation into master observer variable to create file.
void output::record_output_time_step(double currentTime, simulation &currentSim,
                                     parameters &p) {
  // This function records data about B cells population and affinity versus
  // time, this does not include Plasma cells.

  // Bcell data

  for (unsigned int i = 0; i < currentSim.ListB_cell.size(); i++) {
    B_cell *Bcell = currentSim.ListB_cell.at(i);
    // Check integrity of the list

    if(Bcell->cell_type > 2)//Elena: Make sure you only record B cells!
    {
        cerr << "Error, wrong cell in BC list, Type: " <<Bcell->cell_type <<"; State: "<< Bcell->cell_state << endl;
        exit(1);
    }
    if (Bcell->cell_state > 7)
    {
        cerr << "Error, wrong cell in BC list, Type: " <<Bcell->cell_type <<"; State: "<< Bcell->cell_state << endl;
        exit(1);
    }

      Bcell_counts[Bcell->cell_state].add(1);
      
    // Bcell Affinity
    Bcell->setMyAffinity(p);
    Bcell_affinity.add(Bcell->MyAffinity);
    Bcell_BLIMP1.add(Bcell->BLIMP1);
    Bcell_BCL6.add(Bcell->BCL6);
    Bcell_IRF4.add(Bcell->IRF4);
    Bcell_nMut.add(Bcell->myBCR.nMutFromGermline);
    Recycling_CC_nDiv.add(Bcell->Recycling_divisions); //Elena: lymphoma:
  }

  FILE *Bcell_time_data;
  string folder1 = outputFolder + "/Bcell_time.csv";
  char *s1 = const_cast<char *>(folder1.c_str());
  Bcell_time_data = fopen(s1, "a");
  if (Bcell_time_data == NULL){cerr<<"Output::Error file: "<<folder1<<" empty;\n"<<endl;} //Elena: Manage possible errors.
  static bool tmp = true;
  if (tmp) {
    tmp = false;
    fprintf(Bcell_time_data,"%s",
            "time,founder,unselected,contact_FDC,FDC_selected,contact_TC,Selected_by_TC,recycled,apoptosis,nDiv_TCselected,std_nDiv_TCselected,nMut,std_nMut,affinity,std_affinity,BLIMP1,std_BLIMP1,BCL6,std_BCL6,IRF4,std_IRF4\n");
  }

  fprintf(Bcell_time_data, "%f,", currentTime);
  for (int i = 0; i < 8; i++) {
    fprintf(Bcell_time_data, "%f,", Bcell_counts[i].sum);
  }
    fprintf(Bcell_time_data, "%.16G,%.16G,", Recycling_CC_nDiv.average(),
            Recycling_CC_nDiv.stddev());
    fprintf(Bcell_time_data, "%.16G,%.16G,", Bcell_nMut.average(),
             Bcell_nMut.stddev());
    fprintf(Bcell_time_data, "%.16G,%.16G,", Bcell_affinity.average(),
            Bcell_affinity.stddev());
     fprintf(Bcell_time_data, "%.16G,%.16G,", Bcell_BLIMP1.average(),
              Bcell_BLIMP1.stddev());
    fprintf(Bcell_time_data, "%.16G,%.16G,", Bcell_BCL6.average(),
             Bcell_BCL6.stddev());
    fprintf(Bcell_time_data, "%.16G,%.16G\n", Bcell_IRF4.average(),
             Bcell_IRF4.stddev());
    fclose(Bcell_time_data);  //#Recheck take care of bins in gle file

  clear_fileds();
}

void output::write_event(cell *Cellx, stringstream &sim_output) {
  sim_output << Cellx->event.str() << endl;
}

void output::write_event_2file(stringstream &sim_output) {
  FILE *event_data;
  string folder1 = outputFolder + "/event_data.csv";
  event_data = fopen(folder1.c_str(), "a");

  static bool tmp = false;

  if (not(tmp)) {
    fprintf(event_data, "%s",
            "ID,Born_time,MID,States,Affinity,N_of_Ags,N_of_divisions,N_of_"
            "Mutations,delta_aff,FDC_interaction_nums,FDC_interaction_time_avg,"
            "TC_interaction_time,TC_signaling_time,FDC_selected,Selected_by_TC,"
            "Death_time\n");
    tmp = true;
  }
  //Elena: to record signals add two columns: Interacting with fdc, Interacting with TC. add1 if yes and 0 if no
  //Do not record every time step but only during interaction (see event.cpp)
  //Problem I see with this is

  fprintf(event_data, "%s", sim_output.str().c_str());
  fclose(event_data);  //#Recheck take care of bins in gle file
}
// for B cells

//cell states: 0-founder,
//1-unselected,
//2-contact_FDC,
//3-FDC_selected,
//4-contact_TC,
//5-TC_selected,
//6-recycled,
//7-apoptosis,
//8-TC_free,
//9-TC_connected,
//10-Plasma_Out,
//11-Plasma_in_GC,
//12-cell_state_counter
void output::close_event(B_cell *Cellx, stringstream &sim_output, double time) {
  Cellx->event << Cellx->cell_state << "," << Cellx->MyAffinity << ","
               << Cellx->retained_Ag << "," << Cellx->total_number_of_divisions
               << "," << Cellx->myBCR.nMutFromGermline << ","
               << Cellx->delta_Affinity << ",";
  Cellx->event << Cellx->nFDCcontacts << ",";

  if (Cellx->nFDCcontacts == 0) {
    Cellx->event << Cellx->fdc_interaction_time_history << ",";
  } else {
    Cellx->event << double(Cellx->fdc_interaction_time_history /
                           double(Cellx->nFDCcontacts))
                 << ",";
  }

  Cellx->event << Cellx->Tc_interaction_history.first << ","
               << Cellx->Tc_interaction_history.second << ",";

  Cellx->event << Cellx->Selected_by_FDC << "," << Cellx->Selected_by_TC << "," << time;
}

void output::Plasma_output(double currentTime, simulation &currentSim,
                           parameters &p) {
  // Danial: This function only writes down the data of plasma cells at the end
  // of the simulation.
  /*The order is
   1-Time of production(differetiation)
   2-ID
   3-ID of mother B-cell
   4-Total number of divisions
   5-Total amount of Ag
   6-Affinity
   7-Total number of mutations
   */

  FILE *Plasma_cells_data;
  string folder1 = outputFolder + "/Plasma_cells_data.csv";

  //    char *ss1 = const_cast<char*>(folder1.c_str());;

  Plasma_cells_data = fopen(folder1.c_str(), "a");

  //Elena: It is easyer to know from data files what were actually seeing if we use similar names to fields. What do you think?
  fprintf(Plasma_cells_data, "%s",
          "ID,Born_time,MID,GMID,Affinity,Zpos,Zpolarity,retained_Ag,BLIMP1,BCL6,IRF4,N_of_divisions,N_of_Mutations\n");

  // Plasma data
  for (int j = 0; j < currentSim.ListP_cell.size(); j++) {
    Plasma_cell *Plasma = currentSim.ListP_cell.at(j);

    fprintf(Plasma_cells_data, "%d,%f,%d,%d,%.16G,%f,%f,%f,%f,%f,%f,%d,%d\n",
            Plasma->ID, Plasma->birth_time, Plasma->MID,Plasma->GMID, Plasma->MyAffinity,
            Plasma->position.Z, Plasma->polarity.Z, Plasma->retained_Ag, Plasma->BLIMP1,Plasma->BCL6,Plasma->IRF4, Plasma->total_number_of_divisions,
            Plasma->myBCR.nMutFromGermline);
  }
  fclose(Plasma_cells_data);  //#Recheck take care of bins in gle file
}

void output::Memory_output(double currentTime, simulation &currentSim,
                           parameters &p) {
  // Danial: This function only writes down the data of plasma cells at the end
  // of the simulation.
  /*The order is
   1-Time of production(differetiation)
   2-ID
   3-ID of mother B-cell
   4-Total number of divisions
   5-Total amount of Ag
   6-Affinity
   7-Total number of mutations
   */

  FILE *Memory_cells_data;
  string folder2 = outputFolder + "/Memory_cells_data.csv";

  //    char *ss1 = const_cast<char*>(folder1.c_str());;

  Memory_cells_data = fopen(folder2.c_str(), "a");

  //Elena: It is easyer to know from data files what were actually seeing if we use similar names to fields. What do you think?
  fprintf(Memory_cells_data, "%s",
          "ID,Born_time,MID,GMID,Affinity,Zpos,Zpolarity,retained_Ag,BLIMP1,BCL6,IRF4,N_of_divisions,N_of_Mutations\n");

  // Plasma data
  for (int j = 0; j < currentSim.ListM_cell.size(); j++) {
    Memory_cell *Memory = currentSim.ListM_cell.at(j);

    fprintf(Memory_cells_data, "%d,%f,%d,%d,%.16G,%f,%f,%f,%f,%f,%f,%d,%d\n",
            Memory->ID, Memory->birth_time, Memory->MID,Memory->GMID, Memory->MyAffinity,
            Memory->position.Z, Memory->polarity.Z, Memory->retained_Ag, Memory->BLIMP1,Memory->BCL6,Memory->IRF4, Memory->total_number_of_divisions,
            Memory->myBCR.nMutFromGermline);
  }
  fclose(Memory_cells_data);  //#Recheck take care of bins in gle file
}

