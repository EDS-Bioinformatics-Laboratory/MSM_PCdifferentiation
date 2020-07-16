#include "events.h"
//#include "network.h"
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;


/// Philippe 2018-04-24


events::events(string outputFolder, string outhustory, string deadhistory)
{
    outfname = string(outputFolder + outhustory); //Elena: events: define output path using output folder
    currentEventsTime = 0;
}
string events::nameEvent(int typeEvent){
    switch(typeEvent){
    case event_unselected :{return string("unselected");}
    case event_catch_Ag:{return string("catch_FDC");}
    case event_FDC_selected: {return string("FDC_selected");}
    case event_die:{return string("die");}
    case event_start_contact_TC:{return string("start_contact_TC");}
    case event_stop_contact_TC: {return string("stop_contact_TC");}
    case event_start_signaling_TC: {return string("start_signaling_TC");}
    case event_stop_signaling_TC:{return string("stop_signaling_TC");}
    case event_recycling:{return string("recycling");}
    case event_divide: {return string("divide");}
    case event_born:{return string("born");}
    case event_become_memory:{return string("become_memory");} //Indicate when memory cells exit simulation
    case event_become_plasma:{return string("become_plasma");} //Indicate when plasma cells exit simulation
    case event_become_out:{return string("become_out");} //Indicate when plasma and memory cells exit simulation
    case NB_types_events:{return string("Not_an_event");}
    }
    return string("Not_an_event");
}

void events::extendStorage(int ID)
{
    static bool firstTime = true; 

    if(firstTime)
    {
        storage.resize(1000000, NULL);
        firstTime = false;

        eventsOut.open( outfname , ios::out);
        if (!eventsOut) cerr << "ERROR! eventsOutput did not open";

        stringstream res;
        res << "time\tID\tMID\tcellType\tcellState\tBCL6\tBLIMP1\tIRF4\tEvent\tAffinity\tTotalNdivisions\tnMutationsfromGermline\tIamAgHigh"<< endl;

        eventsOut << res.str();
    }

//    if(ID > (int) storage.size())
//    {
//        storage.resize(ID+1,NULL);
//    }

   if(!storage[ID]) storage[ID] = new string();

}

void events::writeEvents()
{

    if(eventsOut) {
        eventsOut.close();
        cerr<<"Event output files written to: "<< outfname <<" ; "<< endl;
    }
    else {cerr<<"Error: writeEvents() no eventOutput file!"<<endl;}

}

void events::recordEvent(B_cell* bc, int typeEvent, double t)
{
    if(t < 0) t = currentEventsTime;

    extendStorage(bc->ID);

    stringstream res2;

    res2 << t << "\t" << bc->ID << "\t" <<bc->MID <<"\t" <<bc->cell_type<<"\t" <<bc->cell_state<<"\t" << bc->BCL6<< "\t" << bc->BLIMP1<< "\t" << bc->IRF4<<"\t"<<nameEvent(typeEvent) << "\t"<< bc->MyAffinity<< "\t" << bc->total_number_of_divisions<< "\t" << bc->myBCR.nMutFromGermline<< "\t"<< bc->IamHighAg<<"\t" <<  endl;
    
    //    storage[bc->ID]->append(res2.str());

    //    eventsOut  << *(storage[bc->ID]) << endl;
        
         eventsOut  << res2.str() << endl;
}

void events::eventSetTime(double t)
{
    currentEventsTime = t;
}
