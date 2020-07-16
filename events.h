#ifndef i_events
#define i_events
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "cell.h"

using namespace std;

/// Philippe 2018-04-25
// Elena: events: Define different events to record output
int get_new_ID();
enum {event_born, event_divide, event_die, event_unselected, event_catch_Ag, event_FDC_selected, event_start_contact_TC, event_stop_contact_TC, event_start_signaling_TC, event_stop_signaling_TC, event_recycling, event_become_memory,event_become_plasma, event_become_out, NB_types_events};
class events
{
public:
    events(string outputFolder, string outhustory, string deadhistory);
    vector<string*> storage;

    string outfname;
//    string deadfname;
//    ofstream eventsDead;
    ofstream eventsOut;

    double currentEventsTime;

    string nameEvent(int typeEvent);
    void extendStorage(int ID);
    void recordEvent(B_cell* bc, int typeEvent, double t = -1);
    void eventSetTime(double newTime);
    void writeEvents();
};
#endif

