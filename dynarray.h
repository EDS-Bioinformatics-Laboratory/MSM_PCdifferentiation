#ifndef i_dynarray
#define i_dynarray

// typedef int vartyp;

/*
 * Erzeugt eine Array aus Variablen des Typs vartyp
 * Bugs: ###
 *
 * Gebrauchsanweisung:
 *
 * - Erzeugung durch    dynarray feldname(Groesse,control,fixed);
 *
 * - Groesse gibt die anfaengliche Laenge des Feldes an.
 * - control=0 bedeutet keine Veraenderung der Variablen belegt.
 *           Sie wird dann immer gleich der maximalen Laenge.
 * - control=1 bedeutet, dass nur die explizit beschriebenen Felder
 *           als belegt deklariert sind. Lesen ausserhalb dieses
 *           Bereichs ist nicht moeglich.
 * - fixed=1 bedeutet, dass die Groesse des Feldes nicht variabel ist.
 * - fixed=0 bedeutet, dass ein Schreiben ausserhalb des definierten
 *         Bereichs zu einer automatischen Vergroesserung des Feldes
 *         fuehrt. (Steuerung auch durch setstep!)
 *
 * - Indizierung laeuft von 0 bis Groesse-1
 *
 * - Der Speicher wird beim Verlassen der entsprechenden Routine
 * automatisch geloescht
 *
 * - Die Kenndaten eines erzeugten Arrays sind mit den Funktionen
 * dynarray::  benutzt / lang
 * abrufbar
 *
 * - Das Verhalten des Feldes beim Ueberschreiten der definierten
 * Laenge (dynarray::lang()) wird durch die lokale Variable
 * dynarray::step bestimmt.
 * Diese wird durch dynarray::setstep(wert) gesetzt.
 * Default-Wert ist 10% der Laenge des Feldes.
 * step=0 bedeutet keine Moeglichkeit der Verlaengerung des Feldes
 * step>0 bedeutet automatische Erweiterung des Feldes in
 *        Schritten der Groesse step
 * Die Moeglichkeit der automatischen Felderweiterung erlaubt es,
 * die Felder zunaechst knapp zu definieren und damit Speicher
 * zu sparen. Allerdings kann das haeufige Verlaengern des Feldes
 * zu Geschwindigkeitsverlusten fuehren.
 *
 * - Die Benutzung eines erzeugten Arrays erfolgt durch die
 * Funktionen
 * dynarray:: write, add, read, erase, insert
 * (siehe Deklaration weiter unten)
 *
 */

#include <iostream>
#include <cstdlib>
#include <new>
using namespace std;

const char DYN_DOKU = 0;
// const int DOKU_LIMIT=8200;
const int DOKU_LIMIT = 255000;

template <class vartyp>

class dynarray {
  public:
   dynarray(); // default Kontruktor

   dynarray(long int laenge, char control, char fixed); // Konstruktor
   /* Erzeugt einen Array von LAENGE Feldelementen
    * und initialisiert dieses  */

   // Self-definition
   dynarray(dynarray<vartyp> &tmp);

   ~dynarray();
   /* Freigeben des Feldes */

   char Vergleich(dynarray<vartyp> &a);
   // Liefert 1, falls a.field==field, sonst 0

   int hal;

   long int benutzt();
   /* Liefert die Zahl der benutzten Feldelemente in einem Feld zurueck */

   void delete_all();
   /* Loescht alle Eintraege, indem belegt=0 gesetzt wird */

   long int lang();
   /* Liefert die Laenge des Feldes zurueck */

   void setstep(const int &wert);
   /* Setzt den Wert von step auf wert. wert=0 unterdrueckt
    * moegliche Felderweiterungen. */

   void setcontrol(char wert);
   /* Setzt die Kontrolle des Indexbereichs ueber Kontrolle
    * ausser Kraft (Vorbereitung fuer deleatoren) */

   char write(const long int &wo, const vartyp &was);
   /* Schreibt den Wert WAS an der Stelle WO in das Feld.
    * Liegt die Stelle WO ausserhalb des Arrays, wird er verlaengert.
    * Es ist in der aufrufenden Routine darauf zu achten, dass WO nicht
    * mehr als IND_step oberhalb der Array-Groesse liegt. Der Grundgedanke
    * ist, dass WO maximal eine Stelle oberhalb von belegt liegt. */

   long add(const vartyp &was);
   // Wie write, nur wird immer am Ende des Feldes angehaengt

   vartyp &read(const long int &wo);
   // Liest den Wert an der Stelle wo und gibt ihn zurueck

   char erase(const long int &wo);
   // Loescht den Wert an der Stelle wo und rueckt die Werte dahinter auf

   char erase_jump(const long int &wo);
   // Wie erase nur nicht nachgeruecken sondern letzten Eintrag vorsetzten

   char insert(const long int &wo, const vartyp &was);
   // Setzt den Wert was an der Stelle wo ein und schiebt alle
   // Werte ab wo um eins nach hinten

   void show();
   // Gibt alle Element des Feldes auf dem Bildschirm aus.

   long int find(const vartyp &was);
   long int find(const vartyp &was, long int von, long int bis);
   // Liefert die erste Position im Bereich von-bis, an der was steht
   // Kommt was nicht vor, wird find=-1

   long int howoften(const vartyp &was, long int von, long int bis);
   // Liefert, wieoft was im Bereich von-bis vorkommt

   long int gemeinsam(dynarray<vartyp> &vergl);
   /* Liefert die Zahl der gemeinsamen Elemente in den Feldern eins und zwei.
    * Konzipiert fuer Felder, in denen Elemente nicht doppelt vorkommen.
    * Doppelte Elemente werden in allen Kombinationen mehrfach gezaehlt. */

   char kopie(dynarray<vartyp> &ziel, long int ab, long int bis, long int nach);
   /* Kopiert von dem Feld die Werte ab "ab" bis "bis"
    * in das Feld ziel und beginnt an der Position "nach".
    * Wenn ziel.step entsprechend gesetzt ist, wird dabei
    * das Feld ziel --- falls notwendig --- automatisch erweitert. */

   // Zuweisungs-Operator
   dynarray<vartyp>&operator =(dynarray<vartyp> &x);

   // Index-Operator
   vartyp&operator [](long int i);

  private:
   // Dokumentation
   // char DYN_DOKU;

   // Feld von Zeigern
   vartyp * ptrfield;

   // Feld-Initialisierung
   void def_ptrfield(const long int &laenge);

   // Speicher des Feldes freigeben
   void destruct();

   // Lokal: Zugriff auf die Feldelemente
   // Lesen an Position n
   // vartyp& field(long int n);
   // Schreiben an Position n
   void sfield(const long int n, const vartyp &was);

   // Zahl der belegten Felder
   long int belegt;

   // Maximale Zahl der Felder
   long int max;

   // Laenge des Typs in Byte
   int typelong;

   // Die Zahl der Schritte bei einer Erweiterung der Heap-Groesse
   int step;

   // =1 falls Kontrolle durch belegt
   char Kontrolle;

   // Vergroessert das vorhandene Feld auf die Groesse BISWO
   // Gibt char = 0 zurueck wenn durchgefuehrt sonst char = 1
   char enlarge(const long int &biswo);
};

/*dynarray<vartyp> operator=(const dynarray<vartyp>& x1,
 *             const dynarray<vartyp>& x2);
 * Ort operator-(Ort r1, const Ort& r2);
 * double operator*(const Ort& r1, const Ort& r2);
 * char operator==(const Ort& a, const Ort& b);
 */

// =====================================================
// Implemantation:
// =====================================================

// Konstruktor: ========================================

template <class vartyp>
void dynarray<vartyp>::def_ptrfield(const long int &laenge) {
   // Dokumentation der Speicherbelegung ? ja=1 nein=0
   // DYN_DOKU=0;
   // Laenge des Typs in Byte
   typelong = sizeof(vartyp);

   // Rueckmeldung:
   /// Philippe 2017-12-4
   /*if ((DYN_DOKU == 1) || (typelong * laenge > DOKU_LIMIT)) {
      cout << "Belege Speicher: " << typelong * laenge << " Byte ";
      cout << "fuer " << laenge << " Elemente mit je "
           << typelong << " Byte!\n";
   }*/

   // Initialisiere das Feld
   ptrfield = new vartyp[laenge];
}
template <class vartyp>
dynarray<vartyp>::dynarray() {
// Default Konstruktor
// belegt entsprechend control deklarieren
   belegt = 0;
   Kontrolle = 1;
   max = 100;
   step = 100;
   // Die eigentliche Felddimensionierung:
   def_ptrfield(100);
   //   cout<<"Run default dynarray-constructor ...\n";
}
template <class vartyp>
dynarray<vartyp>::dynarray(long int laenge, char control, char fixed) {
// Konstruktor
// belegt entsprechend control deklarieren
   if (control == 1) { belegt = 0; } else { belegt = laenge; }
   Kontrolle = control;
   max = laenge;
   if (laenge > 30) {
      if (laenge <= 1000) { step = int (laenge / 10); } else { step = 100; }
   } else { step = 3; }
   // Falls Laenge fest, enlarge unterdruecken
   if (fixed == 1) { step = 0; }
   // Die eigentliche Felddimensionierung:
   def_ptrfield(laenge);
}
template <class vartyp>
dynarray<vartyp>::dynarray(dynarray<vartyp> &tmp) {
   // Self-Konstruktor
   // belegt entsprechend control deklarieren
   belegt = tmp.belegt;
   Kontrolle = tmp.Kontrolle;
   max = tmp.max;
   step = tmp.step;
   // Die eigentliche Felddimensionierung:
   def_ptrfield(max);
}
// Destruktor: ================================

template <class vartyp>
void dynarray<vartyp>::destruct() {
   // Rueckmeldung:
   if ((DYN_DOKU == 1) || (typelong * max > DOKU_LIMIT)) {
      //// Philippe 2017-12-4cout << "Gebe Speicher frei: " << typelong * max << " Byte!\n";
      // char c;
      // cin>>c;
   }
   // Diese Zeile ruft den destruktor der vartyp-Klasse explizit auf:
   // for (int i=0; i<max; ++i) ptrfield[i].~vartyp();
//   delete[] ptrfield;   //tem
}
template <class vartyp>
dynarray<vartyp>::~dynarray() {
   // cout<<"in ~dynarray() ...\n";
   destruct();
}
// Feldzugriffe (private): ================================

/*
 * template<class vartyp>
 * vartyp& dynarray<vartyp>::field(long int n) {
 * return ptrfield[n];
 * }*/

template <class vartyp>
void dynarray<vartyp>::sfield(const long int n, const vartyp &was) {
   ptrfield[n] = was;
}
// Verwendung: ====================================

// Vergleich
template <class vartyp>
char dynarray<vartyp>::Vergleich(dynarray<vartyp> &a) {
   long int n;
   char back = 1;
   for (n = 0; n < belegt; n++) {
      if (a.ptrfield[n] != ptrfield[n]) { back = 0; }
   }
   return back;
}
// benutzt
template <class vartyp>
long int dynarray<vartyp>::benutzt() { return belegt; }

// delete_all
template <class vartyp>
void dynarray<vartyp>::delete_all() { belegt = 0; }

// lang
template <class vartyp>
long int dynarray<vartyp>::lang() { return max; }

// setstep
template <class vartyp>
void dynarray<vartyp>::setstep(const int &wert) { step = wert; }

// setcontrol
template <class vartyp>
void dynarray<vartyp>::setcontrol(char wert) {
   Kontrolle = wert;
   belegt = max;
}
// enlarge (private)
template <class vartyp>
char dynarray<vartyp>::enlarge(const long int &biswo) {
   char c;
   dynarray<vartyp> tmp(max,Kontrolle,0);
   long int newmax = biswo + step;
   if (newmax > max) {
      long int n;
      for (n = 0; n < belegt; n++) {
         tmp.sfield(n,ptrfield[n]);
      }
      long int altbelegt = belegt;
      int altstep = step;
      char altKontrolle = Kontrolle;
      // Gib den Speicher von dem alten Feld ptrfield frei:
      destruct();
      // Neue Feld-Konfiguration laden
      // cout<<"in enlarge ... oldmax="<<max<<"; newmax="<<newmax<<"\n";
      // cin>>c;
      max = newmax;
      step = altstep;
      Kontrolle = altKontrolle;
      // Initialisiere das neue Feld der Laenge newmax
      def_ptrfield(newmax);
      // Laden der alten Daten
      for (n = 0; n < altbelegt; n++) {
         sfield(n,tmp.ptrfield[n]);
      }
      if (Kontrolle == 1) { belegt = altbelegt; } else { belegt = max; }
      c = 0;
   } else { c = 1; }
   return c;
}
// write
template <class vartyp>
char dynarray<vartyp>::write(const long int &wo, const vartyp &was) {
   // Kontrolle, ob WO im erlaubten Bereich
   if ((wo >= max) && (step > 0)) { enlarge(wo); }
   if (wo < max) {
      // Ordne was an der Stelle wo zu:
      sfield(wo,was);
      // Anpassung von belegt
      if (wo >= belegt) { belegt = wo + 1; }
      return 0;
   }
   cerr << "dynarray::write Feldindex ausserhalb des erlaubten Bereichs!\n";
   return 1;
}
// add
template <class vartyp>
long dynarray<vartyp>::add(const vartyp &was) {
   if (Kontrolle == 1) {
      // Kontrolle, ob WO im erlaubten Bereich
      if ((belegt == max) && (step > 0)) { enlarge(max); }
      if (belegt < max) {
         // Ordne am Ende des Feldes zu:
         sfield(belegt,was);
         // Anpassung von belegt
         belegt++;
         return belegt - 1;
      } else {
         cout << "dynarray::add Das Feld ist bereits voll!\n";
         return 1;
      }
   }
   cerr << "dynarray::add arbeitet nur fuer Kontrolle=1!\n";
   return -1;
}
// read
template <class vartyp>
vartyp&dynarray<vartyp>::read(const long int &wo) {
   // Falls im erlaubten Bereich
   if (wo < belegt) { return ptrfield[wo]; }
   cerr << "dynarray::read Feldindex ausserhalb des belegten Bereichs!\n";
   return ptrfield[0]; // dummy
}
// erase
template <class vartyp>
char dynarray<vartyp>::erase(const long int &wo) {
   if (wo < belegt) {
      long int n;
      for (n = wo; n < belegt - 1; n++) {
         sfield(n,ptrfield[n + 1]);
      }
      if ((belegt > 0) && (Kontrolle == 1)) { belegt = belegt - 1; }
      return 0;
   }
   cerr << "dynarray::erase Feldindex ausserhalb des belegten Bereichs!\n";
   char c;
   cin >> c;
   return 1;
}
// erase_jump
template <class vartyp>
char dynarray<vartyp>::erase_jump(const long int &wo) {
   if (wo < belegt) {
      sfield(wo,ptrfield[belegt - 1]);
      if ((belegt > 0) && (Kontrolle == 1)) { belegt = belegt - 1; }
      return 0;
   }
   cerr << "dynarray::erase Feldindex ausserhalb des belegten Bereichs!\n";
   char c;
   cin >> c;
   return 1;
}
// insert
template <class vartyp>
char dynarray<vartyp>::insert(const long int &wo, const vartyp &was) {
   // Kontrolle, ob ein Wert noch rein passt.
   if ((belegt == max) && (step > 0)) { enlarge(max); }
   if ((belegt < max) || (Kontrolle == 0)) {
      // Verschieben der Werte hinter WO}
      long int anfang = belegt - 1;
      // Falls keine Kontrolle ist belegt=max. Dann den letzten Wert nicht
      // schieben!
      if (Kontrolle == 0) { anfang--; }
      for (long int n = anfang; n >= wo; n--) {
         sfield(n + 1,ptrfield[n]);
      }
      // Eintragen des neuen Werts
      sfield(wo,was);
      if (Kontrolle == 1) { belegt++; }
      return 0;
   }
   cerr << "dynarray::insert Das Feld ist bereits voll!\n";
   return 1;
}
// show
template <class vartyp>
void dynarray<vartyp>::show() {
   cout << belegt << " field elements:\n";
   for (long n = 0; n < belegt; ++n) {
      cout << ptrfield[n] << "; ";
   }
   cout << "\n";
}
// find
template <class vartyp>
long int dynarray<vartyp>::find(const vartyp &was, long int von, long int bis) {
   long int found = -1;
   if (von == -1) { von = 0; bis = belegt - 1; }
   long int n = von;
   while (n <= bis && found == -1) {
      if (ptrfield[n] == was) { found = n; } else { n++; }
   }
   return found;
}
template <class vartyp>
long int dynarray<vartyp>::find(const vartyp &was) {
   return find(was,-1,-1);
}
// howoften
template <class vartyp>
long int dynarray<vartyp>::howoften(const vartyp &was, long int von, long int bis) {
   long int found = -1;
   long int zahl = 0;
   if (von == -1) { von = 0; bis = belegt - 1; }
   found = find(was,von,bis);
   while (found >= 0) {
      zahl++;
      found++;
      if (found <= bis) { found = find(was,found,bis); } else { found = -1; }
   }
   return zahl;
}
// gemeinsam
template <class vartyp>
long int dynarray<vartyp>::gemeinsam(dynarray<vartyp> &vergl) {
   long int zahl = 0;
   long int n;
   long int m;
   vartyp wert1;
   vartyp wert2;
   for (n = 0; n < belegt; n++) {
      wert1 = read(n);
      for (m = 0; m < vergl.belegt; m++) {
         wert2 = vergl.read(m);
         if (wert1 == wert2) { zahl++; }
      }
   }
   return zahl;
}
// kopie
template <class vartyp>
char dynarray<vartyp>::kopie(dynarray<vartyp> &ziel,
                             long int ab, long int bis, long int nach) {
   // Check der Begrenzungen
   char ok = 0;
   if (ab == -1) { ab = 0; bis = belegt - 1; }
   if ((bis > belegt) || (ab > bis) || (nach > ziel.belegt)
       || ((nach + bis - ab >= ziel.max) && (ziel.step == 0)))
   { ok = 1; } else {
      if (nach + bis - ab >= ziel.max) { ziel.enlarge(nach + bis - ab + 1); }
      long int n;
      for (n = 0; n <= bis - ab; n++) {
         ziel.write(nach + n,ptrfield[ab + n]);
      }
   }
   return ok;
}
template <class vartyp>
dynarray<vartyp>&dynarray<vartyp>::operator =(dynarray<vartyp> &x) {
   int n;

   // Altes Feld loeschen
   destruct();
   // Eckdaten von x uebernehmen
   max = x.max;
   step = x.step;
   belegt = x.belegt;
   Kontrolle = x.Kontrolle;
   // Initialisiere das neue Feld der Laenge newmax
   def_ptrfield(max);
   // Daten uebertragen
   for (n = 0; n < max; n++) {
      sfield(n,x.ptrfield[n]);
   }

   return *this;
}
template <class vartyp>
vartyp&dynarray<vartyp>::operator [](long int i) {
   if ((i >= 0) && (i < belegt)) { return ptrfield[i]; }

    cerr << "dynarray::[" << i << "] Index nicht im erlaubten Bereich (0.." << belegt - 1 << ")!\n";
   exit(1);
   return ptrfield[0];
}
#endif
