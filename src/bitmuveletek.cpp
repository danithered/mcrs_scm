/*Ezek a fuggvenyek bitmuveletek csinalnak:
 * 	beolvasnak megadott bitet
 * 	modositanak megadott bitet
 * 	kiirnak egy szamot kettes szamrendszerben
*/

#include "bitmuveletek.h"

using namespace std;

//INT SZAMOK
int olvas (int ertek, int bit) {
	/* kiolvassa egy int megadott bitenek erteket
	 * 
	 * ertek: az olvasott int
	 * bit: hanyadik bitet olvassa
	 */
	return((ertek & ( 1 << bit )) >> bit);
}

int ir (int ertek, int bit, int mit) {
  /*int szamokat szerkeszt binaris bit szinten
   * 
   * ertek: a szerkesztett int
   * bit: hanyadik bitet szerkeszti. FIGYELEM: a szamozas 0-tol kezdodik!!!
   * mit: mit irjon be a megadott bitre. Lehetseges ertekek:
   * 		0: 0
   * 		1: 1
   * 		-1: torli az egeszet, nullazza
   * 		minden mas: visszadja az eredeti erteket
  */
	switch (mit) {
		case 0:
			ertek = (~(1 << bit)) & ertek;
			break;
		case 1:
			ertek = (1 << bit) | ertek;
			break;
		case -1:
			return(0);
	}
	return(ertek);
}

int kiirbit (int ertek, int hatar) {
  /* egy intet kiir bitek szerint. Egyreszt kiirja consolera, masreszt visszadja intbe
   * 
   * ertek: az olvasott int
   * hatar: hanyadik bit-ig olvassa. FIGYELEM: a szamozas 0-tol kezdodik!!!
   * 
   * MEGJEGYZES: nem neztem meg meg, hogy mi van akkor, ha tul magas hatarszamot adok meg neki, lehet osszedol.
   */
	int bitszam=0, bit, intbe=0;
	
	//printf("\n");
	for (bitszam=(hatar-1); bitszam >= 0; bitszam--) {
		bit= olvas(ertek, bitszam);
		//printf("%d\t", bit);
		intbe+= bit*pow(10, bitszam);
		
	}
	//printf("\n");
	return(intbe);
}

//POINTEREK
int olvasP (int *ertek, int bit) {
	/* kiolvassa egy int pointer megadott bitenek erteket
	 * 
	 * ertek: az olvasott int pointer
	 * bit: hanyadik bitet olvassa
	 */
	
	return((*ertek & ( 1 << bit )) >> bit);
}

int irP (int *ertek, int bit, int mit) {
  /* int pointer szamokat szerkeszt binaris bit szinten
   * 
   * ertek: a szerkesztett int
   * bit: hanyadik bitet szerkeszti. FIGYELEM: a szamozas 0-tol kezdodik!!!
   * mit: mit irjon be a megadott bitre. Lehetseges ertekek:
   * 		0: 0
   * 		1: 1
   * 		-1: torli az egeszet, nullazza
   * 		minden mas: visszadja az eredeti erteket
   * 
   * Kimenetek:
   * 	0: lefutott, 0-t irt be
   * 	1: lefutott, 1-t irt be
   * 	-1: lefutott, nullazott
   * 	-2: nem tortent semmi
  */
	switch (mit) {
		case 0:
			*ertek = (~(1 << bit)) & *ertek;
			return (0);
		case 1:
			*ertek = (1 << bit) | *ertek;
			return (1);
		case -1:
			*ertek = 0;
			return(-1);
	}
	return(-2);
}

int kiirbitP (int *ertek, int hatar) {
  /* egy int pointert kiir bitek szerint. Egyreszt kiirja consolera, masreszt visszadja intbe
   * 
   * ertek: az olvasott int
   * hatar: hanyadik bit-ig olvassa. FIGYELEM: a szamozas 0-tol kezdodik!!!
   * 
   * MEGJEGYZES: nem neztem meg meg, hogy mi van akkor, ha tul magas hatarszamot adok meg neki, lehet osszedol.
   */
	int bitszam=0, bit, intbe=0;
	
	cout << endl;
	for (bitszam=(hatar-1); bitszam >= 0; bitszam--) {
		bit= olvas(*ertek, bitszam);
		cout << bit << "\t";
		intbe+= bit*pow(10, bitszam);
		
	}
	cout << endl;
	return(intbe);
}
