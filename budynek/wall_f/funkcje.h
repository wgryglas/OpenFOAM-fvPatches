#ifndef funkcje_h
#define funkcje_h
#include <stdio.h>
#include <stdlib.h>
#include "struktury.h"
#define ERROR_ALLOC printf("ALOKACJA NIE UDALA SIE\n")
int translation(struct Wymiary *Wsk_St_wymiary, unsigned int I, unsigned int J, unsigned int K, int typ);
void czytaj_rozmiar(struct Wymiary *Wsk_St_wymiary, struct IOdata *Wsk_St_IOdata);

void inicjalizacja(struct Wymiary *Wsk_St_Wymiary, struct Plyn *Wsk_St_plyn, struct s_Predkosc *Wsk_St_predkosc,struct Cisnienie *Wsk_St_cisnienie,struct IOdata *Wsk_St_IOdata);
#endif
