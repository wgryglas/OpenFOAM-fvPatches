
#ifndef obliczenia_h
#define obliczenia_h
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struktury.h"
#include "funkcje.h"
void obliczanie_predkosci(short int indeks, unsigned int iter, unsigned int i, unsigned int j, unsigned k, struct Wymiary *Wsk_St_Wymiary, struct n_Predkosc *Wsk_St_tym_predkosc, struct s_Predkosc *Wsk_St_predkosc, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct IOdata *Wsk_St_IOdata);
void obliczanie_cisnienia(unsigned int iter, unsigned int i, unsigned int j, unsigned k, struct Wymiary *Wsk_St_Wymiary, struct s_Predkosc *Wsk_St_predkosc, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn);
#endif
