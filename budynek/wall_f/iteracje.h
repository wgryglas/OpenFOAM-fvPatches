#ifndef iteracje_h
#define iteracje_h
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "src/solvers.h"
#include "f_struktury.h"
#include "obliczenia.h"
#include "funkcje.h"

void iteracja(struct Wymiary *Wsk_St_Wymiary, struct s_Predkosc *Wsk_St_s_predkosc, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct IOdata *Wsk_St_IOdata);
void iteracja_predkosci(struct Wymiary *Wsk_St_Wymiary, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct s_Predkosc *Wsk_St_s_predkosc, struct IOdata *Wsk_St_IOdata);
void iteracja_cisnienia( struct Wymiary *Wsk_St_Wymiary, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct s_Predkosc *Wsk_St_s_predkosc);
#endif
