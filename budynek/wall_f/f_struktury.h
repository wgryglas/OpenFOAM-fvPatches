#ifndef f_struktury_h
#define f_struktury_h
#include <stdio.h>
#include <stdlib.h>
#include "struktury.h"

#define ERROR_ALLOC printf("ALOKACJA NIE UDALA SIE\n")
struct n_Predkosc *alokacja_n_p(const double *dini, unsigned int *ilosc_wekt_IJK);
struct s_Predkosc *alokacja_s_p( struct Wymiary *Wsk_St_wymiary);
struct Cisnienie *alokacja_c(const struct Wymiary *Wsk_St_wymiary);
struct Plyn *alokacja_plyn(const struct Wymiary *Wsk_St_wymiary);
void del_st_n_Predkosc(struct n_Predkosc *Wsk_St_predkosc);
void del_st_s_predkosc(struct s_Predkosc *Wsk_St_predkosc);
void del_st_cisnienie(struct Cisnienie *Wsk_St_cisnienie);
#endif