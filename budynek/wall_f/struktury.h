#include <stdio.h>
#include <stdlib.h>
#ifndef struktury_h
#define struktury_h
struct IOdata
{
    double predkosc;
    double grad_cisnienia;
    double grad_ut;
    double cell_dist;

};

struct Wymiary
{
	unsigned int Ilosc_Elem_X; // ilosc elementow na kierunek X
	unsigned int Ilosc_Elem_Y;
	unsigned int Ilosc_Elem_Z;
	unsigned int NI; // ilość nodow na kierunek N-S zawierających się w domenie- nieuwzględnia lini I=0 i I="N+1" poza domeną
	unsigned int NJ; // Na kierunek E-W
	unsigned int NK; // Na kierunek T-B
	double Wymiar_Elem_X; // dlugosc boku elementu
	double Wymiar_Elem_Y;
	double Wymiar_Elem_Z;
	double Ax; // pow. elementu 
	double Ay;
	double Az;
	unsigned int ilosc_wekt_u;
	unsigned int ilosc_u_IJK[3]; // ilosc wektorow u Na kierunku N=S
	unsigned int ilosc_wekt_v;
	unsigned int ilosc_v_IJK[3];
	unsigned int ilosc_wekt_w;
	unsigned int ilosc_w_IJK[3];
	unsigned int ilosc_p_IJK[3];
	unsigned int ilosc_wezlow; // ilosc wezlow cisnienia 
	unsigned int ilosc_d_u; // I*K*2 bo inlet i outlet
	unsigned int ilosc_d_v;
	unsigned int ilosc_d_w;
	double Min;
	unsigned int ilosc_ap_u;
	unsigned int ilosc_ap_v;
	unsigned int ilosc_ap_w;
	unsigned int ilosc_ap_p;
	double p_relax;
	double u_relax;
	double residual;
	unsigned int iter_max;
	short int wymiar;

	double *d_u;
	double *d_v;
	double *d_w;
};
struct n_Predkosc // struktura niosąca wartość prędkości w n-tej iteracji
{
	double *as;
	double *an;
	double *ae;
	double *aw;
	double *at;
	double *ab;
	double *ap;
	double *predkosc;
	double *d;
	double **A_wspol;
	const double *d_ini;
	const double *predkosc_cal;
	unsigned int *ilosc_wekt_IJK;
	const unsigned int *ilosc_wekt; // rozmiar tablic ap i predkosc; 
	const unsigned int *ilosc_wekt_d;
};
struct s_Predkosc // struktura niosąca wartość prędkości z poprzedniej iteracji
{
	double *u;
	double *v;
	double *w;
    double *s_u;
    double *s_v;
    double *s_w;
	double *ap_u;
	double *ap_v;
	double *ap_w;
	 const unsigned int *ilosc_wekt_u;
	const unsigned int *ilosc_wekt_v;
	const unsigned int *ilosc_wekt_w;
};
struct Cisnienie
{
	double *p_correction;
	double *p;
	double *as;
	double *an;
	double *ae;
	double *aw;
	double *at;
	double *ab;
	double *ap;
	double *d;

	unsigned int ilosc_wezlow;


};
struct Plyn
{
	double mi;
	double rho;
};
#endif
