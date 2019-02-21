#include "f_struktury.h"
struct n_Predkosc *alokacja_n_p(const double *dini, unsigned int *ilosc_wekt_IJK)
{
	int i;
	struct n_Predkosc *Wsk_St_predkosc = (struct n_Predkosc*)malloc(sizeof(struct n_Predkosc));
	if (Wsk_St_predkosc == NULL) { ERROR_ALLOC; return NULL; }

    Wsk_St_predkosc->an = (double*)calloc(ilosc_wekt_IJK[0], sizeof(double));
	if (Wsk_St_predkosc->an == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->as = (double*)calloc(ilosc_wekt_IJK[0], sizeof(double));
	if (Wsk_St_predkosc->as == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->aw = (double*)calloc(ilosc_wekt_IJK[0], sizeof(double));
	if (Wsk_St_predkosc->aw == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->ae = (double*)calloc(ilosc_wekt_IJK[0], sizeof(double));
	if (Wsk_St_predkosc->ae == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->at = (double*)calloc(ilosc_wekt_IJK[0] , sizeof(double));
	if (Wsk_St_predkosc->at == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->ab = (double*)calloc(ilosc_wekt_IJK[0], sizeof(double));
	if (Wsk_St_predkosc->ab == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->d = (double*)calloc(ilosc_wekt_IJK[0] , sizeof(double));
	if (Wsk_St_predkosc->d == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->ap = (double*)calloc(ilosc_wekt_IJK[0]  , sizeof(double));
	if (Wsk_St_predkosc->ap == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->predkosc = (double*)calloc(ilosc_wekt_IJK[0] , sizeof(double));
    if (Wsk_St_predkosc->predkosc == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

	
	//Wsk_St_predkosc->A_wspol = (double**)calloc(ilosc_wekt_IJK[0],sizeof(double*));
	//for ( i = 0; i < ilosc_wekt_IJK[0]; i++)
		//Wsk_St_predkosc->A_wspol[i] = (double*)calloc(ilosc_wekt_IJK[0], sizeof(double));

	
	Wsk_St_predkosc->ilosc_wekt_IJK = ilosc_wekt_IJK;

	Wsk_St_predkosc->d_ini = dini;
	return Wsk_St_predkosc;
}

struct s_Predkosc *alokacja_s_p( struct Wymiary *Wsk_St_wymiary)
{
	struct s_Predkosc *Wsk_St_predkosc = (struct s_Predkosc*)malloc(sizeof(struct s_Predkosc));
	if (Wsk_St_predkosc == NULL) { ERROR_ALLOC; return NULL; }

	Wsk_St_predkosc->u = (double*)calloc(Wsk_St_wymiary->ilosc_wekt_u,  sizeof(double));
	if (Wsk_St_predkosc->u == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }
	
	Wsk_St_predkosc->v = (double*)calloc(Wsk_St_wymiary->ilosc_wekt_v,  sizeof(double));
	
	if (Wsk_St_predkosc->v == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }
	Wsk_St_predkosc->w = (double*)calloc(Wsk_St_wymiary->ilosc_wekt_w, sizeof(double));
	if (Wsk_St_predkosc->w == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->s_u = (double*)calloc(Wsk_St_wymiary->ilosc_wekt_u,  sizeof(double));
    if (Wsk_St_predkosc->s_u == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

    Wsk_St_predkosc->s_v = (double*)calloc(Wsk_St_wymiary->ilosc_wekt_v,  sizeof(double));

    if (Wsk_St_predkosc->s_v == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }
    Wsk_St_predkosc->s_w = (double*)calloc(Wsk_St_wymiary->ilosc_wekt_w, sizeof(double));
    if (Wsk_St_predkosc->s_w == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }



	Wsk_St_predkosc->ap_u = (double*)calloc(Wsk_St_wymiary->ilosc_ap_u, sizeof(double));
	if (Wsk_St_predkosc->u == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

	Wsk_St_predkosc->ap_v = (double*)calloc(Wsk_St_wymiary->ilosc_ap_v, sizeof(double));
	if (Wsk_St_predkosc->v == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

	Wsk_St_predkosc->ap_w = (double*)calloc(Wsk_St_wymiary->ilosc_ap_w, sizeof(double));
	if (Wsk_St_predkosc->w == NULL){ free(Wsk_St_predkosc); ERROR_ALLOC;  return NULL; }

	Wsk_St_predkosc->ilosc_wekt_u = &Wsk_St_wymiary->ilosc_wekt_u;
	Wsk_St_predkosc->ilosc_wekt_v = &Wsk_St_wymiary->ilosc_wekt_v;
	Wsk_St_predkosc->ilosc_wekt_w = &Wsk_St_wymiary->ilosc_wekt_w;
	return Wsk_St_predkosc;
}
struct Cisnienie *alokacja_c(const struct Wymiary *Wsk_St_wymiary)
{
	struct Cisnienie *Wsk_St_cisnienie = (struct Cisnienie*)malloc(sizeof(struct Cisnienie));
	if (Wsk_St_cisnienie == NULL) { ERROR_ALLOC; return NULL; }

	Wsk_St_cisnienie->p_correction = (double*)calloc(Wsk_St_wymiary->ilosc_wezlow,sizeof (double));
	if (Wsk_St_cisnienie->p_correction == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->p = (double*)calloc(Wsk_St_wymiary->ilosc_wezlow,sizeof(double));
	if (Wsk_St_cisnienie->p == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }


	Wsk_St_cisnienie->an = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0] * Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->an == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->as = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0]* Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->as == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->aw = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0] * Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->aw == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->ae = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0] * Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->ae == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->at = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0] * Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->at == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->ab = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0] * Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->ab == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->d = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0] * Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->d == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }

	Wsk_St_cisnienie->ap = (double*)calloc(Wsk_St_wymiary->ilosc_p_IJK[0] * Wsk_St_wymiary->ilosc_p_IJK[1] , sizeof(double));
	if (Wsk_St_cisnienie->ap == NULL){ free(Wsk_St_cisnienie); ERROR_ALLOC;  return NULL; }


	Wsk_St_cisnienie->ilosc_wezlow = Wsk_St_wymiary->ilosc_wezlow;
	return Wsk_St_cisnienie;
}
struct Plyn *alokacja_plyn(const struct Wymiary *Wsk_St_wymiary)
{
	struct Plyn *Wsk_St_plyn = (struct Plyn*)malloc(sizeof(struct Plyn));
	if (Wsk_St_plyn == NULL) { ERROR_ALLOC; return NULL; }
    Wsk_St_plyn->mi = 17e-6;
    Wsk_St_plyn->rho = 1.2;
	return Wsk_St_plyn;

}
void del_st_n_Predkosc(struct n_Predkosc *Wsk_St_predkosc)
{
	int i;
	if (Wsk_St_predkosc != NULL)
	{
		free(Wsk_St_predkosc->an);
		free(Wsk_St_predkosc->as);
		free(Wsk_St_predkosc->aw);
		free(Wsk_St_predkosc->ae);
		free(Wsk_St_predkosc->at);
		free(Wsk_St_predkosc->ab);
		free(Wsk_St_predkosc->d);
		//for (i = 0; i < Wsk_St_predkosc->ilosc_wekt_IJK[0]; i++)
		//	free(Wsk_St_predkosc->A_wspol[i]);
//		free(Wsk_St_predkosc->A_wspol);
		//free(Wsk_St_predkosc->predkosc);
		free(Wsk_St_predkosc);

	}
}
void del_st_s_predkosc(struct s_Predkosc *Wsk_St_predkosc)
{
	if (Wsk_St_predkosc != NULL)
	{
		free(Wsk_St_predkosc->ap_u);
		free(Wsk_St_predkosc->ap_v);
		free(Wsk_St_predkosc->ap_w);
		free(Wsk_St_predkosc->u);
		free(Wsk_St_predkosc->v);
		free(Wsk_St_predkosc->w);
		free(Wsk_St_predkosc);
	}
}
void del_st_cisnienie(struct Cisnienie *Wsk_St_cisnienie)
{
	if (Wsk_St_cisnienie != NULL)
	{
		free(Wsk_St_cisnienie->p_correction);
		free(Wsk_St_cisnienie);
	}
}
