#define _CRT_SECURE_NO_DEPRECATE 
#include "funkcje.h"
void czytaj_rozmiar(struct Wymiary *Wsk_St_wymiary, struct IOdata *Wsk_St_IOdata)
{
    double Wielk_domeny_X=Wsk_St_IOdata->cell_dist;
	double Wielk_domeny_Y;
	double Wielk_domeny_Z;


	/*printf("Podaj rozmiar domeny XxYxZ[m]:");
	scanf("%fx%fx%f", &Wielk_domeny_X, &Wielk_domeny_Y, &Wielk_domeny_Z);
	printf("%f %f %f\n", Wielk_domeny_X, Wielk_domeny_Y, Wielk_domeny_Z);

	printf("Podaj ilosc elementow na krawedz XxYxZ");
	scanf("%dx%dx%d", &Wsk_St_wymiary->Ilosc_Elem_X, &Wsk_St_wymiary->Ilosc_Elem_Y, &Wsk_St_wymiary->Ilosc_Elem_Z);
	printf("\n");*/
	/*printf("podaj parametry relaxacji cisnieniexpredkosc");
	scanf("%fx%f", &Wsk_St_wymiary->p_relax, &Wsk_St_wymiary->u_relax);
	printf("\n"); */ 
    /*Wielk_domeny_X = 0.005;*/ Wielk_domeny_Y = 0.1; Wielk_domeny_Z = 0.1;
    
	Wsk_St_wymiary->p_relax = 0.4;
    Wsk_St_wymiary->u_relax = 0.6;
    Wsk_St_wymiary->residual = 0.003;
    Wsk_St_wymiary->iter_max = 1000;
    
    Wsk_St_wymiary->wymiar = 1;
	
    Wsk_St_wymiary->Ilosc_Elem_X = 12; Wsk_St_wymiary->Ilosc_Elem_Y =0; 	Wsk_St_wymiary->Ilosc_Elem_Z = 0;
    Wsk_St_wymiary->Wymiar_Elem_X =(double)(Wielk_domeny_X / Wsk_St_wymiary->Ilosc_Elem_X);
	Wsk_St_wymiary->Wymiar_Elem_Y = 1;
	Wsk_St_wymiary->Wymiar_Elem_Z = 1;
	

	Wsk_St_wymiary->NI = 3 + (Wsk_St_wymiary->Ilosc_Elem_X - 1) * 2;
	Wsk_St_wymiary->NJ = 3;
	Wsk_St_wymiary->NK = 3;

	 
	Wsk_St_wymiary->Ax = Wsk_St_wymiary->Wymiar_Elem_X;
	Wsk_St_wymiary->Ay = 0;
	Wsk_St_wymiary->Az = 0;
	Wsk_St_wymiary->Min = 0;

	if (Wsk_St_wymiary->wymiar >1)
	{
		Wsk_St_wymiary->Ax = Wsk_St_wymiary->Wymiar_Elem_Y;
		Wsk_St_wymiary->Ay = Wsk_St_wymiary->Wymiar_Elem_X;
		Wsk_St_wymiary->NJ = 3 + (Wsk_St_wymiary->Ilosc_Elem_Y - 1) * 2;
		Wsk_St_wymiary->ilosc_d_v = ((Wsk_St_wymiary->NJ - 3) / 2)*((Wsk_St_wymiary->NK - 1) / 2) * 2;
        Wsk_St_wymiary->Wymiar_Elem_Y = 1;//(double)(Wielk_domeny_Y / Wsk_St_wymiary->Ilosc_Elem_Y);
	}
	
	if (Wsk_St_wymiary->wymiar == 3)
	{
		Wsk_St_wymiary->Ilosc_Elem_Z = 150;
		Wsk_St_wymiary->Wymiar_Elem_Z = (double)(Wielk_domeny_Z / Wsk_St_wymiary->Ilosc_Elem_Z);
		Wsk_St_wymiary->ilosc_wekt_w = ((Wsk_St_wymiary->NI + 1) / 2)*((Wsk_St_wymiary->NJ + 3) / 2)*((Wsk_St_wymiary->NK + 1) / 2);
		Wsk_St_wymiary->ilosc_wekt_w = ((Wsk_St_wymiary->NI + 1) / 2)*((Wsk_St_wymiary->NJ + 3) / 2)*((Wsk_St_wymiary->NK + 1) / 2);
		Wsk_St_wymiary->ilosc_w_IJK[0] = (Wsk_St_wymiary->NI + 1) / 2;
		Wsk_St_wymiary->ilosc_w_IJK[1] = (Wsk_St_wymiary->NJ + 3) / 2;
		Wsk_St_wymiary->ilosc_w_IJK[2] = (Wsk_St_wymiary->NK + 1) / 2;
		Wsk_St_wymiary->ilosc_ap_w = (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ - 1)*(Wsk_St_wymiary->NK - 3) / 8;
		Wsk_St_wymiary->ilosc_d_w = ((Wsk_St_wymiary->NJ - 1) / 2)*((Wsk_St_wymiary->NK - 3) / 2) * 2;
		Wsk_St_wymiary->Ax = Wsk_St_wymiary->Wymiar_Elem_Y*Wsk_St_wymiary->Wymiar_Elem_Z;
		Wsk_St_wymiary->Ay = Wsk_St_wymiary->Wymiar_Elem_X*Wsk_St_wymiary->Wymiar_Elem_Z;
		Wsk_St_wymiary->Az = Wsk_St_wymiary->Wymiar_Elem_X*Wsk_St_wymiary->Wymiar_Elem_Y;
		Wsk_St_wymiary->NK = 3 + (Wsk_St_wymiary->Ilosc_Elem_Z - 1) * 2;
	}
    // -4*X wynika z braku wektorów w kacie domeny
    Wsk_St_wymiary->ilosc_wekt_u = ((Wsk_St_wymiary->NI - 1) / 2)*((Wsk_St_wymiary->NJ + 3) / 2)*((Wsk_St_wymiary->NK + 3) / 2);//-4*((Wsk_St_wymiary->NI - 1) / 2);
    Wsk_St_wymiary->ilosc_wekt_v = ((Wsk_St_wymiary->NI + 1) / 2)*((Wsk_St_wymiary->NJ + 1) / 2)*((Wsk_St_wymiary->NK + 3) / 2);//-4*((Wsk_St_wymiary->NI + 1) / 2);
    Wsk_St_wymiary->ilosc_wekt_w = ((Wsk_St_wymiary->NI + 1) / 2)*((Wsk_St_wymiary->NJ + 3) / 2)*((Wsk_St_wymiary->NK + 1) / 2);//-4*((Wsk_St_wymiary->NI + 1) / 2);

	Wsk_St_wymiary->ilosc_u_IJK[0] = (Wsk_St_wymiary->NI - 1) / 2;
	Wsk_St_wymiary->ilosc_u_IJK[1] = (Wsk_St_wymiary->NJ + 3) / 2;
	Wsk_St_wymiary->ilosc_u_IJK[2] = (Wsk_St_wymiary->NK + 3) / 2;

	Wsk_St_wymiary->ilosc_v_IJK[0] = (Wsk_St_wymiary->NI + 1) / 2;
	Wsk_St_wymiary->ilosc_v_IJK[1] = (Wsk_St_wymiary->NJ + 1) / 2;
	Wsk_St_wymiary->ilosc_v_IJK[2] = (Wsk_St_wymiary->NK + 3) / 2;

	Wsk_St_wymiary->ilosc_w_IJK[0] = (Wsk_St_wymiary->NI + 1) / 2;
	Wsk_St_wymiary->ilosc_w_IJK[1] = (Wsk_St_wymiary->NJ + 3) / 2;
	Wsk_St_wymiary->ilosc_w_IJK[2] = (Wsk_St_wymiary->NK + 1) / 2;

	Wsk_St_wymiary->ilosc_ap_u = (Wsk_St_wymiary->NI - 1)*(Wsk_St_wymiary->NJ - 1)*(Wsk_St_wymiary->NK - 1) / 8;
	Wsk_St_wymiary->ilosc_ap_v = (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ - 3)*(Wsk_St_wymiary->NK - 1) / 8;
	Wsk_St_wymiary->ilosc_ap_w = (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ - 1)*(Wsk_St_wymiary->NK - 3) / 8;
	Wsk_St_wymiary->ilosc_ap_p = (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ - 1)*(Wsk_St_wymiary->NK - 1) / 8;

	Wsk_St_wymiary->ilosc_d_u = ((Wsk_St_wymiary->NJ - 1) / 2)*((Wsk_St_wymiary->NK - 1) / 2) * 2; // Inne wartosci niz dla wektorow uvw bo 
	Wsk_St_wymiary->ilosc_d_v = ((Wsk_St_wymiary->NJ - 3) / 2)*((Wsk_St_wymiary->NK - 1) / 2) * 2; // tam sa wziete pod uwage wektory w scianach-zerowe
	Wsk_St_wymiary->ilosc_d_w = ((Wsk_St_wymiary->NJ - 1) / 2)*((Wsk_St_wymiary->NK - 3) / 2) * 2;

    Wsk_St_wymiary->ilosc_wezlow = ((Wsk_St_wymiary->NI - 3) / 2)*((Wsk_St_wymiary->NJ + 3) / 2)*((Wsk_St_wymiary->NK + 3) / 2); // cisnieniowych

    Wsk_St_wymiary->ilosc_p_IJK[0] = (Wsk_St_wymiary->NI -3) / 2;
	Wsk_St_wymiary->ilosc_p_IJK[1] = (Wsk_St_wymiary->NJ + 3) / 2;
	Wsk_St_wymiary->ilosc_p_IJK[2] = (Wsk_St_wymiary->NK + 3) / 2;
	/*printf("Podaj dlugosci krawedzi elementu w XxYxZ[mm]");
	scanf("%lfx%lfx%lf", &Wsk_St_wymiary->Wymiar_Elem_X, &Wsk_St_wymiary->Wymiar_Elem_Y, &Wsk_St_wymiary->Wymiar_Elem_Z);
	printf("\n");

	Wsk_St_wymiary->Ilosc_Elem_X = (int)(Wielk_domena_X / Wsk_St_wymiary->Wymiar_Elem_X);
	Wsk_St_wymiary->Ilosc_Elem_Y = (int)(Wielk_domena_Y / Wsk_St_wymiary->Wymiar_Elem_Y);
	Wsk_St_wymiary->Ilosc_Elem_Z = (int)(Wielk_domena_Z / Wsk_St_wymiary->Wymiar_Elem_Z);

	if ((int)Wielk_domena_X % (int)Wsk_St_wymiary->Wymiar_Elem_X != 0)
	{
	}*/

}
void inicjalizacja(struct Wymiary *Wsk_St_Wymiary, struct Plyn *Wsk_St_plyn, struct s_Predkosc *Wsk_St_predkosc,struct Cisnienie *Wsk_St_cisnienie, struct IOdata *Wsk_St_IOdata)
{
	unsigned int i,k, j,l,m;
	float Mout=0;
	//inicjalizacja inlet
	Wsk_St_Wymiary->d_u = (double*)calloc(Wsk_St_Wymiary->ilosc_d_u, sizeof(double));
	if (Wsk_St_Wymiary->d_u == NULL){ free(Wsk_St_Wymiary); ERROR_ALLOC; }

	Wsk_St_Wymiary->d_v = (double*)calloc(Wsk_St_Wymiary->ilosc_d_v, sizeof(double));
	if (Wsk_St_Wymiary->d_v == NULL){ free(Wsk_St_Wymiary); ERROR_ALLOC; }

	Wsk_St_Wymiary->d_w = (double*)calloc(Wsk_St_Wymiary->ilosc_d_w, sizeof(double));
	if (Wsk_St_Wymiary->d_w == NULL){ free(Wsk_St_Wymiary); ERROR_ALLOC; }

    for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
    for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
    {
        Wsk_St_Wymiary->d_u[translation(Wsk_St_Wymiary, 2, j, k, 4)] = 0.0;
        Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, 2, j, k, 0)] = Wsk_St_Wymiary->d_u[translation(Wsk_St_Wymiary, 2, j, k, 4)];
//		Wsk_St_Wymiary->Min += Wsk_St_Wymiary->d_u[translation(Wsk_St_Wymiary, 2, j, k, 4)];
//		//printf("I=\t\tPREDKOSC UN=%f\t\tPREDKOSC US=\n", Wsk_St_Wymiary->d_u[translation(Wsk_St_Wymiary, 2, j, k, 4)]);
        for (i = 4; i < Wsk_St_Wymiary->NI - 2; i += 2){
            Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)] =0.1; //(double) (rand() % (int) Wsk_St_IOdata->predkosc);
        }
          Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, Wsk_St_Wymiary->NI-1, j, k, 0)] =Wsk_St_IOdata->predkosc;
     //   for (i = 1; i < Wsk_St_Wymiary->NI + 1; i += 2){
       //             Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i, j, k, 3)] = 1.0;// (rand() % 10);

               // }
    }







	Wsk_St_Wymiary->Min *= Wsk_St_plyn->rho*Wsk_St_Wymiary->Ax;
}
int translation(struct Wymiary *Wsk_St_wymiary, unsigned int I, unsigned int J, unsigned int K, int typ)
{
	// Typy 1-u, 2-v, 3-w, 4- p
	switch (typ)
	{
	case 0:
		return (Wsk_St_wymiary->NI - 1)*J / 4 + (Wsk_St_wymiary->NI - 1)*(Wsk_St_wymiary->NJ + 3) / 8 * K + (I - 2) / 2;
	case 1:
		return (Wsk_St_wymiary->NI + 1)*(J - 1) / 4 + (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ + 1) / 8 * K + (I - 1) / 2;
	case 2:
		return (Wsk_St_wymiary->NI + 1)*J / 4 + (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ + 3) / 8 * (K - 1) + (I - 1) / 2;
	case 3:
		return (Wsk_St_wymiary->NI + 1)*J / 4 + (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ + 3) / 8 * K + (I - 1) / 2;

	case 4:
		return ((Wsk_St_wymiary->NJ - 1) / 2 * (K - 2) / 2 + (J - 2) / 2); // translacja d_u
	case 5:
		return ((Wsk_St_wymiary->NJ - 3) / 2 * (K - 2) / 2 + (J - 3) / 2); // translacja d_v
	case 6:
		return ((Wsk_St_wymiary->NJ - 1) / 2 * (K - 3) / 2 + (J - 2) / 2); // translacja d_w
	case 7:
		return (Wsk_St_wymiary->NI - 1)*(J - 2) / 4 + (Wsk_St_wymiary->NI - 1)*(Wsk_St_wymiary->NJ - 1) / 8 * (K - 2) + (I - 2) / 2;//ap_u
	case 8:
		return (Wsk_St_wymiary->NI + 1)*(J - 3) / 4 + (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ - 3) / 8 * (K - 2) + (I - 1) / 2;//ap_v
	case 9:
		return (Wsk_St_wymiary->NI + 1)*(J - 2) / 4 + (Wsk_St_wymiary->NI + 1)*(Wsk_St_wymiary->NJ - 1) / 8 * (K - 3) + (I - 1) / 2;//ap_w
	default:
		printf("Podano zły typ!!\n");
		break;
	}
	return 0;
}
