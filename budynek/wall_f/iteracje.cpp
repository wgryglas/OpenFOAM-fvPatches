#include "iteracje.h"


void iteracja(struct Wymiary *Wsk_St_Wymiary, struct s_Predkosc *Wsk_St_s_predkosc, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct IOdata *Wsk_St_IOdata)
{
	//
	//FILE *f;
	//f = fopen("plik3.txt", "w");
	struct n_Predkosc *Wsk_St_n_predkosc;
	double *u, *v, *w, *p;
	double loc_residual[] = { 0, 0, 0, 0 };;
	double loc_residual_k[] = { 1, 1, 1, 1 };;
	double norm_residual[] = { 0, 0, 0, 0 };;
    double blad;
	short int stop[] = { 0, 0, 0, 0 };
	u = (double*)calloc(Wsk_St_Wymiary->ilosc_ap_u, sizeof(double));
	v = (double*)calloc(Wsk_St_Wymiary->ilosc_ap_v, sizeof(double));
	w = (double*)calloc(Wsk_St_Wymiary->ilosc_ap_w, sizeof(double));
	p = (double*)calloc(Wsk_St_Wymiary->ilosc_ap_p, sizeof(double));
    unsigned int petla = 0, i, j, k, iter;

	if (Wsk_St_Wymiary->wymiar < 2)
	stop[1] = 1;
	if (Wsk_St_Wymiary->wymiar < 3)
	stop[2] = 1;

while ( petla != Wsk_St_Wymiary->iter_max)
	{

		//obliczam predkosc
        for(i=0;i<Wsk_St_Wymiary->ilosc_wekt_u;i++)
            Wsk_St_s_predkosc->s_u[i]=Wsk_St_s_predkosc->u[i];

        for(i=0;i<Wsk_St_Wymiary->ilosc_wekt_v;i++)
            Wsk_St_s_predkosc->s_v[i]=Wsk_St_s_predkosc->v[i];

        for(i=0;i<Wsk_St_Wymiary->ilosc_wekt_w;i++)
            Wsk_St_s_predkosc->s_w[i]=Wsk_St_s_predkosc->w[i];

        iteracja_predkosci( Wsk_St_Wymiary, Wsk_St_cisnienie, Wsk_St_plyn, Wsk_St_s_predkosc,Wsk_St_IOdata);


//        obliczam cisnienie
//        iteracja_cisnienia( Wsk_St_Wymiary, Wsk_St_cisnienie, Wsk_St_plyn, Wsk_St_s_predkosc);

		//>>>>>>>>>>>>>>>>>>>U<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		iter = 0;
           blad = 0;
		for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
		for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
		for (i = 4; i < Wsk_St_Wymiary->NI; i += 2)
		{
			
            u[iter] = Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)];// +
//                      Wsk_St_Wymiary->Ax / Wsk_St_s_predkosc->ap_u[translation(Wsk_St_Wymiary, i, j, k, 7)] *
                 //    (Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i - 1, j, k, 3)] -
                   //   Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i + 1, j, k, 3)]);


            Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)] = u[iter];//Wsk_St_Wymiary->u_relax*u[iter] +
                  // (1 - Wsk_St_Wymiary->u_relax)*Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)];

            if(blad<fabs(Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)] - Wsk_St_s_predkosc->s_u[translation(Wsk_St_Wymiary, i, j, k, 0)] ))
            blad=fabs(Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)] - Wsk_St_s_predkosc->s_u[translation(Wsk_St_Wymiary, i, j, k, 0)]);
//            printf("\t\t\tBLAD=%.8f\n",blad);

            iter++;
		}
        petla++;
//        printf("PEtla=%d\n",petla);
        if(blad<=0.0001 || (petla == Wsk_St_Wymiary->iter_max))
            break;

        /* OBLICZANIE VW
		//>>>>>>>>>>>>>>>>>>>V<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        iter = 0;
		for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
		for (j = 3; j < Wsk_St_Wymiary->NJ - 1; j += 2)
		for (i = 3; i < Wsk_St_Wymiary->NI + 1; i += 2)
		{			
			v[iter] = Wsk_St_s_predkosc->v[translation(Wsk_St_Wymiary, i, j, k, 1)] + Wsk_St_Wymiary->Ay / Wsk_St_s_predkosc->ap_v[translation(Wsk_St_Wymiary, i, j, k, 8)] * (Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j - 1, k, 3)] - Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j + 1, k, 3)]);
            Wsk_St_s_predkosc->v[translation(Wsk_St_Wymiary, i, j, k, 1)] =v[iter];// + (1 - Wsk_St_Wymiary->u_relax)*Wsk_St_s_predkosc->v[translation(Wsk_St_Wymiary, i, j, k, 1)];
			iter++;
			//printf("v %d\n", u[iter - 1]);
        }
		//>>>>>>>>>>>>>>>>>>>W<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        iter = 0;
		for (k = 3; k < Wsk_St_Wymiary->NK - 1; k += 2)
		for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
		for (i = 3; i < Wsk_St_Wymiary->NI + 1; i += 2)
		{
			w[iter] = Wsk_St_s_predkosc->w[translation(Wsk_St_Wymiary, i, j, k, 2)] + Wsk_St_Wymiary->Az / Wsk_St_s_predkosc->ap_w[translation(Wsk_St_Wymiary, i, j, k, 9)] * (Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j, k - 1, 3)] - Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j, k + 1, 3)]);
			Wsk_St_s_predkosc->w[translation(Wsk_St_Wymiary, i, j, k, 2)] = Wsk_St_Wymiary->u_relax*w[iter] + (1 - Wsk_St_Wymiary->u_relax)*Wsk_St_s_predkosc->w[translation(Wsk_St_Wymiary, i, j, k, 2)];
			iter++;
			//printf("w %d\n", u[iter - 1]);
        }*/
		//>>>>>>>>>>>>>>>>>>>P<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		iter = 0;

		for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
		for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
		for (i = 1; i < Wsk_St_Wymiary->NI + 1; i += 2)
		{
            p[iter] = Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i, j, k, 3)] + Wsk_St_Wymiary->p_relax*Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j, k, 3)];


			Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i, j, k, 3)] = p[iter];
			iter++;
		}

/*//<<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>OBLICZANIE RESIDUA<<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
       //>>>>>>>>>>>>>>>>>>>>>U<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (stop[0] == 0)
		{
			Wsk_St_n_predkosc = alokacja_n_p(Wsk_St_Wymiary->d_u, Wsk_St_Wymiary->ilosc_u_IJK);
			
			for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
			for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
			{
				iter = 0;
				for (i = 4; i < Wsk_St_Wymiary->NI - 2; i += 2)
				{
					obliczanie_predkosci( 0, iter, i, j, k, Wsk_St_Wymiary, Wsk_St_n_predkosc, Wsk_St_s_predkosc, Wsk_St_cisnienie, Wsk_St_plyn);

                    loc_residual[0] += fabs(Wsk_St_n_predkosc->d[iter] + Wsk_St_n_predkosc->an[iter] *
                                            Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i + 2, j, k, 0)] +
                                            Wsk_St_n_predkosc->as[iter] * Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i - 2, j, k, 0)] -
                                            Wsk_St_n_predkosc->ap[iter] * Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)]);
					iter++;
				}

			}
			if (petla == 3 && loc_residual[0] != 0)
				loc_residual_k[0] = loc_residual[0];
			norm_residual[0] = loc_residual[0] / loc_residual_k[0];
			if (norm_residual[0] <= Wsk_St_Wymiary->residual)
				stop[0] = 1;
			del_st_n_Predkosc(Wsk_St_n_predkosc);
        }

        //>>>>>>>>>>>>>>>>>>>>>V<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (stop[1] == 0)
		{
			Wsk_St_n_predkosc = alokacja_n_p(Wsk_St_Wymiary->d_v, Wsk_St_Wymiary->ilosc_v_IJK);
			
			for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
			for (j = 3; j < Wsk_St_Wymiary->NJ-1; j += 2)
			{
				iter = 0;
				for (i = 3; i < Wsk_St_Wymiary->NI - 1; i += 2)
				{
					obliczanie_predkosci( 1, iter, i, j, k, Wsk_St_Wymiary, Wsk_St_n_predkosc, Wsk_St_s_predkosc, Wsk_St_cisnienie, Wsk_St_plyn);
					loc_residual[1] += fabs(Wsk_St_n_predkosc->d[iter] + Wsk_St_n_predkosc->an[iter] * Wsk_St_s_predkosc->v[translation(Wsk_St_Wymiary, i + 2, j, k, 1)] + Wsk_St_n_predkosc->as[iter] * Wsk_St_s_predkosc->v[translation(Wsk_St_Wymiary, i - 2, j, k, 1)] - Wsk_St_n_predkosc->ap[iter] * Wsk_St_s_predkosc->v[translation(Wsk_St_Wymiary, i, j, k, 1)]);
					iter++;
				}

			}
			if (petla == 3 && loc_residual[1] != 0)
				loc_residual_k[1] = loc_residual[1];
			norm_residual[1] = loc_residual[1] / loc_residual_k[1];
			if (norm_residual[1] <= Wsk_St_Wymiary->residual)
				stop[1] = 1;
			del_st_n_Predkosc(Wsk_St_n_predkosc);
        }
//>>>>>>>>>>>>>>>>>>>>>W<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<		
        if (stop[2] == 0 )
        {
			Wsk_St_n_predkosc = alokacja_n_p(Wsk_St_Wymiary->d_w, Wsk_St_Wymiary->ilosc_w_IJK);
			
			for (k = 3; k < Wsk_St_Wymiary->NK-1; k += 2)
			for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
			{
				iter = 0;
				for (i = 3; i < Wsk_St_Wymiary->NI - 1; i += 2)
				{
					obliczanie_predkosci( 2, iter, i, j, k, Wsk_St_Wymiary, Wsk_St_n_predkosc, Wsk_St_s_predkosc, Wsk_St_cisnienie, Wsk_St_plyn);
					loc_residual[2] += fabs(Wsk_St_n_predkosc->d[iter] + Wsk_St_n_predkosc->an[iter] * Wsk_St_s_predkosc->w[translation(Wsk_St_Wymiary, i + 2, j, k, 2)] + Wsk_St_n_predkosc->as[iter] * Wsk_St_s_predkosc->w[translation(Wsk_St_Wymiary, i - 2, j, k, 2)] - Wsk_St_n_predkosc->ap[iter] * Wsk_St_s_predkosc->w[translation(Wsk_St_Wymiary, i, j, k, 2)]);
					iter++;
				}

			}
			if (petla == 3 && loc_residual[2] != 0)
				loc_residual_k[2] = loc_residual[2];
			norm_residual[2] = loc_residual[2] / loc_residual_k[2];
			if (norm_residual[2] <= Wsk_St_Wymiary->residual)
				stop[2] = 1;
			del_st_n_Predkosc(Wsk_St_n_predkosc);
        }
//>>>>>>>>>>>>>>>>>>>>>P<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<		
        if (stop[3] == 0)
		{				
			for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
			for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
			{
				iter = 0;
				for (i = 3; i < Wsk_St_Wymiary->NI - 1; i += 2)
				{
					obliczanie_cisnienia( iter, i, j, k, Wsk_St_Wymiary, Wsk_St_s_predkosc, Wsk_St_cisnienie, Wsk_St_plyn);
					
					loc_residual[3] += fabs(Wsk_St_cisnienie->d[iter] - Wsk_St_cisnienie->an[iter] * Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i+2, j , k, 3)] - Wsk_St_cisnienie->as[iter] * Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i-2, j , k, 3)] - Wsk_St_cisnienie->ap[iter] * Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j , k, 3)]);
					iter++;
				}
			}
			if (petla == 3 && loc_residual[3]!=0)
				loc_residual_k[3] = loc_residual[3];
			norm_residual[3] = loc_residual[3] / loc_residual_k[3];
			if (norm_residual[3] <= Wsk_St_Wymiary->residual)
				stop[3] = 1;
            }
        //printf("\n\n\tpetla=%d \t u=%f \t v=%f \t w=%f \t p=%f \n", petla, norm_residual[0], norm_residual[1], norm_residual[2], norm_residual[3]);
//<<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
        if ((stop[0] == 1 && stop[1] == 1 && stop[2] == 1) || (petla == Wsk_St_Wymiary->iter_max))
        break;
       petla++;*/
	}


	free(u); free(v); free(w); free(p);




}
void iteracja_predkosci(struct Wymiary *Wsk_St_Wymiary, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct s_Predkosc *Wsk_St_s_predkosc, struct IOdata *Wsk_St_IOdata)
{
	struct n_Predkosc *Wsk_St_n_predkosc;
    unsigned int i, j, k,l,m, iter, licznik;
    short unsigned int  mod_i, mod_j, mod_k, indeks;
	double *cp;
	double *dp;
	double *A_wspol;

    for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
    for (indeks = 0; indeks <Wsk_St_Wymiary->wymiar; indeks++)
	{

        if (indeks == 0){
			mod_i = 0; mod_j = 0;     mod_k = k;
			Wsk_St_n_predkosc = alokacja_n_p(Wsk_St_Wymiary->d_u, Wsk_St_Wymiary->ilosc_u_IJK);
		}
        else if (indeks == 1){
			mod_i = 1; mod_j = 1; mod_k = k;
			Wsk_St_n_predkosc = alokacja_n_p(Wsk_St_Wymiary->d_v, Wsk_St_Wymiary->ilosc_v_IJK);
		}
        else{
			mod_i = 1; mod_j = 0;     mod_k = k + 1; if (mod_k == Wsk_St_Wymiary->NK) continue;
			Wsk_St_n_predkosc = alokacja_n_p(Wsk_St_Wymiary->d_w, Wsk_St_Wymiary->ilosc_w_IJK);
		}
		iter = 0;


        for (j = 2+mod_j; j < Wsk_St_Wymiary->NJ-mod_j; j += 2){

            for (i = 2 - mod_i; i < Wsk_St_Wymiary->NI + mod_i; i += 2){
                obliczanie_predkosci(indeks, iter, i, j, mod_k, Wsk_St_Wymiary, Wsk_St_n_predkosc, Wsk_St_s_predkosc, Wsk_St_cisnienie, Wsk_St_plyn,Wsk_St_IOdata);
				iter++;
			}

            A_wspol = (double*)calloc(
                        Wsk_St_n_predkosc->ilosc_wekt_IJK[0] *
                        Wsk_St_n_predkosc->ilosc_wekt_IJK[0]
                        ,sizeof(double));

            for (l = 0; l < Wsk_St_n_predkosc->ilosc_wekt_IJK[0]; l++)
            for (m = 0; m < Wsk_St_n_predkosc->ilosc_wekt_IJK[0]; m++)
            {

                if (l == m) A_wspol[l*Wsk_St_n_predkosc->ilosc_wekt_IJK[0] + m] = Wsk_St_n_predkosc->ap[l];
                if (m == l - 1)A_wspol[l*Wsk_St_n_predkosc->ilosc_wekt_IJK[0] + m] = Wsk_St_n_predkosc->as[l];
                if (m == l + 1)A_wspol[l*Wsk_St_n_predkosc->ilosc_wekt_IJK[0] + m] = Wsk_St_n_predkosc->an[l];
            }


            alglib::real_2d_array  algorytm;
            algorytm.setcontent(Wsk_St_n_predkosc->ilosc_wekt_IJK[0],
                Wsk_St_n_predkosc->ilosc_wekt_IJK[0], A_wspol);
            alglib::real_1d_array wynik;
            alglib::real_1d_array D;
            D.setcontent(Wsk_St_n_predkosc->ilosc_wekt_IJK[0] , Wsk_St_n_predkosc->d);
            alglib::ae_int_t E;
            alglib::densesolverreport REP;
            alglib::rmatrixsolve(algorytm, Wsk_St_n_predkosc->ilosc_wekt_IJK[0] , D, E, REP, wynik);

            Wsk_St_n_predkosc->predkosc= wynik.getcontent();

            switch (indeks) {//przepisywanie obliczonej predkosci do tablicy glownej
            case 0:
                licznik = 0;
                for (i = 2; i < Wsk_St_Wymiary->NI; i += 2)
                {
                    Wsk_St_s_predkosc->ap_u[translation(Wsk_St_Wymiary, i, j, mod_k, 7)] = Wsk_St_n_predkosc->ap[licznik];
                    Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, mod_k, 0)] = Wsk_St_n_predkosc->predkosc[licznik];
                    //fprintf(f, "\i=%d\tj=%d\tk=%d\tcos=%f\t\n", iter, licznik, k, Wsk_St_n_predkosc->ap[licznik]);
                    licznik++;
                }
                break;
            case 1:
                licznik = 0;
                for (i = 1; i < Wsk_St_Wymiary->NI + 1; i += 2)
                {
                    Wsk_St_s_predkosc->ap_v[translation(Wsk_St_Wymiary, i, j, mod_k, 8)] = Wsk_St_n_predkosc->ap[licznik];
                    Wsk_St_s_predkosc->v[translation(Wsk_St_Wymiary, i, j, mod_k, 1)] = Wsk_St_n_predkosc->predkosc[licznik];
                    licznik++;
                }
                break;
            case 2:
                licznik = 0;
                for (i = 1; i < Wsk_St_Wymiary->NI + 1; i += 2)
                {
                    Wsk_St_s_predkosc->ap_w[translation(Wsk_St_Wymiary, i, j, mod_k, 9)] = Wsk_St_n_predkosc->ap[licznik];
                    Wsk_St_s_predkosc->w[translation(Wsk_St_Wymiary, i, j, mod_k, 2)] = Wsk_St_n_predkosc->predkosc[licznik];
                    //fprintf(f, "\ti=%d\tj=%d\tk=%d\titer=%d\tdp=%f\tcp=\n", i, mod_j, mod_k, iter, Wsk_St_s_predkosc->ap_w[translation(Wsk_St_Wymiary, i, mod_j, mod_k, 9)]);
                    licznik++;
                }
                break;
            default:
                printf("indeks=%d out of border -> 0-u,1-v,2-w",indeks);
            }

            free(A_wspol);
        }

			del_st_n_Predkosc(Wsk_St_n_predkosc);


		
	}
}

void iteracja_cisnienia(struct Wymiary *Wsk_St_Wymiary, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct s_Predkosc *Wsk_St_s_predkosc)
{
	FILE *f;
//    f = fopen("plik332.txt","w");
    unsigned int i, j, k,l,m, iter, licznik;
	double *cp, *dp, *temp, *A_wspol;

	for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
	for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
	{
//         printf("FLAGA_2.1\n");
		iter = 0;
        for (i = 3; i < Wsk_St_Wymiary->NI - 1; i += 2)
		{
			obliczanie_cisnienia(iter, i, j, k, Wsk_St_Wymiary, Wsk_St_s_predkosc, Wsk_St_cisnienie, Wsk_St_plyn);
			iter++;
             printf("FLAGA_2.2_%d\n",iter);
		}

		// OBLICZANIE przeganianiem
		//>>>>>>>>>>>>>>

		A_wspol = (double*)calloc(Wsk_St_Wymiary->ilosc_p_IJK[0] * Wsk_St_Wymiary->ilosc_p_IJK[0], sizeof(double));

        for (l = 0; l < Wsk_St_Wymiary->ilosc_p_IJK[0];l++)
                            printf("\tcisnienie=%f\n",  Wsk_St_cisnienie->ap[l]);

		for (l = 0; l < Wsk_St_Wymiary->ilosc_p_IJK[0]; l++){
			for (m = 0; m < Wsk_St_Wymiary->ilosc_p_IJK[0]; m++)
			{
				if (l == m) A_wspol[l*Wsk_St_Wymiary->ilosc_p_IJK[0] + m] = Wsk_St_cisnienie->ap[l];
//				fprintf(f,"%lf", Wsk_St_cisnienie->an[l]);
				if (m == l - 1)A_wspol[l*Wsk_St_Wymiary->ilosc_p_IJK[0] + m] = Wsk_St_cisnienie->as[l];
                if (m == l + 1)A_wspol[l*Wsk_St_Wymiary->ilosc_p_IJK[0] + m] = Wsk_St_cisnienie->an[l];
				//fprintf(f,"%lf", A_wspol[l*Wsk_St_Wymiary->ilosc_p_IJK[0] + m]);
			}
			//fprintf(f,"\n");
		}
//        printf("ILOSC=%d",Wsk_St_Wymiary->ilosc_p_IJK[0]);

        for (l = 0; l < Wsk_St_Wymiary->ilosc_p_IJK[0]; l++){
            for (m = 0; m < Wsk_St_Wymiary->ilosc_p_IJK[0]; m++){
                fprintf(f,"\t%.8f\t",A_wspol[l*Wsk_St_Wymiary->ilosc_p_IJK[0]+m] );}
                fprintf(f,"\n");
        }
		fclose(f);
        printf("FLAGA_2.4\n");
		alglib::real_2d_array  algorytm;
		algorytm.setcontent(Wsk_St_Wymiary->ilosc_p_IJK[0], Wsk_St_Wymiary->ilosc_p_IJK[0], A_wspol);
		alglib::real_1d_array wynik;
		alglib::real_1d_array D;
		D.setcontent(Wsk_St_Wymiary->ilosc_p_IJK[0], Wsk_St_cisnienie->d);
		alglib::ae_int_t E;
		alglib::densesolverreport REP;
		alglib::rmatrixsolve(algorytm, Wsk_St_Wymiary->ilosc_p_IJK[0], D, E, REP, wynik);


		temp = wynik.getcontent();
//                for (l = 0; l < Wsk_St_Wymiary->ilosc_p_IJK[0];l++)
//                    printf("\tcisnienie=%f\n", temp[l]);// Wsk_St_cisnienie->d[l]);
		//	for (m = 0; m < iter; m++)
				//printf("\tD=%lf\n", temp[m]);
		//	printf("\n");


//        printf("FLAGA_2.5\n");
		//cp = (double*)calloc(iter, sizeof(double));
		//dp = (double*)calloc(iter, sizeof(double));
		//temp = (double*)calloc(iter, sizeof(double));

		//cp[0] = 0;
		//dp[0] = 0;

		//for (i = 1; i < iter; i++)
		//{
		//	cp[i] = Wsk_St_cisnienie->an[i] / (Wsk_St_cisnienie->ap[i] - Wsk_St_cisnienie->as[i] * cp[i - 1]);
		//	dp[i] = (Wsk_St_cisnienie->d[i] + Wsk_St_cisnienie->as[i] * dp[i - 1]) / (Wsk_St_cisnienie->ap[i] - Wsk_St_cisnienie->as[i] * cp[i - 1]);
		//	//fprintf(f, "\ti=%d\tj=%d\tk=%d\titer=%d\tdp=%f\tcp=%f\n", i, j, k, iter, Wsk_St_cisnienie->d[i], cp[i]);
		//}

		//temp[iter - 1] = dp[iter - 1];

		//for (i = iter - 2; i >= 0; i--)
		//{
		//	temp[i] = dp[i] + cp[i] * temp[i + 1]; // zapisywanie wyników
		//	//fprintf(f, "\ti=%d\tj=%d\tk=%d\tin=%d\t%f\n", i, j, k, iter, temp[i]);
		//}
		licznik = 0;
        for (i = 3; i < Wsk_St_Wymiary->NI - 1; i += 2){
            Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j, k, 3)] = temp[licznik];
            licznik++;
        }
		free(A_wspol);
		//printf("licznik=%d iter=%d\n", licznik, iter);
//        printf("FLAGA_2.6\n");
	}
	
    printf("FLAGA_2.7\n");
		//<<<<<<<<<<<<<<<<<<<<
		//free(cp);
//	fclose(f);
		//free(temp);
	}

