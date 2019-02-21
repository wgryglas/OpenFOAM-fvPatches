#include "obliczenia.h"

void obliczanie_predkosci(short int indeks, unsigned int iter, unsigned int i, unsigned int j, unsigned k, struct Wymiary *Wsk_St_Wymiary, struct n_Predkosc *Wsk_St_tym_predkosc, struct s_Predkosc *Wsk_St_predkosc, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn, struct IOdata *Wsk_St_IOdata)
{
    double Fn=0, Fs=0, Fw=0, Fe=0, Ft=0, Fb=0, Dn=0, Ds=0, Dw=0, De=0, Dt=0, Db=0;
	double un, us, vw=0, ve=0, wt=0, wb=0;
	double Mout=0;
	double cisnienie=0;
	double *predkosc;
    unsigned int l, m;
	short int mod_a[6], mod_b[6], mod_c[6], mod_d[6], mod_e[6], mod_f[6];
    short unsigned int start, end;


	switch (indeks)
	{
    /*SWITCH HELP
     * MOD MODYFIKUJE OBECNE POLOZENIE DO POLOZENIA POSZUKIWANYCH WEKTOROW W KOMORKACH SASIEDNICH.
     * dla u
     * mod[0] -> un+/-1->uN & uP, mod[1] -> un+/-1->uS & uP, mod[2] -> vw->v(I-1,J-1) & v(I+1,J-1), mod[3] -> ve->v(I-1,J+1) & v(I+1,J+1), itp.
     * kolejne litery abc<-IJK uN; def<IJK uS itp.
     * UZYWANE W FUNCKJI translation
     * USTALA ODPOWIEDNIE CISNIENIE DLA PRAWEJ STRONY ->d ORAZ MASKUJE WARTOSC PREDKOSCI Z POP. KROKU JAKO predkosc <- DODAWANE DO d DLA TDMA
     */
	case 0:
		start = 2; end = Wsk_St_Wymiary->NI - 1;
		mod_a[0] = 2; mod_a[1] = -2; mod_a[2] = 1; mod_a[3] = 1; mod_a[4] = 1; mod_a[5] = 1;
		mod_b[0] = 0; mod_b[1] = 0; mod_b[2] = 1; mod_b[3] = -1; mod_b[4] = 0; mod_b[5] = 0;
		mod_c[0] = 0; mod_c[1] = 0; mod_c[2] = 0; mod_c[3] = 0; mod_c[4] = 1; mod_c[5] = -1;
		mod_d[0] = 0; mod_d[1] = 0; mod_d[2] = -1; mod_d[3] = -1; mod_d[4] = -1; mod_d[5] = -1;
		mod_e[0] = 0; mod_e[1] = 0; mod_e[2] = 1; mod_e[3] = -1; mod_e[4] = 0; mod_e[5] = 0;
		mod_f[0] = 0; mod_f[1] = 0; mod_f[2] = 0; mod_f[3] = 0; mod_f[4] = 1; mod_f[5] = -1;
		
        cisnienie = (Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i + 1, j, k, 3)] - Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i - 1, j, k, 3)]) *Wsk_St_Wymiary->Ax;
		predkosc = Wsk_St_predkosc->u;
		break;
	case 1:
		start = 1; end = Wsk_St_Wymiary->NI;
		mod_a[0] = 1; mod_a[1] = -1; mod_a[2] = 0; mod_a[3] = 0; mod_a[4] = 0; mod_a[5] = 0;
		mod_b[0] = -1; mod_b[1] = -1; mod_b[2] = 2; mod_b[3] = -2; mod_b[4] = -1; mod_b[5] = -1;
		mod_c[0] = 0; mod_c[1] = 0; mod_c[2] = 0; mod_c[3] = 0; mod_c[4] = 1; mod_c[5] = -1;
		mod_d[0] = 1; mod_d[1] = -1; mod_d[2] = 0; mod_d[3] = 0; mod_d[4] = 0; mod_d[5] = 0;
		mod_e[0] = 1; mod_e[1] = 1; mod_e[2] = 0; mod_e[3] = 0; mod_e[4] = 1; mod_e[5] = 1;
		mod_f[0] = 0; mod_f[1] = 0; mod_f[2] = 0; mod_f[3] = 0; mod_f[4] = 1; mod_f[5] = -1;
		cisnienie = (Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i, j + 1, k, 3)] - Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i, j - 1, k, 3)]) *Wsk_St_Wymiary->Ay;
		predkosc = Wsk_St_predkosc->v;
		break;
	case 2:
		start = 1; end = Wsk_St_Wymiary->NI;
		mod_a[0] = 1; mod_a[1] = -1; mod_a[2] = 0; mod_a[3] = 0; mod_a[4] = 0; mod_a[5] = 0;
		mod_b[0] = 0; mod_b[1] = 0; mod_b[2] = 1; mod_b[3] = -1; mod_b[4] = 0; mod_b[5] = 0;
		mod_c[0] = -1; mod_c[1] = -1; mod_c[2] = 1; mod_c[3] = 1; mod_c[4] = 2; mod_c[5] = -2;
		mod_d[0] = 1; mod_d[1] = -1; mod_d[2] = 0; mod_d[3] = 0; mod_d[4] = 0; mod_d[5] = 0;
		mod_e[0] = 0; mod_e[1] = 0; mod_e[2] = 1; mod_e[3] = -1; mod_e[4] = 0; mod_e[5] = 0;
		mod_f[0] = 1; mod_f[1] = 1; mod_f[2] = -1; mod_f[3] = -1; mod_f[4] = 0; mod_f[5] = 0;
		cisnienie = (Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i, j, k + 1, 3)] - Wsk_St_cisnienie->p[translation(Wsk_St_Wymiary, i, j, k - 1, 3)]) *Wsk_St_Wymiary->Az;
		predkosc = Wsk_St_predkosc->w;
		break;
	default:
		break;
	}


	if (i == start)
	{
        Wsk_St_tym_predkosc->ap[iter] = 1.0;
		Wsk_St_tym_predkosc->d[iter] = Wsk_St_tym_predkosc->d_ini[translation(Wsk_St_Wymiary, i, j, k, indeks + 4)];

	}
	else if (i == end){

        Wsk_St_tym_predkosc->ap[iter] = 1.0;
        if (indeks == 0){
//            for (l = 2; l < Wsk_St_Wymiary->NK; l += 2)
//            for (m = 2; m < Wsk_St_Wymiary->NJ; m += 2)
//            {
//                Mout += Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i, l, m, 0)];

//            }
//            if (Mout != 0){
//                Wsk_St_tym_predkosc->d[iter] = Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i - 2, j, k, 0)] * Wsk_St_Wymiary->Min / (Wsk_St_plyn->rho*Wsk_St_Wymiary->Ax*Mout);
//            }
//            else
                Wsk_St_tym_predkosc->d[iter] = Wsk_St_IOdata->predkosc;//Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i - 2, j, k, 0)];
        }
			if (indeks == 1 )
				Wsk_St_tym_predkosc->d[iter] = Wsk_St_predkosc->v[translation(Wsk_St_Wymiary, i - 2, j, k, 1)];
			if ( indeks == 2)
				Wsk_St_tym_predkosc->d[iter] = Wsk_St_predkosc->w[translation(Wsk_St_Wymiary, i - 2, j, k, 2)];
	}
	else{

        un = (Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i + mod_a[0], j + mod_b[0], k + mod_c[0], 0)] -
                Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i + mod_d[0], j + mod_e[0], k +mod_f[0], 0)]) / Wsk_St_Wymiary->Wymiar_Elem_X;

        us = (-Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i + mod_a[1], j + mod_b[1], k + mod_c[1], 0)] +
                Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i + mod_d[1], j + mod_e[1], k + mod_f[1], 0)]) / Wsk_St_Wymiary->Wymiar_Elem_X;

		if (Wsk_St_Wymiary->wymiar >1){		
            vw = (Wsk_St_predkosc->v[translation(Wsk_St_Wymiary, i + mod_a[2], j + mod_b[2], k + mod_c[2], 1)] +
                    Wsk_St_predkosc->v[translation(Wsk_St_Wymiary, i + mod_d[2], j + mod_e[2], k +
                    mod_f[2], 1)]) / 2;
            ve = (Wsk_St_predkosc->v[translation(Wsk_St_Wymiary, i + mod_a[3], j + mod_b[3], k + mod_c[3], 1)] +
                    Wsk_St_predkosc->v[translation(Wsk_St_Wymiary, i + mod_d[3], j + mod_e[3], k +
                    mod_f[3], 1)]) / 2;
		}
		if (Wsk_St_Wymiary->wymiar==3){
            wt = (Wsk_St_predkosc->w[translation(Wsk_St_Wymiary, i + mod_a[4], j + mod_b[4], k + mod_c[4], 2)] +
                    Wsk_St_predkosc->w[translation(Wsk_St_Wymiary, i + mod_d[4], j + mod_e[4],k + mod_f[4], 2)]) / 2;
            wb = (Wsk_St_predkosc->w[translation(Wsk_St_Wymiary, i + mod_a[5], j + mod_b[5], k + mod_c[5], 2)] +
                    Wsk_St_predkosc->w[translation(Wsk_St_Wymiary, i + mod_d[5], j + mod_e[5],k + mod_f[5], 2)]) / 2;
		}

		

			

		
//        Fn = Wsk_St_plyn->rho*Wsk_St_Wymiary->Ax*un;//Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i , j , k , 0)];
//        Fs = Wsk_St_plyn->rho*Wsk_St_Wymiary->Ax*us;//Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i-2, j, k, 0)];
//		Fw = Wsk_St_plyn->rho*Wsk_St_Wymiary->Ay*vw;
//        Fe = Wsk_St_plyn->rho*Wsk_St_Wymiary->Ay*ve;
//		Ft = Wsk_St_plyn->rho*Wsk_St_Wymiary->Az*wt;
//		Fb = Wsk_St_plyn->rho*Wsk_St_Wymiary->Az*wb;

        Dn = Wsk_St_Wymiary->Ax/pow(Wsk_St_Wymiary->Wymiar_Elem_X,2)*(Wsk_St_plyn->mi/Wsk_St_plyn->rho+pow(0.42,2)*pow((i)*Wsk_St_Wymiary->Wymiar_Elem_X,2)*un);
        //Wsk_St_plyn->mi * Wsk_St_Wymiary->Ax / Wsk_St_Wymiary->Wymiar_Elem_X;
        Ds = Wsk_St_Wymiary->Ax/pow(Wsk_St_Wymiary->Wymiar_Elem_X,2)*(Wsk_St_plyn->mi/Wsk_St_plyn->rho+pow(0.42,2)*pow((i-2)*Wsk_St_Wymiary->Wymiar_Elem_X,2)*us);
        //Wsk_St_plyn->mi * Wsk_St_Wymiary->Ax / Wsk_St_Wymiary->Wymiar_Elem_X;
//        Dw =  Wsk_St_plyn->mi * Wsk_St_Wymiary->Ay / Wsk_St_Wymiary->Wymiar_Elem_Y;
//        De =  Wsk_St_plyn->mi * Wsk_St_Wymiary->Ay / Wsk_St_Wymiary->Wymiar_Elem_Y;
//        Dt =  Wsk_St_plyn->mi * Wsk_St_Wymiary->Az / Wsk_St_Wymiary->Wymiar_Elem_Z;
//        Db =  Wsk_St_plyn->mi * Wsk_St_Wymiary->Az / Wsk_St_Wymiary->Wymiar_Elem_Z;
//        printf("I=%d\tuN=%.8f\tDn=%.8f\n",i,un,Dn);
//        printf("I=%d\tus=%.8f\tDs=%.8f\n",i,us,Ds);


		if (k == 2 ){
            Wsk_St_tym_predkosc->ab[iter] =Db*2;
		}
        else{			
            Wsk_St_tym_predkosc->ab[iter] = Db + Fb / 2;
		}
		if (k == Wsk_St_Wymiary->NK - 1 ){			
            Wsk_St_tym_predkosc->at[iter] = Dt*2;
		}
		else{
            Wsk_St_tym_predkosc->at[iter] = Dt - Ft / 2;
		}

		if (j == 2){
            Wsk_St_tym_predkosc->ae[iter] = De*2;
		}
		else{
            Wsk_St_tym_predkosc->ae[iter] = De + Fe;
		}
        if (j == Wsk_St_Wymiary->NJ - 1){
            Wsk_St_tym_predkosc->aw[iter] = Dw*2;
		}
		else{
            Wsk_St_tym_predkosc->aw[iter] = Dw - Fw / 2;
		}

        Wsk_St_tym_predkosc->as[iter] = (Ds);// + Fs);
        Wsk_St_tym_predkosc->an[iter] = (Dn);// - Fn / 2);

        Wsk_St_tym_predkosc->ap[iter] = -Wsk_St_tym_predkosc->an[iter] - Wsk_St_tym_predkosc->as[iter] +
                                         Wsk_St_tym_predkosc->aw[iter] + Wsk_St_tym_predkosc->ae[iter] +
                                         Wsk_St_tym_predkosc->at[iter] + Wsk_St_tym_predkosc->ab[iter] +
                                         Fn - Fs + Fw - Fe + Ft - Fb;

        Wsk_St_tym_predkosc->d[iter] = Wsk_St_IOdata->grad_cisnienia;

//    PRAWA STRONA
        if (Wsk_St_Wymiary->wymiar > 1)
            Wsk_St_tym_predkosc->d[iter] += Wsk_St_tym_predkosc->aw[iter] * predkosc[translation(Wsk_St_Wymiary, i, j + 2, k, indeks)] +
            Wsk_St_tym_predkosc->ae[iter] * predkosc[translation(Wsk_St_Wymiary, i, j - 2, k, indeks)];
		
		if (Wsk_St_Wymiary->wymiar==3)
			Wsk_St_tym_predkosc->d[iter] += Wsk_St_tym_predkosc->at[iter] * predkosc[translation(Wsk_St_Wymiary, i, j, k + 2, indeks)] +
			Wsk_St_tym_predkosc->ab[iter] * predkosc[translation(Wsk_St_Wymiary, i, j, k - 2, indeks)];
		
	}




}
void obliczanie_cisnienia( unsigned int iter, unsigned int i, unsigned int j, unsigned k, struct Wymiary *Wsk_St_Wymiary, struct s_Predkosc *Wsk_St_predkosc, struct Cisnienie *Wsk_St_cisnienie, struct Plyn *Wsk_St_plyn)
{
	if (i==1)	
	{
		Wsk_St_cisnienie->ap[iter] = 1;
		Wsk_St_cisnienie->d[iter] = 0;
	}
	else if (i==Wsk_St_Wymiary->NI)
	{
		Wsk_St_cisnienie->ap[iter] = 1;
		Wsk_St_cisnienie->d[iter] = 0;

	}
    else{
        if (k == 2){
			Wsk_St_cisnienie->ab[iter] = 0;
		}
        else{
			Wsk_St_cisnienie->ab[iter] = Wsk_St_plyn->rho*Wsk_St_Wymiary->Az*(Wsk_St_Wymiary->Az / Wsk_St_predkosc->ap_w[translation(Wsk_St_Wymiary, i, j, k - 1, 9)]);
			//fprintf(f, "\ti=%d\tj=%d\tk=%d\titer=%d\tdp=%f\tcp=\n", i, j, k-1, iter, Wsk_St_predkosc->ap_w[translation(Wsk_St_Wymiary, i, j, k - 1, 9)]);
		}
        if (k == Wsk_St_Wymiary->NK - 1){
			Wsk_St_cisnienie->at[iter] = 0;
		}
        else{
			if (Wsk_St_Wymiary->wymiar==3)
			Wsk_St_cisnienie->at[iter] = Wsk_St_plyn->rho*Wsk_St_Wymiary->Az*(Wsk_St_Wymiary->Az / Wsk_St_predkosc->ap_w[translation(Wsk_St_Wymiary, i, j, k + 1, 9)]);			
		}

        if (j == 2){
			Wsk_St_cisnienie->ae[iter] = 0;
		}
        else{
			Wsk_St_cisnienie->ae[iter] = Wsk_St_plyn->rho*Wsk_St_Wymiary->Ay*(Wsk_St_Wymiary->Ay / Wsk_St_predkosc->ap_v[translation(Wsk_St_Wymiary, i, j - 1, k, 8)]);
			
		}
        if (j == Wsk_St_Wymiary->NJ - 1){
			Wsk_St_cisnienie->aw[iter] = 0;
		}
        else{
			Wsk_St_cisnienie->aw[iter] = Wsk_St_plyn->rho*Wsk_St_Wymiary->Ay*(Wsk_St_Wymiary->Ay / Wsk_St_predkosc->ap_v[translation(Wsk_St_Wymiary, i, j + 1, k, 8)]);
			
		}

        if (i == 3){
            Wsk_St_cisnienie->as[iter] = 0;
        }
        else{
             Wsk_St_cisnienie->as[iter] =-(Wsk_St_plyn->rho*Wsk_St_Wymiary->Ax*(Wsk_St_Wymiary->Ax / Wsk_St_predkosc->ap_u[translation(Wsk_St_Wymiary, i - 1, j, k, 7)]));

        }
        if (j == Wsk_St_Wymiary->NI - 2){
            Wsk_St_cisnienie->an[iter] = 0;
        }
        else{
            Wsk_St_cisnienie->an[iter] =- (Wsk_St_plyn->rho*Wsk_St_Wymiary->Ax*(Wsk_St_Wymiary->Ax / Wsk_St_predkosc->ap_u[translation(Wsk_St_Wymiary, i + 1, j, k, 7)]));

        }




        Wsk_St_cisnienie->ap[iter] =- Wsk_St_cisnienie->an[iter] - Wsk_St_cisnienie->as[iter] +
                                      Wsk_St_cisnienie->aw[iter] + Wsk_St_cisnienie->ae[iter] +
                                      Wsk_St_cisnienie->at[iter] + Wsk_St_cisnienie->ab[iter];

        Wsk_St_cisnienie->d[iter] = Wsk_St_plyn->rho*Wsk_St_Wymiary->Ax*(Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i - 1, j, k, 0)] -
                                                                         Wsk_St_predkosc->u[translation(Wsk_St_Wymiary, i + 1, j, k, 0)]);


		if (Wsk_St_Wymiary->wymiar > 1){
			
            Wsk_St_cisnienie->d[iter] += Wsk_St_plyn->rho*Wsk_St_Wymiary->Ay*(Wsk_St_predkosc->v[translation(Wsk_St_Wymiary, i, j - 1, k, 1)] -
                                                                              Wsk_St_predkosc->v[translation(Wsk_St_Wymiary, i, j + 1, k, 1)]);
            Wsk_St_cisnienie->d[iter] += Wsk_St_cisnienie->aw[iter] * Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j + 2, k, 3)] +
                                         Wsk_St_cisnienie->ae[iter] * Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j - 2, k, 3)];
		}
			
		if (Wsk_St_Wymiary->wymiar==3){
            Wsk_St_cisnienie->d[iter] += Wsk_St_plyn->rho*Wsk_St_Wymiary->Az*(Wsk_St_predkosc->w[translation(Wsk_St_Wymiary, i, j, k - 1, 2)] -
                                                                              Wsk_St_predkosc->w[translation(Wsk_St_Wymiary, i, j, k + 1, 2)]);
            Wsk_St_cisnienie->d[iter] += Wsk_St_cisnienie->at[iter] * Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j, k + 2, 3)] +
                                         Wsk_St_cisnienie->ab[iter] * Wsk_St_cisnienie->p_correction[translation(Wsk_St_Wymiary, i, j, k - 2, 3)];
		}
    }
}
