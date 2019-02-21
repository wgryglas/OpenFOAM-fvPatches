#include "wall_f.h"

void wall_function(struct IOdata *Wsk_St_IOdata)
{

	int i, k, j;

	struct Wymiary *Wsk_St_Wymiary = (struct Wymiary*)malloc(sizeof(struct Wymiary));
      if (Wsk_St_Wymiary == NULL) ERROR_ALLOC;


    time_t start,end;


    czytaj_rozmiar(Wsk_St_Wymiary,Wsk_St_IOdata);
//            printf("NI=%d\tNk=%d\tNJ=%d\n ii %d", Wsk_St_Wymiary->NI, Wsk_St_Wymiary->NJ, Wsk_St_Wymiary->NK);

	struct Plyn *Wsk_St_plyn = alokacja_plyn(Wsk_St_Wymiary);
    struct s_Predkosc *Wsk_St_s_predkosc = alokacja_s_p(Wsk_St_Wymiary);
    struct Cisnienie *Wsk_St_cisnienie = alokacja_c(Wsk_St_Wymiary);

    //Wsk_St_IOdata->grad_cisnienia=0;//0.0027324;
   // Wsk_St_IOdata->predkosc=0;// 0.25;

    //inicjalizacja
        inicjalizacja(Wsk_St_Wymiary, Wsk_St_plyn, Wsk_St_s_predkosc,Wsk_St_cisnienie, Wsk_St_IOdata);

//            for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
//            for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
//            for (i = 2; i < Wsk_St_Wymiary->NI; i += 2)
//                printf("\t\t\tU_START=%f\n", Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)]);

    //obliczenia
    start = time(NULL);
        iteracja(Wsk_St_Wymiary, Wsk_St_s_predkosc, Wsk_St_cisnienie, Wsk_St_plyn,Wsk_St_IOdata);
    end = time(NULL);
//        postporc(Wsk_St_IOdata,Wsk_St_Wymiary,Wsk_St_s_predkosc,Wsk_St_plyn);
//    printf("The loop used %f seconds.\n", difftime(end, start));

//        for (k = 2; k < Wsk_St_Wymiary->NK; k += 2)
//		for (j = 2; j < Wsk_St_Wymiary->NJ; j += 2)
//        for (i = 2; i < Wsk_St_Wymiary->NI; i += 2)

//                printf("\t\t\tUPOLICZONE=%f\n", Wsk_St_s_predkosc->u[translation(Wsk_St_Wymiary, i, j, k, 0)]);





}







