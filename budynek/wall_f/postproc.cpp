#include "postproc.h"

void postproc(struct IOdata *Wsk_St_IOdata,struct Wymiary *Wsk_St_Wymiary,struct s_Predkosc *Wsk_St_predkosc, struct Plyn *Wsk_St_plyn )
{

    double shearstress;
    double ushear;
    double grad_wall;



    grad_wall=(Wsk_St_predkosc->u[translation(Wsk_St_Wymiary,4,2,2,0)]-Wsk_St_predkosc->u[translation(Wsk_St_Wymiary,2,2,2,0)])/Wsk_St_Wymiary->Wymiar_Elem_X;
//    shearstress=fabs();
    ushear=sqrt(shearstress/Wsk_St_plyn->rho);
    Wsk_St_IOdata->grad_ut=ushear/Wsk_St_plyn->mi*Wsk_St_plyn->rho*Wsk_St_Wymiary->Wymiar_Elem_Y/2;


}
