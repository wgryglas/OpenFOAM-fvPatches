/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "budynekFvPatchField.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::budynekFvPatchField::budynekFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchVectorField(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{
   // initFraction(dir);
}


Foam::budynekFvPatchField::budynekFvPatchField
(
    const budynekFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchVectorField(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)
{
    if (&iF && mapper.hasUnmapped())
    {
        WarningIn
        (
            "budynekFvPatchField<Type>::budynekFvPatchField\n"
            "(\n"
            "    const budynekFvPatchField<Type>&,\n"
            "    const fvPatch&,\n"
            "    const DimensionedField<Type, volMesh>&,\n"
            "    const fvPatchFieldMapper&\n"
            ")\n"
        )   << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


Foam::budynekFvPatchField::budynekFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchVectorField(p, iF, dict),
    refValue_("refValue", dict, p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{

    dict.lookup("direction") >> dir;
    initFraction(dir);
    evaluate();
}



Foam::budynekFvPatchField::budynekFvPatchField
(
    const budynekFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchVectorField(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{ }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::budynekFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchVectorField::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    valueFraction_.autoMap(m);
}



void Foam::budynekFvPatchField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    transformFvPatchVectorField::rmap(ptf, addr);

    const budynekFvPatchField& dmptf =
        refCast<const budynekFvPatchField>(ptf);

    refValue_.rmap(dmptf.refValue_, addr);
    refGrad_.rmap(dmptf.refGrad_, addr);
    valueFraction_.rmap(dmptf.valueFraction_, addr);
}



Foam::tmp<Foam::vectorField >
Foam::budynekFvPatchField::snGrad() const
{
    const Field<vector> pif(this->patchInternalField());

    tmp<Field<vector> > normalValue = transform(valueFraction_, refValue_);

    tmp<Field<vector> > gradValue = pif + refGrad_/this->patch().deltaCoeffs();

    tmp<Field<vector> > transformGradValue =
        transform(I - valueFraction_, gradValue);

    return
        (normalValue + transformGradValue - pif)*
        this->patch().deltaCoeffs();
}


void Foam::budynekFvPatchField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    updateBoundaryValues();

    tmp<Field<vector> > normalValue = transform(valueFraction_, refValue_);

    tmp<Field<vector> > gradValue =
        this->patchInternalField() + refGrad_/this->patch().deltaCoeffs();

    tmp<Field<vector> > transformGradValue =
        transform(I - valueFraction_, gradValue);


    Field<vector>::operator=(normalValue + transformGradValue);


    transformFvPatchVectorField::evaluate();

}


Foam::tmp<Foam::vectorField >
Foam::budynekFvPatchField::snGradTransformDiag() const
{
    vectorField diag(valueFraction_.size());


    diag.replace
    (
        vector::X,
        sqrt(mag(valueFraction_.component(symmTensor::XX)))
    );
    diag.replace
    (
        vector::Y,
        sqrt(mag(valueFraction_.component(symmTensor::YY)))
    );
    diag.replace
    (
        vector::Z,
        sqrt(mag(valueFraction_.component(symmTensor::ZZ)))
    );

 return transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
}


void Foam::budynekFvPatchField::updateBoundaryValues()
{
    struct IOdata *Wsk_St_IOdata = (struct IOdata*)malloc(sizeof(struct IOdata));

  //  std::cout <<"OBLICZENIA MOJE\n" << std::endl;
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const volScalarField& p = db().lookupObject<volScalarField>("p");
    const volVectorField& U = db().lookupObject<volVectorField>("U");
    const labelUList& fc = patch().faceCells();
    const scalarField f_cell_coeff = mesh.boundary()[patch().index()].deltaCoeffs();

            //    int coos = patch().boundaryMesh().findPatchID("inlet");
            //    std::cout <<"COSSSSSSSSSSSSS" << coos<< std::endl;
            //   const scalarField f_cell_coeff2 = mesh.boundary()[coos].deltaCoeffs();
            //   double AA=0;
            //   double BB=0;
            //    forAll(f_cell_coeff2,cellID)
            //    {
            //        AA=AA+U[fc[cellID]].x();
            //        BB++;
            //    }
            //    std::cout <<"COS2S=" << AA<< std::endl;
            //    AA=AA/BB;
            //    std::cout <<"COSSSSS3=" << AA<< std::endl;
            //    std::cout <<"COSSSSSSSSS4S=" << BB<< std::endl;
    //FatalIOError.exit();
    



// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//volVectorField gradp(fvc::grad(p)); <---- To zostalo zamienione na to cos ponizej. "grad(p)" 
//                                          oznacza tylko nazwę pola pod jakim zostanie zcashowany 
//                                          gradient. OF ma taką możliwość, że jak jakaś operacja 
//                                          zostanie wykonana to ją cachuje i jeśli znów gdzieś ktoś
//                                          chce takie pole policzyć to on z pamięci je przepisuje. 
//                                          W tym wypadku to nie wiem czy oni dokładnie pod taką nazwą cachują
//                                          ale z naszego punktu widzenia nie musi być to poprawnie zrobione. 
    volVectorField gradp(  fv::gaussGrad<scalar>(p.mesh()).calcGrad(p,"grad(p)") );
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    forAll(f_cell_coeff,cellID)
    {
        Wsk_St_IOdata->cell_dist = 1/f_cell_coeff[cellID];
        Wsk_St_IOdata->grad_cisnienia = gradp.boundaryField().boundaryInternalField()[patch().index()][cellID][0];
        Wsk_St_IOdata->predkosc = U[fc[cellID]].x();

        wall_function(Wsk_St_IOdata);
        refGrad_[cellID].x()=Wsk_St_IOdata->grad_ut;
        // FatalIOError.exit();
    }



    //TODO
   // std::cout<<"wielkosc gradP  " <<gradp.boundaryField().boundaryInternalField()[patch().index()].size()<< std::endl;
   // std::cout<<"valfract 1"<<valueFraction_[1].yy()<< std::endl;
//    patch().boundaryMesh()
//     std::cout<< meshObject[0].info()<<std::endl;
   // std::cout<<refGrad_[0].x()<<std::endl;

}

void Foam::budynekFvPatchField::initFraction(Foam::label dir)
{

    SymmTensor<scalar> t;
    switch (dir)
    {
    case 0: //x direction refvalue
        t.xx()=1;   t.yy()=0;   t.zz()=0;   t.xy()=0;   t.xz()=0;    t.yz()=0;
        break;
    case 1:
        t.xx()=0;   t.yy()=1;   t.zz()=0;   t.xy()=0;   t.xz()=0;    t.yz()=0;
        break;
    case 2:
        t.xx()=0;   t.yy()=0;   t.zz()=1;   t.xy()=0;   t.xz()=0;    t.yz()=0;
        break;
    default:
        FatalErrorIn("budynekFvPatchField::initFraction")
                <<"allowed only main directions, 1,2,3 are respectivly x,y,z"
                <<abort(FatalError);
        break;
    }

   valueFraction_ = t;
}



void Foam::budynekFvPatchField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
   os.writeKeyword("direction") << dir << token::END_STATEMENT << nl;
    refValue_.writeEntry("refValue",os);
}


// ************************************************************************* //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        budynekFvPatchField
    )
}
