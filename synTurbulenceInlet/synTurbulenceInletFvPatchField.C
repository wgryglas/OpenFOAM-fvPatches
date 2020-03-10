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

#include "synTurbulenceInletFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "Time.H"
#include "Pstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include <string>
namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
synTurbulenceInletFvPatchField::synTurbulenceInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    flucts_(p.size(),vector::zero),
    referenceField_(p.size()),
    curTimeIndex_(-1),
    synTurb_(db()),
    corelate_(false)
{}


synTurbulenceInletFvPatchField::synTurbulenceInletFvPatchField
(
    const synTurbulenceInletFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flucts_(ptf.flucts_,mapper),
    referenceField_(ptf.referenceField_, mapper),
    curTimeIndex_(-1),
    synTurb_(ptf.synTurb_),
    corelate_(false)
{
}


synTurbulenceInletFvPatchField::synTurbulenceInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    flucts_(p.size(),vector::zero),
    referenceField_("referenceField", dict, p.size()),
    curTimeIndex_(-1),
    synTurb_(db(), dict),
    corelate_(false)
{
    if (dict.found("value")) {
        flucts_ = (*this) - referenceField_;
        corelate_=true;
    }
    synTurb_.setRefVelocity(average(mag(referenceField_)));
}


synTurbulenceInletFvPatchField::synTurbulenceInletFvPatchField
(
    const synTurbulenceInletFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    flucts_(ptf.flucts_),
    referenceField_(ptf.referenceField_),
    curTimeIndex_(-1),
    synTurb_(ptf.synTurb_),
    corelate_(false)
{
}


synTurbulenceInletFvPatchField::synTurbulenceInletFvPatchField
(
    const synTurbulenceInletFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    flucts_(ptf.flucts_),
    referenceField_(ptf.referenceField_),
    curTimeIndex_(-1),
    synTurb_(ptf.synTurb_),
    corelate_(false)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void synTurbulenceInletFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    flucts_.autoMap(m);
    referenceField_.autoMap(m);
    synTurb_.setRefVelocity(average(mag(referenceField_)));
}


void synTurbulenceInletFvPatchField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const synTurbulenceInletFvPatchField& tiptf =
        refCast<const synTurbulenceInletFvPatchField>(ptf);

    flucts_.rmap(tiptf.flucts_, addr);
    referenceField_.rmap(tiptf.referenceField_, addr);
    synTurb_.setRefVelocity(average(mag(referenceField_)));
}


void synTurbulenceInletFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        vectorField& patchField =*this;
        const vectorField& faceCenters = this->patch().Cf();

        synTurb_.setTimeStep(db().time().deltaT().value());
        synTurb_.computeNewFluctuations(faceCenters, flucts_, corelate_);

        if(!corelate_)
            corelate_=true;

        patchField = referenceField_ + flucts_;

//        vectorField& patchField = *this;

//        vectorField randomField(this->size());

//        forAll(patchField, facei)
//        {
//            ranGen_.randomise(randomField[facei]);
//        }

//        // Correction-factor to compensate for the loss of RMS fluctuation
//        // due to the temporal correlation introduced by the alpha parameter.
//        scalar rmsCorr = sqrt(12*(2*alpha_ - sqr(alpha_)))/alpha_;

//        patchField =
//            (1 - alpha_)*patchField
//          + alpha_*
//            (
//                referenceField_
//              + rmsCorr*cmptMultiply
//                (
//                    randomField - 0.5*vector::one,
//                    fluctuationScale_
//                )*mag(referenceField_)
//            );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void synTurbulenceInletFvPatchField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    referenceField_.writeEntry("referenceField", os);
    synTurb_.write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        synTurbulenceInletFvPatchField
    );
}
