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

const word FIXED_TYPE = "fixed";
const word INTERPOLATED_TYPE = "interpolated";

PatchFunction1Types::MappedFile<vector>* newVelocityMapper(const fvPatch& p, const dictionary& dict) {
    return new PatchFunction1Types::MappedFile<vector>(
                p.patch(),
                "velocityMapper",
                dict,
                "U",
                true
            );
}
PatchFunction1Types::MappedFile<vector>* cloneVelocityMapper(const fvPatch& p, const PatchFunction1Types::MappedFile<vector>& other) {
    return new PatchFunction1Types::MappedFile<vector>(other, p.patch());
}
PatchFunction1Types::MappedFile<vector>* cloneVelocityMapper(const PatchFunction1Types::MappedFile<vector>& other) {
    return new PatchFunction1Types::MappedFile<vector>(other);
}


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
    synTurb_(db(), p),
    refType_(FIXED_TYPE),
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
    refType_(ptf.refType_),
    curTimeIndex_(-1),
    synTurb_(ptf.synTurb_),
    corelate_(false)
{
    if(refType_ == INTERPOLATED_TYPE) {
        velocityMapper_.set(cloneVelocityMapper(p, ptf.velocityMapper_()));
    }
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
    referenceField_(p.size()),
    refType_(dict.lookupOrDefault("referenceType", FIXED_TYPE)),
    curTimeIndex_(-1),
    synTurb_(db(), dict, iF.mesh(), p),
    corelate_(false)
{
    vectorField& patchField = *this;

    if(refType_ == INTERPOLATED_TYPE) {
        velocityMapper_.set(newVelocityMapper(p, dict));
        referenceField_ = velocityMapper_->value(db().time().timeOutputValue());
    }
    else {
        referenceField_ = tmp<vectorField>(new vectorField("referenceField", dict, p.size()));
    }

    //synTurb_.setRefVelocity(average(mag(referenceField_)));

    if (dict.found("value")) {
        //line below might be unecessary
        //fixedValueFvPatchVectorField::operator=(vectorField("value", dict, p.size()));
        patchField = vectorField("value", dict, p.size());
        flucts_ = patchField - referenceField_;
        if(flucts_.size() > 0) {
            //Pout << "Avg flucts" << average(flucts_) << ", avg field " << average(patchField) <<", avg ref " << average(referenceField_) <<endl;
            if(average(mag(flucts_+vector(SMALL, SMALL, SMALL))) > 3*SMALL)
                corelate_ = true;
        }
    }
    else {
        synTurb_.computeNonuniformFlucts(referenceField_, db().time().timeOutputValue(), flucts_, false);
        fixedValueFvPatchVectorField::operator=(referenceField_ + flucts_);
    }
}


synTurbulenceInletFvPatchField::synTurbulenceInletFvPatchField
(
    const synTurbulenceInletFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    flucts_(ptf.flucts_),
    referenceField_(ptf.referenceField_),
    refType_(ptf.refType_),
    curTimeIndex_(-1),
    synTurb_(ptf.synTurb_),
    corelate_(false)
{
    if(refType_ == INTERPOLATED_TYPE) {
        velocityMapper_.set(cloneVelocityMapper(ptf.velocityMapper_()));
    }
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
    refType_(ptf.refType_),
    curTimeIndex_(-1),
    synTurb_(ptf.synTurb_),
    corelate_(false)
{
    if(refType_ == INTERPOLATED_TYPE) {
        velocityMapper_.set(cloneVelocityMapper(ptf.velocityMapper_()));
    }
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

    synTurb_.autoMap(m);
    //synTurb_.setRefVelocity(average(mag(referenceField_)));
    if(refType_ == INTERPOLATED_TYPE) {
        velocityMapper_().autoMap(m);
    }

    //Info <<"autoMap called" << endl;
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

    synTurb_.rmap(tiptf.synTurb_, addr);
    //synTurb_.setRefVelocity(average(mag(referenceField_ + VSMALL)));


    if(refType_ == INTERPOLATED_TYPE) {
        velocityMapper_().rmap(tiptf.velocityMapper_(), addr);
    }

    //Info <<" rmap called" << endl;
}

#include "OPstream.H"
#include "IPstream.H"
void synTurbulenceInletFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //Pout << "Updating boundary values of "<<this->patch().name()<<", ref field is " << average(mag(referenceField_)) <<" wtih size "<< size() <<", is my proc id " << Pstream::myProcNo() << endl;

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        vectorField& patchField = *this;

        if(refType_ == INTERPOLATED_TYPE) {
            referenceField_ = velocityMapper_->value( db().time().timeOutputValue() );
            //synTurb_.setRefVelocity(average(mag(referenceField_+ VSMALL)));
        }

        synTurb_.setTimeStep(db().time().deltaT().value());

        //synTurb_.computeNewFluctuations(flucts_, corelate_);
        synTurb_.computeNonuniformFlucts(referenceField_, db().time().timeOutputValue(), flucts_, corelate_);

//        Pout << "Fluctuations updated on " << patch().name() << endl;

        if( ! corelate_ )
            corelate_ = true;

        patchField = referenceField_ + flucts_;

        if(synTurb_.isPrintingStats()){
            Pout <<"Reference field(avg):"<< average(mag(referenceField_)) <<", Fluctuations (avg)" << average(mag(flucts_)) << endl;
        }

        curTimeIndex_ = this->db().time().timeIndex();

        fixedValueFvPatchVectorField::updateCoeffs();
    }
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
