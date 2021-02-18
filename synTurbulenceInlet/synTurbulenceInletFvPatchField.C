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
#include "turbProperties.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include <string>
namespace Foam
{

const word FIXED_TYPE = "fixed";
const word INTERPOLATED_TYPE = "interpolated";



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
{
//    Pout <<"(p, iF) const. called" << endl;
}

//TODO fix synTurb constructor
synTurbulenceInletFvPatchField::synTurbulenceInletFvPatchField
(
    const synTurbulenceInletFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flucts_(ptf.flucts_, mapper),
    referenceField_(ptf.referenceField_, mapper),
    refType_(ptf.refType_),
    curTimeIndex_(-1),
    synTurb_(ptf.synTurb_),
    corelate_(false)
{
    if(refType_ == INTERPOLATED_TYPE) {
        velocityMapper_.set(cloneMapper(p, ptf.velocityMapper_()));
//        referenceField_ = velocityMapper_->value(db().time().timeOutputValue());
        //vectorField& patchField = *this;
//        patchField=(referenceField_ + flucts_);
    }
//    Pout <<"(ptf, p, iF, mapper) const. called" << endl;
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
        velocityMapper_.set(newMapper<vector>("U", p, dict));
        referenceField_ = velocityMapper_->value(db().time().timeOutputValue());
    }
    else {
        referenceField_ = vectorField("referenceField", dict, p.size());
    }

    //synTurb_.setRefVelocity(average(mag(referenceField_)));

    if (dict.found("value") && patch().size() > 0) {
        //line below might be unecessary
        //fixedValueFvPatchVectorField::operator=(vectorField("value", dict, p.size()));
        patchField = vectorField("value", dict, p.size());
        flucts_ = patchField - referenceField_;

        Pout << "(" << patch().name() << ")" << "Avg flucts" << average(flucts_) << ", avg field " << average(patchField) <<", avg ref " << average(referenceField_) <<endl;
        if(max(mag(flucts_+vector(SMALL, SMALL, SMALL))) > sqrt(3.0)*SMALL + SMALL)
            corelate_ = true;
    }
    else {
        //below causes huge flucts and destroys solution while single BC is split onto more than 1 processor
        //synTurb_.computeNonuniformFlucts(referenceField_, db().time().timeOutputValue(), flucts_, false);
        fixedValueFvPatchVectorField::operator=(referenceField_);
    }

//    Pout <<"(p, iF, dict) const. called" << endl;
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
        velocityMapper_.set(cloneMapper(ptf.velocityMapper_()));
    }
//    Pout <<"(ptf) const. called" << endl;
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
        velocityMapper_.set(cloneMapper(ptf.velocityMapper_()));
    }

//    Pout <<"(ptf, iF) const. called" << endl;
//    if(patch().size()> 0)
//        Pout <<"("<<patch().name()<<"="<<patch().size()<<")"<<"Setting flucts: " << min(flucts_) << max(flucts_) << endl;
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

//    Pout <<"autoMap called" << endl;
}


void synTurbulenceInletFvPatchField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const synTurbulenceInletFvPatchField& tiptf = refCast<const synTurbulenceInletFvPatchField>(ptf);

    flucts_.rmap(tiptf.flucts_, addr);
    referenceField_.rmap(tiptf.referenceField_, addr);

//    synTurb_.rmap(tiptf.synTurb_, addr);
//    //synTurb_.setRefVelocity(average(mag(referenceField_ + VSMALL)));

//    if(refType_ == INTERPOLATED_TYPE && tiptf.velocityMapper_.valid() && velocityMapper_.valid()) {
//        velocityMapper_().rmap(tiptf.velocityMapper_(), addr);
//    }

//    Info <<"(" << patch().name() << "=" << patch().size() <<")" <<" rmap called" << endl;
}

#include "OPstream.H"
#include "IPstream.H"
void synTurbulenceInletFvPatchField::updateCoeffs()
{
    if (this->updated()) {
        return;
    }

//    if(db().time().timeIndex() == 1) {
//        flucts_ = vector::zero;
//    }
//    else {
//        Pout <<" time: " << db().time().timeOutputValue() << "vs " << db().time().value() << endl;
//    }

    //Pout << "Updating boundary values of "<<this->patch().name()<<", ref field is " << average(mag(referenceField_)) <<" wtih size "<< size() <<", is my proc id " << Pstream::myProcNo() << endl;

    if (curTimeIndex_ != this->db().time().timeIndex()) {

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

        if(synTurb_.isPrintingStats()) {
            Pout <<"Reference field(avg):"<< average(mag(referenceField_)) <<", Fluctuations (avg)" << average(mag(flucts_)) << endl;
        }

        curTimeIndex_ = this->db().time().timeIndex();

        fixedValueFvPatchVectorField::updateCoeffs();
    }
}


void synTurbulenceInletFvPatchField::write(Ostream& os) const {
    fvPatchVectorField::write(os);
    if(refType_ == INTERPOLATED_TYPE) {
        os.writeEntry("referenceType", INTERPOLATED_TYPE);
    }
    else {
        referenceField_.writeEntry("referenceField", os);
    }
    synTurb_.write(os);
    //const vectorField& patchField = *this;
    vectorField current = referenceField_ + flucts_;
    current.writeEntry("value", os);
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
