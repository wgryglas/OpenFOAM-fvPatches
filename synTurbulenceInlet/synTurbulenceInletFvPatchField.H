/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::turbulentInletFvPatchField

Group
    grpInletBoundaryConditions

Description
    This boundary condition generates a fluctuating inlet condition by adding
    a random component to a reference (mean) field.

    \f[
        x_p = (1 - \alpha) x_p^{n-1} + \alpha (x_{ref} + s C_{RMS} x_{ref})
    \f]

    where

    \vartable
        x_p     | patch values
        x_{ref} | reference patch values
        n       | time level
        \alpha  | fraction of new random component added to previous time value
        C_{RMS} | RMS coefficient
        s       | fluctuation scale
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        fluctuationScale | RMS fluctuation scale (fraction of mean) | yes |
        referenceField | reference (mean) field | yes        |
        alpha | fraction of new random component added to previous| no| 0.1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            turbulentInlet;
        fluctuationScale 0.1;
        referenceField  uniform 10;
        alpha           0.1;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    turbulentInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef synTurbulenceInletFvPatchField_H
#define synTurbulenceInletFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "randomAngles.H"
#include "synTurbulence.H"
#include "vector.H"
#include "MappedFile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class synTurbulenceInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class synTurbulenceInletFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Fluctuation fields
        vectorField flucts_;

        //- Reference field
        vectorField referenceField_;


        //- Switch for type of refernece field:  fixed | interpolated
        word refType_;
        autoPtr<PatchFunction1Types::MappedFile<vector> > velocityMapper_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- data structure holding information about syntetic turbulence
        // end evaluating velocity fluctuation field
        synTurbulence synTurb_;
        bool corelate_;

public:

    //- Runtime type information
    TypeName("synTurbulenceInlet");


    // Constructors

        //- Construct from patch and internal field
        synTurbulenceInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        synTurbulenceInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given synTurbulenceInletFvPatchField
        //  onto a new patch
        synTurbulenceInletFvPatchField
        (
            const synTurbulenceInletFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        synTurbulenceInletFvPatchField
        (
            const synTurbulenceInletFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new synTurbulenceInletFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        synTurbulenceInletFvPatchField
        (
            const synTurbulenceInletFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new synTurbulenceInletFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the fluctuation scale
            const vectorField& fluctuationField() const
            {
                return flucts_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            vectorField& fluctuationField()
            {
                return flucts_;
            }

            //- Return the reference field
            const vectorField& referenceField() const
            {
                return referenceField_;
            }

            //- Return reference to the reference field to allow adjustment
            vectorField& referenceField()
            {
                return referenceField_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif //synTurbulenceInletFvPatchField_H

// ************************************************************************* //
