    /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::budynekFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Base class for direction-mixed boundary conditions.

SourceFiles
    budynekFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef budynekFvPatchField_H
#define budynekFvPatchField_H

#include "fvPatchFields.H"
#include "transformFvPatchFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include "wall_f/struktury.h"
#include "wall_f/wall_f.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class budynekFvPatchField Declaration
\*---------------------------------------------------------------------------*/
class budynekFvPatchField
:
    public transformFvPatchVectorField//transformFvPatchField<vector>
{
    // Private data

        //- Value field
        vectorField refValue_;

        //- Normal gradient field
        vectorField refGrad_;

        //- Fraction (0-1) of value used for boundary condition
        symmTensorField valueFraction_;

        label dir;


public:

    //- Runtime type information
    TypeName("budynek")


    // Constructors

        //- Construct from patch and internal field
        budynekFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        budynekFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given budynekFvPatchField onto
        //  a new patch
        budynekFvPatchField
        (
            const budynekFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new budynekFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        budynekFvPatchField
        (
            const budynekFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField >
            (
                new budynekFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
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


        // Return defining fields

            virtual Field<vector>& refValue()
            {
                return refValue_;
            }

            virtual const Field<vector>& refValue() const
            {
                return refValue_;
            }

            virtual Field<vector>& refGrad()
            {
                return refGrad_;
            }

            virtual const Field<vector>& refGrad() const
            {
                return refGrad_;
            }

            virtual symmTensorField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const symmTensorField& valueFraction() const
            {
                return valueFraction_;
            }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<vector> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking

            );

            //- Return face-gradient transform diagonal
            virtual tmp<Field<vector> > snGradTransformDiag() const;

            void updateBoundaryValues();

            void initFraction(label dir);


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchVectorField&) {}
        virtual void operator+=(const fvPatchVectorField&) {}
        virtual void operator-=(const fvPatchVectorField&) {}
        virtual void operator*=(const fvPatchVectorField&) {}
        virtual void operator/=(const fvPatchVectorField&) {}

        virtual void operator=(const Field<vector>&) {}
        virtual void operator+=(const Field<vector>&) {}
        virtual void operator-=(const Field<vector>&) {}
        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const vector&) {}
        virtual void operator+=(const vector&) {}
        virtual void operator-=(const vector&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
