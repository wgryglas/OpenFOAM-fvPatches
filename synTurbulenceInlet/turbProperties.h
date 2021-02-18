#ifndef TURBPROPERTIES_H
#define TURBPROPERTIES_H

#include "scalar.H"
#include "vector.H"
#include "vectorField.H"
#include "scalarField.H"
#include "fvMesh.H"
#include "constants.H"
#include "IOobject.H"
#include "tmp.H"
#include "fvPatchField.H"
#include "MappedFile.H"

namespace Foam {

    namespace TurbProperties {
        extern const scalar Cmu;
        extern const word TIME_SCALE_PROP_NAME;
    }

    template<typename T>
    PatchFunction1Types::MappedFile<T>* newMapper(word name, const fvPatch& p, const dictionary& dict)  {
        return new PatchFunction1Types::MappedFile<T>(
                    p.patch(),
                    name + "Mapper",
                    dict,
                    name,
                    true
                );
    }

    template<typename T>
    PatchFunction1Types::MappedFile<T>* cloneMapper(const fvPatch& p, const PatchFunction1Types::MappedFile<T>& other) {
        return new PatchFunction1Types::MappedFile<T>(other, p.patch());
    }

    template<typename T>
    PatchFunction1Types::MappedFile<T>* cloneMapper(const PatchFunction1Types::MappedFile<T>& other) {
        return new PatchFunction1Types::MappedFile<T>(other);
    }

    class synTurbulenceInletFvPatchField;

    struct synTurbulenceParameters: public Foam::tmp<synTurbulenceParameters>::refCount {
        virtual void update(const vectorField& refVelocity, const scalar& timeValue) = 0;
        virtual const scalarField& getTurbLengthScales() const  = 0;
        virtual const scalarField& getUrms() const  = 0;
        virtual const scalarField& getDissipationRates() const  = 0;
        virtual const scalarField& getTimeScales() const = 0;
        virtual ~synTurbulenceParameters() {}
        virtual void write(Ostream& os) const = 0;
        virtual synTurbulenceParameters* clone(const fvPatch& patch) const = 0;

        virtual word type() const = 0;
        virtual void autoMap(const fvPatchFieldMapper&) = 0;
        virtual void rmap(const synTurbulenceParameters&, const labelList&) = 0;
    };
}

#endif // TURBPROPERTIES_H
