#ifndef FIXEDTURBPROPERTIES_H
#define FIXEDTURBPROPERTIES_H

#include "turbProperties.h"
#include "fvPatch.H"

namespace Foam {

    class FixedTurbProperties : public synTurbulenceParameters
    {
        scalar m_tls;
        scalar m_ti;

        scalarField f_tls;
        scalarField f_umrs;
        scalarField f_eps;
        scalarField f_tts;

    public:
        static const word typeName;
        word type() const;

        FixedTurbProperties(const dictionary& dict, const fvPatch& patch);
        FixedTurbProperties(const scalar& intensity, const scalar& lengthScale, const fvPatch& patch);

        ~FixedTurbProperties() {}

        const scalarField& getTurbLengthScales() const { return f_tls; }

        const scalarField& getUrms() const { return f_umrs; }

        const scalarField& getDissipationRates() const { return f_eps; }

        const scalarField& getTimeScales() const { return f_tts; }

        void update(const vectorField& refVelocity, const scalar &timeValue);

        void write(Ostream& os) const;

        synTurbulenceParameters* clone(const fvPatch &patch) const;


        void autoMap(const fvPatchFieldMapper&) {}

        //- Reverse map the given fvPatchField onto this fvPatchField
        void rmap(const synTurbulenceParameters&, const labelList&) {}
    };
}

#endif // FIXEDTURBPROPERTIES_H
