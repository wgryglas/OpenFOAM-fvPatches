#ifndef INTERPOLATEDTURBPROPERTIES_H
#define INTERPOLATEDTURBPROPERTIES_H

#include "synTurbulence.H"


namespace Foam {
    class InterpolatedTurbProperties : public synTurbulenceParameters
    {
        PatchFunction1Types::MappedFile<scalar> kMapper;
        PatchFunction1Types::MappedFile<scalar> omegaMapper;

        scalarField f_tls;
        scalarField f_umrs;
        scalarField f_eps;
        scalarField f_tts;

        scalar m_nu;

        bool overwriteTimeScale;
        autoPtr<PatchFunction1Types::MappedFile<scalar>> timeScaleMapper;


        InterpolatedTurbProperties(const InterpolatedTurbProperties& other, const fvPatch& patch);

        //prevent copy constructors and defualt
        InterpolatedTurbProperties(const InterpolatedTurbProperties& other) = delete;
        InterpolatedTurbProperties() = delete;

    public:
        static const word typeName;
        word type() const;

        InterpolatedTurbProperties(const dictionary& dict, scalar nu, const fvPatch& patch);

        ~InterpolatedTurbProperties() {}

        const scalarField& getTurbLengthScales() const {
            return f_tls;
        }

        const scalarField& getUrms() const {
            return f_umrs;
        }

        const scalarField& getDissipationRates() const {
            return f_eps;
        }

        const scalarField& getTimeScales() const {
            return f_tts;
        }


        void write(Ostream& os) const;

        void update(const vectorField& refVelocity, const scalar& timeValue);

        synTurbulenceParameters* clone(const fvPatch &patch) const;

        void autoMap(const fvPatchFieldMapper&);

        void rmap(const synTurbulenceParameters&, const labelList&);

    };
}


#endif // INTERPOLATEDTURBPROPERTIES_H
