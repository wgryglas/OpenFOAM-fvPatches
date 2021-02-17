#include "interpolatedTurbProperties.h"

#include "dictionary.H"
#include "scalar.H"

/*
 * epsilon =  0.09*k*omega (dyssypacja)
 * lt = k^3/2 / epsilon (skala dlugosci turbulentnej)
 *
 * k = 3/2 up^2 ---> up = sqrt(2/3*k)
*/

namespace Foam {
    const word InterpolatedTurbProperties::typeName = "interpolatedTurbProperites";
    word InterpolatedTurbProperties::type() const { return InterpolatedTurbProperties::typeName; }

    const scalar Cmu = 0.09;

    InterpolatedTurbProperties::InterpolatedTurbProperties(const InterpolatedTurbProperties &other, const fvPatch &patch):
        kMapper( new PatchFunction1Types::MappedFile<scalar>(other.kMapper(), patch.patch()) ),
        omegaMapper( new PatchFunction1Types::MappedFile<scalar>(other.omegaMapper(), patch.patch()) ),
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size()),
        m_nu(other.m_nu)
    {
    }

    InterpolatedTurbProperties::InterpolatedTurbProperties(const dictionary& dict, scalar nu, const fvPatch &patch):
        kMapper(
            new PatchFunction1Types::MappedFile<scalar>(
                patch.patch(),
                "kMapper",
                dict.subDict("k"),
                "k",
                true
            )
        ),
        omegaMapper(
            new PatchFunction1Types::MappedFile<scalar>(
                patch.patch(),
                "omegaMapper",
                dict.subDict("omega"),
                "omega",
                true
            )
        ),
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size()),
        m_nu(nu)
    {
    }

    void InterpolatedTurbProperties::update(const vectorField& refVelocity, const scalar& timeValue) {
        Info << "Starting update with size: " << refVelocity.size() << endl;

//        scalarField k = kMapper->value(timeValue);
        scalarField k(refVelocity.size(), 6.88E-03);
//        Info << "k range: " << min(k) <<" / " << max(k) << nl;

        Info << "After k update" << endl;
        //scalarField omega = omegaMapper->value(timeValue);
        scalarField omega(refVelocity.size(), 100.0);
//        Info << "omeaga range: " << min(omega) <<" / " << max(omega) << nl;

        Info << "After k/omega update" << endl;

        f_eps = Cmu * k * omega;

//        Info <<"k/omage range: " << min(k / f_eps) <<" / " <<max(k / f_eps) << nl;

        Info << "After eps update" << endl;

        f_umrs = sqrt( (2.0/3) * k);

        Info << "After umrs update" << endl;

        scalarField eps_plus = f_eps + SMALL;

        f_tls = pow(k, 3.0/2) / (eps_plus);

        Info << "After ls update" << endl;

        if( max(f_umrs) > SMALL )
            f_tts = k / (f_eps + SMALL);
        else //perhaps error should be thorwn
            f_tts = 1.0;

        //f_tts = sqrt(m_nu / eps_plus );

        Info << "After ts update" << endl;
    }

    synTurbulenceParameters *InterpolatedTurbProperties::clone(const fvPatch &patch) const {
        return new InterpolatedTurbProperties(*this, patch);
    }

    void InterpolatedTurbProperties::autoMap(const fvPatchFieldMapper &mapper) {
        Info << "Calling InterpolatedTurbProperties::autoMap" << nl;

        kMapper->autoMap(mapper);
        omegaMapper->autoMap(mapper);

    }


    void InterpolatedTurbProperties::rmap(const synTurbulenceParameters &pf1, const labelList & addr) {
        Info << "Calling InterpolatedTurbProperties::rmap" << nl;

        const InterpolatedTurbProperties& tiptf = refCast<const InterpolatedTurbProperties>(pf1);
        kMapper->rmap(tiptf.kMapper(), addr);
        omegaMapper->rmap(tiptf.omegaMapper(), addr);
    }

    void InterpolatedTurbProperties::write(Ostream& os) const {
        os.beginBlock("properties");
            os.writeEntry("type", "interpolated");
            os.beginBlock("k");
                kMapper().writeData(os);
            os.endBlock();
            os.beginBlock("omega");
                omegaMapper().writeData(os);
            os.endBlock();
        os.endBlock();
    }

}

