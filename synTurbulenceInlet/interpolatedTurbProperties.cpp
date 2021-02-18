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

    InterpolatedTurbProperties::InterpolatedTurbProperties(const InterpolatedTurbProperties &other, const fvPatch &patch):
        kMapper(other.kMapper, patch.patch()),
        omegaMapper(other.omegaMapper, patch.patch()),
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size()),
        m_nu(other.m_nu),
        overwriteTimeScale(other.overwriteTimeScale)
    {
        Pout << "Cloning " << patch.name() << " with size " << patch.size() << endl;

        if(overwriteTimeScale) {
            if(!other.timeScaleMapper.empty()) {
                timeScaleMapper.set(cloneMapper(patch, other.timeScaleMapper()));
            }
            else if(other.f_tts.size() > 0) {
                f_tts = other.f_tts[0];
            }
        }

//        else {
//            FatalErrorInFunction << "Can't copy time scale with overwritten fixed scale, source has 0 size list of time scale values" << exit(FatalError);
//        }
    }

    InterpolatedTurbProperties::InterpolatedTurbProperties(const dictionary& dict, scalar nu, const fvPatch &patch):
        kMapper(            
            patch.patch(),
            "kMapper",
            dict.subDict("k"),
            "k",
            true
        ),
        omegaMapper(
            patch.patch(),
            "omegaMapper",
            dict.subDict("omega"),
            "omega",
            true
        ),
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size()),
        m_nu(nu),
        overwriteTimeScale(false)
    {
        if(dict.found(TurbProperties::TIME_SCALE_PROP_NAME)) {
            overwriteTimeScale = true;
            if(dict.subDictPtr(TurbProperties::TIME_SCALE_PROP_NAME) != NULL) {
                timeScaleMapper.set(newMapper<scalar>(TurbProperties::TIME_SCALE_PROP_NAME, patch, dict.subDict(TurbProperties::TIME_SCALE_PROP_NAME)));
            }
            else {
                f_tts = dict.lookupType<scalar>(TurbProperties::TIME_SCALE_PROP_NAME);
            }
        }
//        Warning: MappedFile utility requires to run with 0 value (if only 0 value is present as boundaryData) for the first time.
//        Otherwise in parallel execution it crashes.
        kMapper.value(0);
        omegaMapper.value(0);
        if(!timeScaleMapper.empty()) {
            timeScaleMapper().value(0);
        }
    }

    void InterpolatedTurbProperties::update(const vectorField& refVelocity, const scalar& timeValue) {
//        Pout << "Starting update with size: " << refVelocity.size() << endl;

        if(refVelocity.size() == 0) return;

//        Pout << "Time " << timeValue << endl;

        scalarField k = kMapper.value(timeValue);
//        scalarField k(refVelocity.size(), 6.88E-03);

//        Pout << "k range: " << min(k) <<" / " << max(k) << nl;

//        Pout << "After k update" << endl;

        scalarField omega = omegaMapper.value(timeValue);
//        scalarField omega(refVelocity.size(), 100.0);
//        Pout << "omeaga range: " << min(omega) <<" / " << max(omega) << nl;

//        Pout << "After k/omega update" << endl;

        f_eps = TurbProperties::Cmu * k * omega;

//        Info <<"k/omage range: " << min(k / f_eps) <<" / " <<max(k / f_eps) << nl;

//        Pout << "After eps update" << endl;

        f_umrs = sqrt( (2.0/3) * k);

//        Pout << "After umrs update" << endl;

        scalarField eps_plus = f_eps + SMALL;

        f_tls = pow(k, 3.0/2) / (eps_plus);

//        Pout << "After ls update" << endl;

        if(!overwriteTimeScale) {
            Pout << "Using uniform" << endl;
            if( max(f_umrs) > SMALL )
                f_tts = k / (f_eps + SMALL);
            else //perhaps error should be thorwn
                f_tts = 1.0;
            //f_tts = sqrt(m_nu / eps_plus );
        }
        else if(timeScaleMapper.valid()) {
            Pout << "Interpolating time scale" << endl;
            f_tts = timeScaleMapper->value(timeValue);
        }
//        Pout << "After ts update" << endl;
    }

    synTurbulenceParameters *InterpolatedTurbProperties::clone(const fvPatch &patch) const {
        return new InterpolatedTurbProperties(*this, patch);
    }

    void InterpolatedTurbProperties::autoMap(const fvPatchFieldMapper &mapper) {
        Info << "Calling InterpolatedTurbProperties::autoMap" << nl;

        kMapper.autoMap(mapper);
        omegaMapper.autoMap(mapper);
        if(!timeScaleMapper.empty()) {
            timeScaleMapper->autoMap(mapper);
        }
    }


    void InterpolatedTurbProperties::rmap(const synTurbulenceParameters &pf1, const labelList & addr) {
        Info << "Calling InterpolatedTurbProperties::rmap" << nl;

        const InterpolatedTurbProperties& tiptf = refCast<const InterpolatedTurbProperties>(pf1);
        kMapper.rmap(tiptf.kMapper, addr);
        omegaMapper.rmap(tiptf.omegaMapper, addr);
        if(!timeScaleMapper.empty()) {
            timeScaleMapper->rmap(tiptf.timeScaleMapper(), addr);
        }
    }

    void InterpolatedTurbProperties::write(Ostream& os) const {
        os.beginBlock("properties");
            os.writeEntry("type", "interpolated");
            os.beginBlock("k");
                kMapper.writeData(os);
            os.endBlock();
            os.beginBlock("omega");
                omegaMapper.writeData(os);
            os.endBlock();
            if(overwriteTimeScale) {
                if(timeScaleMapper.empty()) {
                    if(f_tts.size() > 0) {
                        os.writeEntry(TurbProperties::TIME_SCALE_PROP_NAME, f_tts[0]);
                    }
                }
                else {
                    os.beginBlock(TurbProperties::TIME_SCALE_PROP_NAME);
                        timeScaleMapper->writeData(os);
                    os.endBlock();
                }
            }
        os.endBlock();
    }

}

