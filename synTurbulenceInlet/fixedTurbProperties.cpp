#include "fixedTurbProperties.h"

#include "scalar.H"

namespace Foam {

    const word FixedTurbProperties::typeName = "fixedTurbProperites";
    word FixedTurbProperties::type() const { return FixedTurbProperties::typeName; }


    FixedTurbProperties::FixedTurbProperties(const dictionary& dict, const fvPatch& patch) :
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size()), overwriteTimeScale(false)
    {
        dict.lookup("turbIntensity") >> m_ti;
        dict.lookup("turbScale") >> m_tls;

        if(dict.found(TurbProperties::TIME_SCALE_PROP_NAME)) {
            overwriteTimeScale = true;
            f_tts = dict.lookupType<scalar>(TurbProperties::TIME_SCALE_PROP_NAME);
        }
    }

    FixedTurbProperties::FixedTurbProperties(const FixedTurbProperties &other, const fvPatch &patch) :
        m_tls(other.m_tls), m_ti(other.m_ti), overwriteTimeScale(other.overwriteTimeScale),
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size())
    {
        if(overwriteTimeScale && other.f_tts.size() > 0) {
            f_tts = other.f_tts[0];
        }
//        else {
//            FatalErrorInFunction << "Can't copy time scale with overwritten fixed scale, source has 0 size list of time scale values" << exit(FatalError);
//        }
    }

    void FixedTurbProperties::update(const vectorField& refVelocity, const scalar &timeValue) {

        // Pout << "Max refVelocity: " << max(refVelocity) << endl;

        f_umrs = m_ti * max(mag(refVelocity));

        scalarField k = 3./2 * (f_umrs * f_umrs);

        f_eps = pow(k, 3./2) / (m_tls + SMALL); //pow(TurbProperties::Cmu, 0.75) *

        if(!overwriteTimeScale) {
            if( max(f_umrs) > SMALL )
                f_tts = k / (f_eps + SMALL);
            else //perhaps error should be thorwn
                f_tts = 1.0;
        }

        f_tls = m_tls;
    }

    void FixedTurbProperties::write(Ostream& os) const {
        os.beginBlock("properties");
        os.writeEntry("type", "fixed");
        os.writeEntry("turbScale", m_tls);
        os.writeEntry("turbIntensity", m_ti);

        if(overwriteTimeScale && f_tts.size() > 0) {
            os.writeEntry(TurbProperties::TIME_SCALE_PROP_NAME, f_tts[0]);
        }

//        os.writeEntry("turbTimeScale", m_tts);
//        os.writeEntry("turbDissiaptionRate", m_eps);
//        os.writeEntry("turbUp", m_up);

        os.endBlock();
    }

    synTurbulenceParameters *FixedTurbProperties::clone(const fvPatch &patch) const {
        return new FixedTurbProperties(*this, patch);
    }
}
