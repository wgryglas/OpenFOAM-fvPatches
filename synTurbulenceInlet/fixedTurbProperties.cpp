#include "fixedTurbProperties.h"

#include "scalar.H"

namespace Foam {
    const scalar Cmu = 0.09;

    const word FixedTurbProperties::typeName = "fixedTurbProperites";
    word FixedTurbProperties::type() const { return FixedTurbProperties::typeName; }


    FixedTurbProperties::FixedTurbProperties(const dictionary& dict, const fvPatch& patch) :
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size())
    {
        dict.lookup("turbIntensity") >> m_ti;
        dict.lookup("turbScale") >> m_tls;
    }

    FixedTurbProperties::FixedTurbProperties(const scalar& intensity, const scalar& lengthScale, const fvPatch &patch) :
        m_tls(lengthScale), m_ti(intensity),
        f_tls(patch.size()), f_umrs(patch.size()), f_eps(patch.size()), f_tts(patch.size())
    {
    }

    void FixedTurbProperties::update(const vectorField& refVelocity, const scalar &timeValue) {
        f_umrs = m_ti * max(mag(refVelocity));

        scalarField k = 3./2 * (f_umrs * f_umrs);

        f_eps = pow(k, 3./2) / (m_tls + SMALL); //pow(Cmu, 0.75) *

        if( max(f_umrs) > SMALL )
            f_tts = k / (f_eps + SMALL);
        else //perhaps error should be thorwn
            f_tts = 1.0;

        f_tls = m_tls;
    }

    void FixedTurbProperties::write(Ostream& os) const {
        os.beginBlock("properties");
        os.writeEntry("type", "fixed");
        os.writeEntry("turbScale", m_tls);
        os.writeEntry("turbIntensity", m_ti);

//        os.writeEntry("turbTimeScale", m_tts);
//        os.writeEntry("turbDissiaptionRate", m_eps);
//        os.writeEntry("turbUp", m_up);

        os.endBlock();
    }

    synTurbulenceParameters *FixedTurbProperties::clone(const fvPatch &patch) const {
        return new FixedTurbProperties(m_ti, m_tls, patch);
    }
}
