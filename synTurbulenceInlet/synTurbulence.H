#ifndef SYNTURBPARAMETERS_H
#define SYNTURBPARAMETERS_H

#include "scalar.H"
#include "label.H"
#include "objectRegistry.H"
#include "randomAngles.H"
#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "vectorField.H"


namespace Foam
{
    class synTurbulence
    {
        const objectRegistry& db_;

        scalar m_dxmin;
        label m_nmodes;
        scalar m_wew1fct;
        scalar m_visc;
        scalar m_qm;
        scalar m_sli;
        scalar m_up;
        scalar m_epsm;
        scalar m_dt;
        scalar m_T;
        scalar m_u_inf;
        scalar m_ti;

        randomAngles angles_;

    private:
        void updateParameters();

    public:
        synTurbulence(const objectRegistry& reg);
        synTurbulence(const objectRegistry& reg, const dictionary& dict);

        //Default copy/assign constructors

        virtual ~synTurbulence(){}

        void write(Ostream& os) const;

        void setRefVelocity(scalar u_inf);
        void setTurbulenceIntensity(scalar ti);
        void setMinDivision(scalar dxmin);
        void setNumModes(label nmodes);
        void setViscosity(scalar vis);
        void setCharacteristicLengthScale(scalar sli);
        void setTimeStep(scalar dt);


        scalar dxmin() const;
        label nmodes() const;
        scalar wew1fct() const;
        scalar visc() const;
        scalar qm() const;
        scalar sli() const;
        scalar up() const;
        scalar epsm() const;
        scalar dt() const;
        scalar T() const;

        void computeNewFluctuations(const vectorField & coords, vectorField & flucts, bool corelate=true);

    };

} //end namespace foam

#endif // SYNTURBPARAMETERS_H
