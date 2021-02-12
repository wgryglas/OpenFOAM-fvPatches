#ifndef SYNTURBPARAMETERS_H
#define SYNTURBPARAMETERS_H

#include "scalar.H"
#include "label.H"
#include "objectRegistry.H"
#include "randomAngles.H"
#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "vectorField.H"
#include "fvMesh.H"
#include "constants.H"

namespace Foam
{
    struct synTurbulenceParameters: tmp<synTurbulenceParameters>::refCount {                     
        virtual void update(const vectorField& coords, const vectorField& refVelocity, const scalar& timeValue) = 0;
        virtual const scalarField& getTurbLengthScales() const  = 0;
        virtual const scalarField& getUrms() const  = 0;
        virtual const scalarField& getDissipationRates() const  = 0;
        virtual const scalarField& getTimeScales() const = 0;
        virtual ~synTurbulenceParameters() {}
        virtual void write(Ostream& os) const = 0;
    };

    class synTurbulence
    {
        const objectRegistry& db_;

        scalar m_dxmin;
        label m_nmodes;
        scalar m_minWavelengthFactor;
        scalar m_visc;
        scalar m_qm;
        scalar m_sli;
        scalar m_up;
        scalar m_epsm;
        scalar m_dt;
        scalar m_T;
        scalar m_u_inf;
        scalar m_ti;

        randomAngles m_angles;

        tmp<synTurbulenceParameters> properites;

        void updateParameters();

        // Coefficient required for von Karman spectral energy
        static const scalar Ce;

        // Kolmogorove wavelength
        scalar kolmogorovWavelength(const scalar& nu, const scalar& eps) const {
            using constant::mathematical::pi;
            return 2*pi * pow(eps / (nu*nu*nu), 0.25);
        }
        tmp<scalarField> kolmogorovWavelength(const scalar &nu, const scalarField &eps) const {
            tmp<scalarField> tmpRes(new scalarField(eps.size(), 0.0));
            scalarField& res = tmpRes.ref();
            forAll(eps, i) {
                res[i] = kolmogorovWavelength(nu, eps[i]);
            }
            return tmpRes;
        }

        // highest value of possible wavelengths under current mesh grid
        scalar kMax() const {
            using namespace constant::mathematical;
            return pi / m_dxmin;
        }

        // the most energetic wavelength
        scalar kE(const scalar& turbLengthScale) const {
            using namespace constant::mathematical;
            return (Ce * 9 * pi) / (55 * turbLengthScale);
        }

        // smallest possible wavelengths
        scalar kMin(const scalar& turbLengthScale) const {
            return kE(turbLengthScale) / m_minWavelengthFactor;
        }

        // Kolmogorov energy spectrum releated to given wavelength k.
        // u_mrs - mean root square of fluctuations, k_e - the most energetic wavelength, k_etha - Kolmogorow wavelength
        scalar E(const scalar& k, const scalar& u_mrs, const scalar& k_e, const scalar& k_etha) const {
            scalar k_ratio = k / k_e;
            scalar k_etha_ratio = k / k_etha;
            return Ce * (u_mrs*u_mrs / k_e) * pow(k_ratio, 4.0) / pow(1 + k_ratio*k_ratio, 17.0/6) * exp(-2*k_etha_ratio*k_etha_ratio);
        }

        // velocity fluctuation ampl. for the fourier mode related to given wavelength k (based on energy contribution from wavelengths in the range k < ... < k + dk)
        // k - given wavelength, dk - wevelengths valus width used for energy range,
        // u_mrs - mrs of total vel. fluctuation, k_max - maximal possible wavelength and k_etha - kolmogorow wavelength
        scalar uAmpl(const scalar& k, const scalar& dk, const scalar& u_mrs, const scalar& k_max, const scalar& k_eta) {
            return sqrt( E(k, u_mrs, k_max, k_eta) * dk );
        }



    public:
        synTurbulence(const objectRegistry& reg);
        synTurbulence(const objectRegistry& reg, const dictionary& dict, const fvMesh& mesh);

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


        /**
         * @return smallest mesh size
         */
        scalar dxmin() const {
            return m_dxmin;
        }


        /**
         * @return number of modes used in fourier series to represent fluctuations
         */
        inline label nmodes() const {
            return m_nmodes;
        }

        scalar wew1fct() const {
            return m_minWavelengthFactor;
        }

        /**
         * @brief visc - returns molecular viscosity value
         * @return
         */
        scalar visc() const {
            return m_visc;
        }

        scalar qm() const {
            return m_qm;
        }

        /**
         * @return turbulent length scale
         */
        scalar sli() const {
            return m_sli;
        }

        scalar up() const {
            return m_up;
        }
        scalar epsm() const {
            return m_epsm;
        }
        scalar dt() const {
            return m_dt;
        }
        scalar T() const {
            return m_T;
        }

        void computeNonuniformFlucts(const vectorField & coords, const vectorField& refVelocity, const scalar &timeValue, vectorField & flucts, bool corelate=true);

        void computeNewFluctuations(const vectorField & coords, vectorField & flucts, bool corelate=true);

    };

} //end namespace foam

#endif // SYNTURBPARAMETERS_H
