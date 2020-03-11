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

        void updateParameters();

        // Coefficient required for von Karman spectral energy
        static const scalar Ce;

        // Kolmogorove wavelength
        scalar kEtha(const scalar& nu, const scalar& eps) const {
            return pow(eps / (nu*nu*nu), 0.25);
        }
        inline scalar kEtha() const {
            return kEtha(visc(), m_epsm);
        }

        // highest value of possible wavelengths under given mesh grid
        scalar kMax() const {
            using namespace constant::mathematical;
            return pi / m_dxmin;
        }

        // the most energetic wavelength
        scalar kE() const {
            using namespace constant::mathematical;
            return (Ce * 9 * pi) / (55 * m_sli);
        }

        // smallest possible wavelengths
        scalar kMin() const {
            return kE() / m_minWavelengthFactor;
        }

        // Kolmogorov energy spectrum
        scalar E(const scalar& k, const scalar& u_rms, const scalar& k_e, const scalar& k_eta) const {
            scalar k_ratio = k / k_e;
            return Ce * u_rms*u_rms / k_e * pow(k_ratio, 4) / (1 + k_ratio*k_ratio, 17.0/6) * exp(-2*(k/k_eta));
        }
        inline scalar E(const scalar& k, const scalar& u_rms) const {
            return E(k, u_rms, kE(), kEtha());
        }
        inline scalar E(const scalar& k) const {
            return E(k, up(), kE(), kEtha());
        }

        // velocity ampl. for the fourier mode related to given wavelength k (based on energy contribution from wavelengths in the range k < ... < k + dk)
        scalar uAmpl(const scalar& k, const scalar& dk, const scalar& u_rms, const scalar& k_max, const scalar& k_eta) {
            return sqrt(E(k, u_rms, k_max, k_eta) * dk);
        }
        inline scalar uAmpl(const scalar& k, const scalar& dk, const scalar& u_rms) {
            return uAmpl(k, dk, u_rms, kE(), kEtha());
        }
        inline scalar uAmpl(const scalar& k, const scalar& dk) {
            return uAmpl(k, dk, up(), kE(), kEtha());
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
        label nmodes() const {
            return m_nmodes;
        }

        scalar wew1fct() const {
            return m_minWavelengthFactor;
        }
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

        void computeNewFluctuations(const vectorField & coords, vectorField & flucts, bool corelate=true);

    };

} //end namespace foam

#endif // SYNTURBPARAMETERS_H
