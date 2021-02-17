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
#include "turbProperties.h"

namespace Foam
{
    class synTurbulence
    {
        const objectRegistry& db_;

        scalar m_dxmin;
        label m_nmodes;
        scalar m_minWavelengthFactor;
        scalar m_visc;

        bool m_stats;


        scalar m_qm;   //Deprecated, used in old implementation, should be removed. Current values are just defualt
        scalar m_sli;  //Deprecated, used in old implementation, should be removed. Current values are just defualt
        scalar m_up;   //Deprecated, used in old implementation, should be removed. Current values are just defualt
        scalar m_epsm; //Deprecated, used in old implementation, should be removed. Current values are just defualt
        scalar m_dt;   //Deprecated, used in old implementation, should be removed. Current values are just defualt
        scalar m_T;    //Deprecated, used in old implementation, should be removed. Current values are just defualt
        scalar m_u_inf;//Deprecated, used in old implementation, should be removed. Current values are just defualt
        scalar m_ti;   //Deprecated, used in old implementation, should be removed. Current values are just defualt

        const fvPatch& m_patch;

        randomAngles m_angles;

        autoPtr<synTurbulenceParameters> properites;

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
        synTurbulence(const objectRegistry& reg, const fvPatch& patch);
        synTurbulence(const objectRegistry& reg, const dictionary& dict, const fvMesh& mesh, const fvPatch& patch);

        synTurbulence(const synTurbulence& other);

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



        bool isPrintingStats() const { return m_stats; }

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

        void computeNonuniformFlucts(const vectorField& refVelocity, const scalar &timeValue, vectorField & flucts, bool corelate=true);

        void computeNewFluctuations(vectorField & flucts, bool corelate=true);


        void autoMap(const fvPatchFieldMapper&);

        void rmap(const synTurbulence& other, const labelList& addrs);
    };

} //end namespace foam

#endif // SYNTURBPARAMETERS_H
