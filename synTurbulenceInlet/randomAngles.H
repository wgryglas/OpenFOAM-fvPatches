#ifndef RANDOMANGLES_H
#define RANDOMANGLES_H

#include "Random.H"
#include "scalarField.H"

namespace Foam {

    /// \brief class holding random data required to generate fluctuation. Data shall be
    /// computed only at root process and send to other processes to keep the same random
    /// numbers. Updating random numbers shall be performed always by the same process, because
    /// method "RecomputeRandomAngles" initilize random number generator and this shall be done
    /// only once for whole computation process
    class randomAngles
    {
        scalarField m_phi;
        scalarField m_psi;
        scalarField m_alpha;
        scalarField m_tetha;

        Random m_random;

    public:
        randomAngles(int nmodes);
        //Default copy/assign constructors

        virtual ~randomAngles(){}

        inline const scalarField& GetPhi() const;
        inline const scalarField& GetPsi() const;
        inline const scalarField& GetAlpha() const;
        inline const scalarField& GetTetha() const;

//        inline scalarField& UpdatePhi();
//        inline scalarField& UpdatePsi();
//        inline scalarField& UpdateAlpha();
//        inline scalarField& UpdateTetha();

        void RecomputeRandomAngles();
    };

    const scalarField& randomAngles::GetPhi() const
    {
        return m_phi;
    }

    const scalarField& randomAngles::GetPsi() const
    {
        return m_psi;
    }

    const scalarField& randomAngles::GetAlpha() const
    {
        return m_alpha;
    }

    const scalarField& randomAngles::GetTetha() const
    {
        return m_tetha;
    }

//    scalarField& randomAngles::UpdatePhi()
//    {
//        return m_phi;
//    }

//    scalarField& randomAngles::UpdatePsi()
//    {
//        return m_psi;
//    }

//    scalarField& randomAngles::UpdateAlpha()
//    {
//        return m_alpha;
//    }

//    scalarField& randomAngles::UpdateTetha()
//    {
//        return m_tetha;
//    }
}

    #endif // RANDOMANGLES_H
