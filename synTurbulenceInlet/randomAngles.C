#include "randomAngles.H"

#include "constants.H"
#include "UPstream.H"
#include "OPstream.H"
#include "IPstream.H"

#include <iostream>

namespace Foam
{

    randomAngles::randomAngles(int nmodes):
         m_phi(nmodes, 0.), m_psi(nmodes, 0.), m_alpha(nmodes, 0.), m_tetha(nmodes, 0.), m_random(label(0))
    {
    }

    void randomAngles::RecomputeRandomAngles()
    {
        if(Pstream::master())
        {
            using constant::mathematical::pi;
            scalar ang;

            forAll(m_phi, id)
            {
                m_phi[id]  = scalar01()*2*pi;
                m_psi[id]   = scalar01()*2*pi;
                m_alpha[id] = scalar01()*2*pi;
                ang = scalar01();
                m_tetha[id] = std::acos(1.-ang/0.5);
            }

            for( int slave=Pstream::firstSlave(); slave<=Pstream::lastSlave(); slave++ )
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << m_phi;
                toSlave << m_alpha;
                toSlave << m_psi;
                toSlave << m_tetha;
            }
        }
        else
        {
            IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo() );
            fromMaster >> m_phi;
            fromMaster >> m_alpha;
            fromMaster >> m_psi;
            fromMaster >> m_tetha;
        }
    }
}
