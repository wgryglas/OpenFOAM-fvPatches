#include "randomAngles.H"

#include "constants.H"
#include "UPstream.H"
#include "OPstream.H"
#include "IPstream.H"

#include <iostream>

namespace Foam
{
    randomAngles::randomAngles(int nmodes):
         m_phi(nmodes, 0.), m_psi(nmodes, 0.), m_alpha(nmodes, 0.), m_tetha(nmodes, 0.), m_random(label(100))
    {
    }

    void randomAngles::RecomputeRandomAngles() {
        using constant::mathematical::pi;
        // Does the communication is necessary here? Since we use the same seed for all the processes than it
        // should produce the same result as long as the number of modes is the same, what will be the case for
        // particular boundary.
//        if(Pstream::master()) {

            forAll(m_phi, id) {
                m_phi[id]  = randomScalar01()*2*pi;
                m_psi[id]   = randomScalar01()*2*pi;
                m_alpha[id] = randomScalar01()*2*pi;
                m_tetha[id] = std::acos(1.0 - 2*randomScalar01());
            }

//            for(label slave=Pstream::firstSlave(); slave<=Pstream::lastSlave(); slave++ ) {
//                OPstream toSlave(Pstream::commsTypes::blocking, slave);
//                toSlave << m_phi;
//                toSlave << m_alpha;
//                toSlave << m_psi;
//                toSlave << m_tetha;
//            }
//        }
//        else {
//            IPstream fromMaster(Pstream::commsTypes::blocking, Pstream::masterNo() );
//            fromMaster >> m_phi;
//            fromMaster >> m_alpha;
//            fromMaster >> m_psi;
//            fromMaster >> m_tetha;
//        }
    }
}
