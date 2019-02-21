#include "synTurbulence.H"

#include "dictionary.H"
#include "Ostream.H"
#include "constants.H"
#include "Pstream.H"

namespace Foam
{
    synTurbulence::synTurbulence(const objectRegistry& reg)
    :
        db_(reg),
        m_dxmin(0.5),
        m_nmodes(150),
        m_wew1fct(2.),
        m_visc(1.),
        m_qm(0.015),
        m_sli(0.05),
        m_up(1.),
        m_epsm(1.),
        m_dt(0.1),
        m_T(1.),
        m_u_inf(1.),
        m_ti(0.1),
        angles_(m_nmodes)
    {
        updateParameters();
    }

    synTurbulence::synTurbulence(const objectRegistry& reg, const dictionary& dict)
    :
        db_(reg),
        m_dxmin( 0.5 ),
        m_nmodes( dict.lookupOrDefault<scalar>("nmodes",150) ),
        m_wew1fct(2.),
        m_qm(0.015),
        m_sli(0.05),
        m_up(1.),
        m_epsm(1.),
        m_dt(0.1),
        m_u_inf(1.),
        angles_(m_nmodes)
    {
        dict.lookup("nu") >> m_visc;
        dict.lookup("dxmin") >> m_dxmin;
        dict.lookup("turbIntensity") >> m_ti;
        dict.lookup("turbScale") >> m_sli;

        if(Pstream::master())
            write(Info);

        updateParameters();
    }

    void synTurbulence::updateParameters()
    {
        //Read from transoportProperties dictionary
//        const dictionary& transportProperties = db_.lookupObject<IOdictionary>("transportProperties");
//        dimensionedScalar viscosity(transportProperties.lookup("nu"));
//        m_visc = viscosity.value();
        m_up = m_ti*m_u_inf;
        m_qm = 1.5*m_up*m_up;
        m_epsm = std::pow(m_qm,1.5)/m_sli;

        if(m_u_inf != 0)
        {
            m_T = m_qm/m_epsm;
        }
        else
        {
            m_T = 0;
        }
    }

    void synTurbulence::setRefVelocity(scalar u_inf)
    {
        m_u_inf = u_inf;
        updateParameters();
    }

    void synTurbulence::setTurbulenceIntensity(scalar ti)
    {
        m_ti = ti;
        updateParameters();
    }

    void synTurbulence::setMinDivision(scalar dxmin)
    {
        m_dxmin = dxmin;
        updateParameters();
    }

    void synTurbulence::setNumModes(label nmodes)
    {
        m_nmodes = nmodes;
        updateParameters();
    }

    void synTurbulence::setViscosity(scalar vis)
    {
        m_visc = vis;
        updateParameters();
    }

    void synTurbulence::setCharacteristicLengthScale(scalar sli)
    {
        m_sli = sli;
        updateParameters();
    }

    void synTurbulence::setTimeStep(scalar dt)
    {
        m_dt = dt;
        updateParameters();
    }


    scalar synTurbulence::dxmin() const
    {
        return m_dxmin;
    }

    label synTurbulence::nmodes() const
    {
        return m_nmodes;
    }

    scalar synTurbulence::wew1fct() const
    {
        return m_wew1fct;
    }

    scalar synTurbulence::visc() const
    {
        return m_visc;
    }

    scalar synTurbulence::qm() const
    {
        return m_qm;
    }

    scalar synTurbulence::sli() const
    {
        return m_sli;
    }

    scalar synTurbulence::up() const
    {
        return m_up;
    }

    scalar synTurbulence::epsm() const
    {
        return m_epsm;;
    }

    scalar synTurbulence::dt() const
    {
        return m_dt;
    }

    scalar synTurbulence::T() const
    {
        return m_T;
    }

    void synTurbulence::computeNewFluctuations(const vectorField &coords, vectorField &flucts, bool corelate)
    {
        angles_.RecomputeRandomAngles();

        if(coords.size() == 0)
            return;

        double amp = 1.452762113;
        double pi = constant::mathematical::pi;
        int n = coords.size();
        // Deklaracje zmiennych potrzebnych w programie
        double wnrn, wnre, wnreta, wnr1, dkl;
        double xc, yc, zc, utrp, vtrp, wtrp;
        double kxi, kyi, kzi, sx, sy, sz, kx, ky, kz, rk;
        double arg, tfunk, e, utn;

        // Zaalokuj i wyzeruj tablice
        double *wnr, *wnrf, *dkn, *kxio, *kyio, *kzio;
        double *sxio, *syio, *szio;

        double a, b;

        a = std::exp(-dt()/T());
        b = std::sqrt(1-a*a);

        //temporal tables for fluctuation
        double *u = new double[n];
        double *v = new double[n];
        double *w = new double[n];

        wnr = new double[nmodes()];
        wnrf = new double[nmodes()+1];
        dkn = new double[nmodes()];
        kxio = new double[nmodes()];
        kyio = new double[nmodes()];
        kzio = new double[nmodes()];
        sxio = new double[nmodes()];
        syio = new double[nmodes()];
        szio = new double[nmodes()];


        // Najwyzsza liczba falowa
        wnrn = 2*pi/dxmin();

        // Liczba falowa zwiazana z modem o najwiekszej energii
        wnre = 9*pi*amp/(55*sli());

        // Liczba falowa uzyta w czlonie lepkim w spektrum von Karmana
        wnreta = std::pow((epsm()/std::pow(visc(), 3.0)), 0.25);

        // Najmniejsza liczba falowa
        wnr1 = wnre/wew1fct();

        // Krok miedzy liczbami falowymi
        dkl = (wnrn-wnr1)/nmodes();

        // Wyznacz rownomierny ciag liczb falowych
        for(int m = 0; m<=nmodes(); ++m)
            wnrf[m] = wnr1+dkl*(double)m;

        // Wyznacz liczby falowe w srodkach podzialow
        for(int m = 0; m<nmodes(); ++m)
        {
            wnr[m] = 0.5*(wnrf[m] + wnrf[m+1]);
            dkn[m] = wnrf[m+1] -wnrf[m];
        }

        // Wektory falowe z losowych katow
        for(int m = 0; m<nmodes(); ++m)
        {
            kxio[m] = std::sin(angles_.GetTetha()[m])*std::cos(angles_.GetPhi()[m]);
            kyio[m] = std::sin(angles_.GetTetha()[m])*std::sin(angles_.GetPhi()[m]);
            kzio[m] = std::cos(angles_.GetTetha()[m]);

            // sigma (s = sigma) z losowych katow. Sigma jest wersorem, ktory
            // wyznacza kierunek syntetycznego wektora predkosci (u, v, w)
            sxio[m]=std::cos(angles_.GetPhi()[m])*std::cos(angles_.GetTetha()[m])*std::cos(angles_.GetAlpha()[m])-std::sin(angles_.GetPhi()[m])*std::sin(angles_.GetAlpha()[m]);
            syio[m]=std::sin(angles_.GetPhi()[m])*std::cos(angles_.GetTetha()[m])*std::cos(angles_.GetAlpha()[m])+std::cos(angles_.GetPhi()[m])*std::sin(angles_.GetAlpha()[m]);
            szio[m]=-std::sin(angles_.GetTetha()[m])*std::cos(angles_.GetAlpha()[m]);
        }

        // Petla po siatce
        for(int k = 0; k<n; ++k)
        {
            xc = coords[k].x();
            yc = coords[k].y();
            zc = coords[k].z();

            // Zainicjalizuj fluktuacje turbulentne zerem
            utrp = 0;
            vtrp = 0;
            wtrp = 0;

            // Petla po wszystkich liczbach falowych
            for(int m = 0; m<nmodes(); ++m)
            {
                kxi = kxio[m];
                kyi = kyio[m];
                kzi = kzio[m];

                sx = sxio[m];
                sy = syio[m];
                sz = szio[m];

                kx = kxi*wnr[m];
                ky = kyi*wnr[m];
                kz = kzi*wnr[m];
                rk = std::sqrt(kx*kx + ky*ky + kz*kz);

                // Jesli liczba falowa rk jest mniejsza niz najwieksza dopuszczalna,
                // wygeneruj fluktuacje
                if(rk < wnrn)
                {
                    arg = kx*xc + ky*yc + kz*zc + angles_.GetPsi()[m];

                    tfunk = std::cos(arg);

                    // Spektrum von Karmana
                    e = amp/wnre*std::pow(wnr[m]/wnre,4.0)/(std::pow(1+std::pow(wnr[m]/wnre,2.0),17./6.))*std::exp(-2*std::pow(wnr[m]/wnreta,2.0));

                    utn = std::sqrt(e*up()*up()*dkn[m]);

                    // Syntetyczne pole predkosci
                    utrp += 2*utn*tfunk*sx;
                    vtrp += 2*utn*tfunk*sy;
                    wtrp += 2*utn*tfunk*sz;
                }
            }

            u[k] = utrp;
            v[k] = vtrp;
            w[k] = wtrp;

        } // Koniec petli po siatce

        if(corelate)
        {
            for(int k = 0; k<n; ++k)
            {
                flucts[k].x() = flucts[k].x()*a + b*u[k];
                flucts[k].y() = flucts[k].y()*a + b*v[k];
                flucts[k].z() = flucts[k].z()*a + b*w[k];
            }
        }
        else
        {
            for(int k = 0; k<n; ++k)
            {
                flucts[k].x() = u[k];
                flucts[k].y() = v[k];
                flucts[k].z() = w[k];
            }
        }


        // Usun tablice
        delete [] u;
        delete [] v;
        delete [] w;

        delete [] wnrf;
        delete [] dkn;
        delete [] kxio;
        delete [] kyio;
        delete [] kzio;
        delete [] sxio;
        delete [] syio;
        delete [] szio;
    }

    void synTurbulence::write(Ostream& os) const
    {
        os.writeKeyword("nu") << m_visc << token::END_STATEMENT << nl;
        os.writeKeyword("dxmin") << m_dxmin << token::END_STATEMENT << nl;
        os.writeKeyword("turbIntensity") << m_ti << token::END_STATEMENT << nl;
        os.writeKeyword("turbScale") << m_sli << token::END_STATEMENT << nl;
    }

} //namespace Foam

