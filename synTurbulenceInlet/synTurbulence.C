#include "synTurbulence.H"

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include <limits>

namespace Foam
{
    const scalar Cmu = 0.09;

    class intensityScaleParameters : public synTurbulenceParameters {
        scalar m_tls;
        scalar m_up;
        scalar m_eps;
        scalar m_tts;
        scalar m_ti;
        tmp<scalarField> f_tls;
        tmp<scalarField> f_umrs;
        tmp<scalarField> f_eps;
        tmp<scalarField> f_tts;

        void computeEpsAndTimeScale() {
            scalar k = 1.5*m_up*m_up;
            m_eps = pow(Cmu, 0.75) * pow(k, 1.5) / m_tls;

            if( m_up > 1e-12 ) {
                m_tts = k / m_eps;
            }
            else {
                m_tts = 1.0;
            }
        }

    public:
        intensityScaleParameters(const dictionary& dict) {

            dict.lookup("turbIntensity") >> m_ti;

            dict.lookup("turbScale") >> m_tls;
        }

        intensityScaleParameters(const scalar& intensity, const scalar& lengthScale):
            m_tls(lengthScale), m_ti(intensity)
        {
        }

        const scalarField& getTurbLengthScales() const {
            return f_tts.ref();
        }

        const scalarField& getUrms() const {
            return f_umrs.ref();
        }

        const scalarField& getDissipationRates() const {
            return f_eps.ref();
        }

        const scalarField& getTimeScales() const {
            return f_tts.ref();
        }

        void update(const vectorField &coords, const vectorField& refVelocity, const scalar &timeValue) {
            if( ! f_tls.valid() ) {
                m_up = m_ti * average(mag(refVelocity));
                computeEpsAndTimeScale();
                f_tls = tmp<scalarField>(new scalarField(coords.size(), m_tls));
                f_umrs = tmp<scalarField>(new scalarField(coords.size(), m_up));
                f_eps = tmp<scalarField>(new scalarField(coords.size(), m_eps));
                f_tts = tmp<scalarField>(new scalarField(coords.size(), m_tts));
            }
        }

        virtual void write(Ostream& os) const {
            os.beginBlock("properties");
            os.writeEntry("type", "fixed");
            os.writeEntry("turbScale", m_tls);
            os.writeEntry("turbTimeScale", m_tts);
            os.writeEntry("turbDissiaptionRate", m_eps);
            os.writeEntry("turbIntensity", m_ti);
            os.writeEntry("turbUp", m_up);
            os.endBlock();
        }
    };

    static scalar smallestEdgeLength(const fvMesh& mesh) {
        const edgeList& edges = mesh.edges();
        const pointField& pp = mesh.points();
        scalar min = std::numeric_limits<scalar>::max();
        forAll (edges, edgei) {
            scalar el = edges[edgei].mag(pp);
            if (el < min) {
                min = el;
            }
        }
        return min;
    }

    const scalar synTurbulence::Ce = 1.452762113;

    synTurbulence::synTurbulence(const objectRegistry& reg)
    :
        db_(reg),
        m_dxmin(0.5),
        m_nmodes(150),
        m_minWavelengthFactor(2.),
        m_visc(1.),
        m_qm(0.015),
        m_sli(0.05),
        m_up(1.),
        m_epsm(1.),
        m_dt(0.1),
        m_T(1.),
        m_u_inf(1.),
        m_ti(0.1),
        m_angles(m_nmodes)
    {
        updateParameters();
    }

    synTurbulence::synTurbulence(const objectRegistry& reg, const dictionary& dict, const fvMesh& mesh)
    :
        db_(reg),
        m_dxmin( 0.5 ),
        m_nmodes( dict.lookupOrDefault<scalar>("nmodes", 150) ),
        m_minWavelengthFactor(2.),
        m_qm(0.015),
        m_sli(0.05),
        m_up(1.),
        m_epsm(1.),
        m_dt(0.1),
        m_u_inf(1.),
        m_angles(m_nmodes)
    {
        dict.lookup("nu") >> m_visc;

        if( dict.found("dxmin") ){
            dict.lookup("dxmin") >> m_dxmin;
        }
        else {
           m_dxmin = smallestEdgeLength(mesh);
        }


        dict.lookup("turbIntensity") >> m_ti;
        dict.lookup("turbScale") >> m_sli;


        const dictionary& props = dict.subDict("properties");
        word propType;
        props.lookup("type") >> propType;


        if(propType == "fixed") {
            properites = tmp<synTurbulenceParameters>(new intensityScaleParameters(props));
        }
        else {
            FatalError << "Properties type not supported";
        }

        //no need, info stream automatically is handled only in master process
        if(Pstream::master())
            write(Info);

        updateParameters();
    }

    void synTurbulence::updateParameters() {
        //Read from transoportProperties dictionary
//        const dictionary& transportProperties = db_.lookupObject<IOdictionary>("transportProperties");
//        dimensionedScalar viscosity(transportProperties.lookup("nu"));
//        m_visc = viscosity.value();
        m_up = m_ti*m_u_inf;
        m_qm = 1.5*m_up*m_up;
        m_epsm = Foam::pow(m_qm, 1.5) / m_sli;

        if(m_u_inf != 0) {
            m_T = m_qm/m_epsm;
        }
        else {
            m_T = 0;
        }
    }

    void synTurbulence::setRefVelocity(scalar u_inf) {
        m_u_inf = u_inf;
        updateParameters();
    }

    void synTurbulence::setTurbulenceIntensity(scalar ti) {
        m_ti = ti;
        updateParameters();
    }

    void synTurbulence::setMinDivision(scalar dxmin) {
        m_dxmin = dxmin;
        updateParameters();
    }

    void synTurbulence::setNumModes(label nmodes) {
        m_nmodes = nmodes;
        updateParameters();
    }

    void synTurbulence::setViscosity(scalar vis) {
        m_visc = vis;
        updateParameters();
    }

    void synTurbulence::setCharacteristicLengthScale(scalar sli) {
        m_sli = sli;
        updateParameters();
    }

    void synTurbulence::setTimeStep(scalar dt) {
        m_dt = dt;
        updateParameters();
    }


    tmp<scalarField> linspace(scalar min, scalar max, label size) {
        tmp<scalarField> ptr = tmp<scalarField>(new scalarField(size, min));
        scalarField& values = ptr.ref();
        scalar dx = (min - max) / (size -1);
        for(label i = 1; i<size; ++i) {
            values[i] += i*dx;
        }
        return ptr;
    }

    void synTurbulence::computeNonuniformFlucts(const vectorField &coords, const vectorField& refVelocity, const scalar& timeValue, vectorField &rFlucts, bool corelate) {
        using constant::mathematical::pi;

        m_angles.RecomputeRandomAngles();
        const scalarField& theta = m_angles.GetTetha();
        const scalarField& phi = m_angles.GetPhi();
        const scalarField& alpha = m_angles.GetAlpha();
        const scalarField& psi = m_angles.GetPsi();

        if(coords.size() == 0)
            return;

        properites->update(coords, refVelocity, timeValue);

        scalar kmin = kMin( max(properites->getTurbLengthScales()) );
        scalar kmax = kMax();

        scalarField kEthaField = kolmogorovWavelength(visc(), properites->getDissipationRates());
        const scalarField& UrmsField = properites->getUrms();

        label N = coords.size();

        vectorField newFlucts(N, vector::zero);


        // compute equaly distributed wave numbers, in the center of discret spacing
        scalar kSpacing = (kmax - kmin) / (nmodes() + 1);
        scalarField wavelengths(nmodes());
        forAll(wavelengths, i) {
            wavelengths[i] = kSpacing/2 + i*kSpacing;
        }

        // Compute wave vectors using random angles
        vectorField wavevectors(nmodes(), vector::zero);
        wavevectors.replace(vector::X, sin(theta)*cos(phi)*wavelengths);
        wavevectors.replace(vector::Y, sin(theta)*sin(phi)*wavelengths );
        wavevectors.replace(vector::Z, cos(theta)*wavelengths );

        // Compute sigma also using random angles. This field will guarantee zero divergence
        // It's a vector that would apply direction to the flucts velocity
        vectorField sigma(nmodes(), vector::zero);
        sigma.replace(vector::X,  cos(phi)*cos(theta)*cos(alpha)-sin(phi)*sin(alpha) );
        sigma.replace(vector::Y,  sin(phi)*cos(theta)*cos(alpha)+cos(phi)*sin(alpha) );
        sigma.replace(vector::Z, -sin(theta)*cos(alpha) );

        //loop over mesh boundary points
        for(int i=0; i < N; ++i) {
            const vector& pntI = coords[i];
            const scalar& urmsI = UrmsField[i];
            const scalar& kEthaI = kEthaField[i];
            vector& fluctsI = newFlucts[i];

            //loop over modes
            for(int m=0; m < nmodes(); ++m) {
                //note: I've removed if condition from original code.
                // The asseration was meaningless and if k spacing would change it might potentially remve highest wavelengths from spectrum
                // Wavevectors are generated as vector pointing to point on the sphere(R=1), thus any wavevector meets condition
                // mag(k) < kmax as k is lengths are generated in range <kmin, kmax> (even (kmin, kmax) )

                //flucts amplitude based on von Karman spectrum
                scalar u = uAmpl(wavelengths[m], kSpacing, urmsI, kmax, kEthaI);

                // fourier series component based on equation:
                // \vec{v_m} = u_m * cos (\vec{k_m} \cdot \vec{x} + psi_m) \cdot \vec{sigma_m}
                // where u_m is m-th amplitude, k is wave vector, x coordinate vector, psi phase shift, sigma direction vector
                fluctsI += (u * cos( (wavevectors[m] & pntI) + psi[m])) * sigma[m];
            }
            fluctsI *= 2;
        }

        if(corelate) {
            //use nonuniform time scales, respectively to given point at boundary
            scalarField a = exp( -dt()/properites->getTimeScales() );
            rFlucts = rFlucts * a + sqrt(1.0 - a*a)*newFlucts;
        }
        else {
            rFlucts = newFlucts;
        }

        Info <<"min/max flucts " << min(mag(rFlucts)) <<"/"<<max(mag(rFlucts)) << endl;
    }



    void synTurbulence::computeNewFluctuations(const vectorField &coords, vectorField &flucts, bool corelate) {
        using namespace constant::mathematical;
        using namespace Foam;

        m_angles.RecomputeRandomAngles();

        if(coords.size() == 0)
            return;

        double amp = 1.452762113;
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

        a = exp(-dt()/T());
        b = sqrt(1-a*a);

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
        wnrn = pi / dxmin();

        // Liczba falowa zwiazana z modem o najwiekszej energii
        wnre = 9*pi*amp/(55*sli());

        // Liczba falowa uzyta w czlonie lepkim w spektrum von Karmana
        wnreta = 2*pi * pow((epsm()/pow(visc(), 3.0)), 0.25);

        // Najmniejsza liczba falowa
        wnr1 = wnre/wew1fct();

        // Krok miedzy liczbami falowymi
        dkl = (wnrn-wnr1)/nmodes();

        // Wyznacz rownomierny ciag liczb falowych
        for(int m = 0; m<=nmodes(); ++m)
            wnrf[m] = wnr1 + dkl*m;

        // Wyznacz liczby falowe w srodkach podzialow
        for(int m = 0; m < nmodes(); ++m) {
            wnr[m] = 0.5*(wnrf[m] + wnrf[m+1]);
            dkn[m] = wnrf[m+1] -wnrf[m];
        }

        // Wektory falowe z losowych katow
        for(int m = 0; m<nmodes(); ++m)
        {
            kxio[m] = sin(m_angles.GetTetha()[m])*cos(m_angles.GetPhi()[m]);
            kyio[m] = sin(m_angles.GetTetha()[m])*sin(m_angles.GetPhi()[m]);
            kzio[m] = cos(m_angles.GetTetha()[m]);

            // sigma (s = sigma) z losowych katow. Sigma jest wersorem, ktory
            // wyznacza kierunek syntetycznego wektora predkosci (u, v, w)
            sxio[m]=cos(m_angles.GetPhi()[m])*cos(m_angles.GetTetha()[m])*cos(m_angles.GetAlpha()[m])-sin(m_angles.GetPhi()[m])*sin(m_angles.GetAlpha()[m]);
            syio[m]=sin(m_angles.GetPhi()[m])*cos(m_angles.GetTetha()[m])*cos(m_angles.GetAlpha()[m])+cos(m_angles.GetPhi()[m])*sin(m_angles.GetAlpha()[m]);
            szio[m]=-sin(m_angles.GetTetha()[m])*cos(m_angles.GetAlpha()[m]);
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
                rk = sqrt(kx*kx + ky*ky + kz*kz);

                // Jesli liczba falowa rk jest mniejsza niz najwieksza dopuszczalna,
                // wygeneruj fluktuacje
                if(rk < wnrn)
                {
                    arg = kx*xc + ky*yc + kz*zc + m_angles.GetPsi()[m];

                    tfunk = cos(arg);

                    // Spektrum von Karmana
                    e = amp/wnre*pow(wnr[m]/wnre,4.0)/(pow(1+pow(wnr[m]/wnre,2.0),17./6.))*exp(-2*pow(wnr[m]/wnreta,2.0));

                    utn = sqrt(e*up()*up()*dkn[m]);

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

        properites->write(os);
    }

} //namespace Foam

