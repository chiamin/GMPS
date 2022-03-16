#include "GaussianMPS.h"
using namespace itensor;
using namespace std;

CMatrix tight_binding_Hamilt (int L, Real t, Real mu, Real damp_fac=1., bool damp_from_right=true, bool verbose=false)
{
    cout << "L = " << L << endl;
    CMatrix H (L,L);
    for(int i = 0; i < L; i++)
    {
        H(i,i) = -mu;
        if (i != L-1)
        {
            int damp_dist = (damp_from_right ? L-2-i : i);
            Real ti = t * pow (damp_fac, damp_dist);
            H(i,i+1) = -ti;
            H(i+1,i) = -ti;
            if (verbose)
                cout << "Hk, t " << i << " = " << ti << endl;
        }
    }
    return H;
}

AutoMPO tight_binding_ampo (const Electron& sites, Real t, Real mu)
{
    int L = length(sites);
    AutoMPO ampo (sites);
    for(int i = 1; i <= L; i++)
    {
        ampo += -mu,"Nup",i;
        ampo += -mu,"Ndn",i;
        if (i != L)
        {
            ampo += -t,"Cdagup",i,"Cup",i+1;
            ampo += +t,"Cup",i,"Cdagup",i+1;
            ampo += -t,"Cdagdn",i,"Cdn",i+1;
            ampo += +t,"Cdn",i,"Cdagdn",i+1;
        }
    }
    return ampo;
}

MPO Make_NMPO (const SiteSet& sites)
{
    AutoMPO ampo (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        ampo += 1.0,"Ntot",i;
    }
    return toMPO (ampo);
}

int main(int argc, char* argv[])
{
    int L = 4;

    Electron sites (L);
    auto ampo = tight_binding_ampo (sites, 1., 0.);
    auto H0 = toMPO (ampo);
    InitState init (sites);
    init.set(1,"UpDn");
    init.set(2,"UpDn");
    auto psi0 = MPS (init);
    auto sweeps = Sweeps(3);
    sweeps.maxdim() = 10,20,40;
    sweeps.cutoff() = 1E-12;
    auto en0 = dmrg (psi0, H0, sweeps, {"Quiet",true});
    cout << "E0 = " << en0 << endl;

    auto H = tight_binding_Hamilt (L, 1., 0.);
    CMatrix _Uik;
    Vector _ens;    
    diagHermitian (H, _Uik, _ens);

    cout << _ens << endl;
    auto U = CMatrix (columns(_Uik,L-2,L));
    cout << (transpose(U)*H*U) << endl;
cout << "U" << U << endl;
    auto psi = GaussianMPS (U, U, 4, 0.01);
    psi.replaceSiteInds (inds(sites));
    cout << "E " << inner(psi,H0,psi) << endl;

    auto Nmpo = Make_NMPO (sites);
    cout << "Npar = " << inner(psi,Nmpo,psi) << endl;
    return 0;
}
