#include "GaussianMPS.h"
#include "ReadInput.h"
using namespace itensor;
using namespace std;

CMatrix read_determinants (const string& file, int N)
{
    vector<vector<Real>> data = readtxt<Real> (file);

    int Nlines = data.size();
    assert (Nlines % N == 0);
    int Npar = Nlines / N;

    cout << "particle number = " << Npar << endl;
    if (Npar > N)
    {
        cout << "Error: Lx and/or Ly not correct. Particle number is larger than the number of sites" << endl;
        throw;
    }

    CMatrix phi (N, Npar);
    int it = 0;
    for(int j = 0; j < Npar; j++)
        for(int i = 0; i < N; i++)
        {
            Real real = data[it][0],
                 imag = data[it][1];
            /*if (imag == 0.)
            {
                phi(i,j) = real;
            }
            else
            {*/
                phi(i,j) = complex<Real> (real,imag);
            //}
            it++;
        }
    return phi;
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
    string infile = argv[1];

    //Read in individual parameters from the input file
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    InputGroup input (infile,"basic");

    auto det_up_file   = input.getString("det_up_file");
    auto det_dn_file   = input.getString("det_dn_file");
    auto block_size = input.getInt("block_size");
    auto occ_crit   = input.getReal("occ_crit");
    auto lx         = input.getInt("Lx");
    auto ly         = input.getInt("Ly");
    auto cutoff     = input.getReal("Cutoff");
    auto maxdim     = input.getInt("MaxDim");
    auto outfile    = input.getString("MPSFileName");
    auto grandcanon = input.getYesNo("grand_canonical",false);
    auto conserveSz = input.getYesNo("conserve_Sz",true);

    auto phi_up = read_determinants (det_up_file, lx*ly);
    auto phi_dn = read_determinants (det_dn_file, lx*ly);

    MPS psi;
    while (true)
    {
        try
        {
            cout << "block size = " << block_size << endl;
            psi = GaussianMPS (phi_up, phi_dn, block_size, occ_crit, {"ConserveNf",!grandcanon,"ConserveSz",conserveSz,"Cutoff",cutoff,"MaxDim",maxdim});
            break;
        }
        catch (const std::underflow_error& e)
        {
            block_size++;
        }
    }

    auto sites = Electron (siteInds(psi));
    auto Nop = Make_NMPO (sites);
    cout << "totN = " << inner(psi,Nop,psi) << endl;

    cout << "max bond dim = " << maxLinkDim(psi) << endl;
    writeToFile (outfile, psi);
    return 0;
}
