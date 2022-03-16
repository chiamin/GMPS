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

int main(int argc, char* argv[])
{
    string infile = argv[1];

    //Read in individual parameters from the input file
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    InputGroup input (infile,"basic");

    auto det_file   = input.getString("det_file");
    auto block_size = input.getInt("block_size");
    auto occ_crit   = input.getReal("occ_crit");
    auto lx         = input.getInt("Lx");
    auto ly         = input.getInt("Ly");
    auto cutoff     = input.getReal("Cutoff");
    auto maxdim     = input.getInt("MaxDim");
    auto outfile    = input.getString("MPSFileName");

    auto phi = read_determinants (det_file, lx*ly);

    auto psi = GaussianMPS (phi, phi, block_size, occ_crit, {"Cutoff",cutoff,"MaxDim",maxdim});
    cout << "max bond dim = " << maxLinkDim(psi) << endl;
    writeToFile (outfile, psi);
    return 0;
}
