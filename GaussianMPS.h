#ifndef __GAUSSIANMPS_H_CMC__
#define __GAUSSIANMPS_H_CMC__
#include "itensor/all.h"
#include "IUtility.h"
#include "GeneralUtility.h"
using namespace itensor;
using namespace std;

using VecT = CVector;
using MatT = CMatrix;

CVector inline
operator*(CMatrixRefc const& A,
          CVectorRefc const& b)
{
    CVector res(nrows(A));
    mult(A,b,makeRef(res));
    return res;
}

template <typename MatT>
void check_unitary (const MatT& M)
{
    auto v1 = column(M,0);
    auto v2 = column(M,1);
    mycheck (abs(norm(v1)-1) < 1e-12, "norm != 1");
    mycheck (abs(norm(v2)-1) < 1e-12, "norm != 1");
    mycheck (abs(v1*v2) < 1e-12, "not ortho");
}

// Get a series of 2x2 matrices R, such that
// [     ] [       ]    [ x ]
// |  R  | | orbit |  = |   |
// [     ] [       ]    [ 0 ]
//
// Return: vector = [R_N-1, R_N-2, ..., R_1]
// R_i applies on the i-th and (i+1)-th elements of <orbit>
//
// If one applies R_N-1, and then R_N-2, ..., and then R_1 to <orbit>,
// the elements of <orbit> will be rotated to zero from bottom to top
// and finally only the first element is 1, and the others are all zero.
template <typename VecType>
vector<auto> Get_rot_mats (VecType orbit_)
{
    auto tmp = orbit_(0);
    typedef typename std::conditional<is_same_v<decltype(tmp),double>, Matrix, CMatrix>::type MatT;
    typedef typename std::conditional<is_same_v<decltype(tmp),double>, Vector, CVector>::type VecT;
    auto orbit = VecT(orbit_);

    // Define the matrix type to Matrix or CMatrix depends on weather <orbit> is a real or complex vector
    //typedef typename std::conditional<is_same_v<decltype(orbit(0)),Cplx>, CMatrix, Matrix>::type MatT;
    vector<MatT> rots;
    for(int j = orbit.size()-2; j >= 0; j--)
    {
        MatT rot (2,2);
        rot(0,0) = iut::conj(orbit(j));
        rot(0,1) = iut::conj(orbit(j+1));
        rot(1,0) = -orbit(j+1);
        rot(1,1) = orbit(j);
        auto sub = subVector (orbit,j,j+2);
        rot *= (1./norm(sub));
        subVector (orbit,j,j+2) &= rot * sub;
        rots.push_back (rot);
    }

    return rots;
}

int get_approx_occ (Real occ, Real crit)
{
    Real occ_crit = 1.-occ;
    Real emp_crit = occ;
    if (occ_crit > crit and emp_crit > crit)
    {
        cout << "Cannot find an occupied or empty orbital" << endl;
        cout << "Occupation number = " << occ << endl;
        cout << "Critiria = " << crit << endl;
        throw;
    }
    if (occ_crit < emp_crit)   // occupied
        return 1;
    else   // empty
        return 0;
}

template <typename MatT>
tuple<int,auto> get_approx_occ (const Vector& occ, const MatT& U, Real crit)
{
    int N = occ.size();
    Real occ_crit = 1.-occ(0);
    Real emp_crit = occ(N-1);
    if (occ_crit > crit and emp_crit > crit)
    {
        cout << "Cannot find an occupied or empty orbital" << endl;
        cout << "Occupation number = " << occ(0) << "  " << occ(N-1) << endl;
        cout << "Critiria = " << crit << endl;
        throw std::underflow_error("block too small");
    }
    int occ_i, orb_i;
    if (occ_crit < emp_crit)   // occupied
    {
        occ_i = 1;
        orb_i = 0;
    }
    else   // empty
    {
        occ_i = 0;
        orb_i = N-1;
    }
    auto tmp = U(0,0);
    typedef typename std::conditional<is_same_v<decltype(tmp),double>, Vector, CVector>::type VecT;
    return make_tuple(occ_i, VecT(column(U,orb_i)));
}

template <typename MatT>
ITensor to_many_body_rot (const Index& s1, const Index& s2, const MatT& rot)
// Spinless (spin polarized) single-body rotational gate. IQIndices must be spinless fermion.
{
    Index sP1 = prime(s1);
    Index sP2 = prime(s2);
    ITensor res(dag(s1),dag(s2),sP1,sP2);

    res.set (s1(1),s2(1),sP1(1),sP2(1),1.0);
    res.set (s1(1),s2(2),sP1(1),sP2(2),rot(0,0));
    res.set (s1(2),s2(1),sP1(2),sP2(1),rot(1,1));
    res.set (s1(1),s2(2),sP1(2),sP2(1),rot(0,1));
    res.set (s1(2),s2(1),sP1(1),sP2(2),rot(1,0));
    res.set (s1(2),s2(2),sP1(2),sP2(2),1.0);

    return res;
}

ITensor Hub_combiner (const Index& iihub, const Index& iiup, const Index& iidn)
// Split or combine between index {0,up,dn,2} and indices {0,up},{0,dn}
{
    Index si = dag(iihub);
    ITensor del (si,iiup,iidn);
    del.set (si=1,iiup=1,iidn=1, 1.);
    del.set (si=2,iiup=2,iidn=1, 1.);
    del.set (si=3,iiup=1,iidn=2, 1.);
    del.set (si=4,iiup=2,iidn=2, 1.);
    return del;
}

Index Spin_polar_index (const string& spin_str)
// <spin_str> can be "up" or "dn"
{
    if (spin_str == "up") {
        return Index (QN({"Sz",0},{"Nf",0,-1}),1,
                      QN({"Sz",1},{"Nf",1,-1}),1);
    }
    else if (spin_str == "dn") {
        return Index (QN({"Sz",0},{"Nf",0,-1}),1,
                      QN({"Sz",-1},{"Nf",1,-1}),1);
    }
    cout << "Error: spin_polar_index: invalid spin_str" << endl;
    throw;
}

ITensor Swap (Index ii1, Index ii2)
{
    ii1 = dag(ii1);
    ii2 = dag(ii2);
    Index ii1_pr = prime(dag(ii1)),
          ii2_pr = prime(dag(ii2));

    auto pfs1 = iut::get_fermion_parity (ii1);
    auto pfs2 = iut::get_fermion_parity (ii2);

    ITensor sw (ii1, ii2, ii1_pr, ii2_pr);
    for(int i = 1; i <= dim(ii1); ++i)
        for(int j = 1; j <= dim(ii2); ++j)
        {
            if (pfs1.at(i) == 1 and pfs2.at(j) == 1)
                sw.set (ii1(i), ii2(j), ii1_pr(i), ii2_pr(j),-1.);
            else
                sw.set (ii1(i), ii2(j), ii1_pr(i), ii2_pr(j),1.);
        }
    return sw;
}

template <typename MatT>
ITensor Rot_hub (const Index& s1, const Index& s2, const MatT& rot_up_mat, const MatT& rot_dn_mat)
// Single-body rotational gate for Hubbard indices.
{
    Index iiup1 = Spin_polar_index ("up"),
          iiup2 = Spin_polar_index ("up"),
          iidn1 = Spin_polar_index ("dn"),
          iidn2 = Spin_polar_index ("dn");

    ITensor rot_up = to_many_body_rot (iiup1, iiup2, rot_up_mat),
            rot_dn = to_many_body_rot (iidn1, iidn2, rot_dn_mat),
            swap   = Swap (iidn1, iiup2);

    rot_up.prime();
    rot_dn.prime();

    ITensor re = swap * rot_up * rot_dn;
    re.mapPrime (1,0);
    swap.mapPrime (0,2);
    re *= swap;
    re.mapPrime (2,1);

    ITensor comb1 = Hub_combiner (s1, iiup1, iidn1),
            comb2 = Hub_combiner (s2, iiup2, iidn2);

    re *= comb1;
    re *= comb2;
    re *= dag(prime(comb1));
    re *= dag(prime(comb2));

    return re;
}

template <typename MatT>
void Apply_gates (MPS& psi, const vector<int>& pos, const vector<MatT>& rots_up, const vector<MatT>& rots_dn, Args const& args)
// Apply the gates to <psi>. <psi> is changed after the function call.
{
    if (pos.size() != rots_up.size() || pos.size() != rots_dn.size()) {
        cout << "Error: Apply_gates: positions and angles sizes not match" << endl;
        cout << "       " << pos.size() << ", " << rots_up.size() << ", " << rots_dn.size() << endl;
        throw;
    }

    auto sites = Electron(siteInds(psi));
    //for(size_t i = 0; i < pos.size(); i++)
    for(int i = pos.size()-1; i >= 0; i--)
    {
        int p = pos.at(i);
        MatT rot_up = conj(transpose(rots_up.at(i)));
        MatT rot_dn = conj(transpose(rots_dn.at(i)));

        ITensor gate = Rot_hub (sites(p), sites(p+1), rot_up, rot_dn);
        psi.position (p, args);
        applyGate (gate, psi, args);
    }
}

template <typename MatT>
void apply_one_body_rot_mat (const vector<MatT>& rots, const vector<int>& pos, MatT& lambda)
{
    for(int i = 0; i < pos.size(); i++)
    {
        int j = pos.at(i);
        auto const& rot = rots.at(i);
        columns(lambda,j-1,j+1) &= columns(lambda,j-1,j+1) * conj(transpose(rot));
        rows(lambda,j-1,j+1) &= rot * MatT(rows(lambda,j-1,j+1));
    }
}

template <typename MatT>
void apply_one_body_rot_mat_left (const vector<MatT>& rots, const vector<int>& pos, MatT& lambda)
{
    for(int i = 0; i < pos.size(); i++)
    {
        int j = pos.at(i);
        auto const& rot = rots.at(i);
        rows(lambda,j-1,j+1) &= rot * MatT(rows(lambda,j-1,j+1));
    }
}

template <typename MatT>
void apply_one_body_rot_mat_right (const vector<MatT>& rots, const vector<int>& pos, MatT& lambda)
{
    for(int i = 0; i < pos.size(); i++)
    {
        int j = pos.at(i);
        auto const& rot = rots.at(i);
        columns(lambda,j-1,j+1) &= columns(lambda,j-1,j+1) * conj(transpose(rot));
    }
}

template <typename VecT>
void apply_one_body_rot_vec (const vector<MatT>& rots, const vector<int>& pos, VecT& v)
{
    for(int i = 0; i < pos.size(); i++)
    {
        int j = pos.at(i);
        auto const& rot = rots.at(i);
        subVector(v,j-1,j+1) &= rot * subVector(v,j-1,j+1);
    }
}

string get_state_str (int occ_up, int occ_dn)
{
    string st;
    if (occ_up and occ_dn)
        st = "UpDn";
    else if (occ_up and !occ_dn)
        st = "Up";
    else if (!occ_up and occ_dn)
        st = "Dn";
    else
        st = "Emp";
    return st;
}

template <typename MatT>
void check_hermitian (const MatT& M)
{
    auto D = conj(transpose(M)) - M;
    Real nn = norm(D);
    mycheck (nn < 1e-12, "check hermitian failed");
}

template <typename MatT>
void check_denmat (int Np, const MatT& lambda)
{
    check_hermitian (lambda);
    MatT U;
    Vector occs;
    diagHermitian (lambda, U, occs);

    Real n = 0.,
         nn = 0.;
    for(int i = 0; i < occs.size(); i++)
    {
        mycheck (occs(i) > 0. or abs(occs(i)) < 1e-12, "negative occupasion");
        n += occs(i);
        if constexpr (is_same_v<decltype(lambda(i,i)),Real>)
            nn += lambda(i,i);
        else
            nn += lambda(i,i).real();
    }
    mycheck (abs(n-Np) < 1e-2, "particle number check failed");
    mycheck (abs(n-nn) < 1e-2, "particle number check failed2");
}

void update_Np (const string& st, int& Np_up, int& Np_dn)
{
    if (st == "UpDn")
    {
        Np_up++;
        Np_dn++;
    }
    else if (st == "Up")
    {
        Np_up++;
    }
    else if (st == "Dn")
    {
        Np_dn++;
    }
}

template <typename MatT>
MPS GaussianMPS (const MatT& phi_up, const MatT& phi_dn, int block_size, Real crit, Args const& args=Args::global())
{
    int N = nrows (phi_up);
    Electron sites (N);

    // Get  one-body density matrix
    MatT lambda_up = phi_up * conj(transpose(phi_up));
    MatT lambda_dn = phi_dn * conj(transpose(phi_dn));

    int Np_up = ncols (phi_up);
    int Np_dn = ncols (phi_dn);
    //check_denmat (Np_up, lambda_up);
    //check_denmat (Np_dn, lambda_dn);

    //Electron sites (N);
    InitState init (sites);

    vector<MatT> rot_up_all, rot_dn_all;
    vector<int> pos_all;
    int Np_up_check=0, Np_dn_check=0;
    for(int i1 = 1; i1 < N; i1++)
    {
        // Get the sub-matrix of lambda
        int i2 = i1 + block_size - 1;
        if (i2 > N)
            i2 = N;
        auto sub_lambda_up = MatT (subMatrix (lambda_up, i1-1, i2, i1-1, i2));
        auto sub_lambda_dn = MatT (subMatrix (lambda_dn, i1-1, i2, i1-1, i2));
        //check_denmat (Np_up-Np_up_check, MatT (subMatrix (lambda_up, i1-1, N, i1-1, N)));
        //check_denmat (Np_dn-Np_dn_check, MatT (subMatrix (lambda_dn, i1-1, N, i1-1, N)));

        // Diagonalize the sub-matrix
        //MatT u_up, u_dn;
        MatT u_up, u_dn;
        Vector occs_up, occs_dn;
        diagHermitian (sub_lambda_up, u_up, occs_up);
        diagHermitian (sub_lambda_dn, u_dn, occs_dn);

        // Find an occupied or empty orbital
        auto [occ_up, orb_up] = get_approx_occ (occs_up, u_up, crit);
        auto [occ_dn, orb_dn] = get_approx_occ (occs_dn, u_dn, crit);

        // Set the product state at site i1
        auto st = get_state_str (occ_up, occ_dn);
        init.set (i1, st);
        update_Np (st, Np_up_check, Np_dn_check);

        // Get the angles to rotate the basis
        auto rots_up = Get_rot_mats (orb_up);
        auto rots_dn = Get_rot_mats (orb_dn);

        // Get the position to apply the gates
        vector<int> pos;
        for(int i = i2-1; i >= i1; i--)
        {
            pos.push_back (i);
        }

        // Rotate lambda
        apply_one_body_rot_mat (rots_up, pos, lambda_up);
        apply_one_body_rot_mat (rots_dn, pos, lambda_dn);

        pos_all.insert (pos_all.end(), pos.begin(), pos.end());
        rot_up_all.insert (rot_up_all.end(), rots_up.begin(), rots_up.end());
        rot_dn_all.insert (rot_dn_all.end(), rots_dn.begin(), rots_dn.end());
    }
    // Last site
    Real occ_up, occ_dn;
    auto tmp = lambda_up(N-1,N-1);
    if constexpr (is_same_v<decltype(tmp),Real>)
    {
        occ_up = lambda_up(N-1,N-1);
        occ_dn = lambda_dn(N-1,N-1);
    }
    else
    {
        occ_up = lambda_up(N-1,N-1).real();
        occ_dn = lambda_dn(N-1,N-1).real();
    }
    int n_up = get_approx_occ (occ_up, crit);
    int n_dn = get_approx_occ (occ_dn, crit);
    auto st = get_state_str (n_up, n_dn);
    init.set (N, st);
    update_Np (st, Np_up_check, Np_dn_check);
    if (Np_up_check != Np_up or Np_dn_check != Np_dn)
    {
        cout << "particle number not match: "<< Np_up_check << " " << Np_up << " | " << Np_dn_check << " " << Np_dn << endl;
        throw;
    }

    // Initialize MPS as a product state
    auto psi = MPS (init);

    // Apply the gates to rotate the MPS to new basis
    cout << "number of gates = " << pos_all.size() << endl;
    Apply_gates (psi, pos_all, rot_up_all, rot_dn_all, args);

    return psi;
}
#endif
