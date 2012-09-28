//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolver.i.hh
 * \brief  LinearSolver inline member definitions
 * \author Jeremy Roberts
 * \date   Sep 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef LINEARSOLVER_I_HH_
#define LINEARSOLVER_I_HH_


namespace callow
{

inline int LinearSolver::solve(const Vector &b, Vector &x)
  {
    Require(x.size() == b.size());
    Require(x.size() == d_A->number_rows());

    d_status = MAXIT;
    solve_impl(b, x);
    if (d_status ==  MAXIT)
    {
      printf("*** %s did not converge within the maximum number of iterations\n",
             d_name.c_str());
    }
    return d_status;
  }

} // end namespace callow

#endif // LINEARSOLVER_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file LinearSolver.i.hh
//---------------------------------------------------------------------------//