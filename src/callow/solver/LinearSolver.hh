//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearSolver.hh
 *  @author robertsj
 *  @date   Sep 13, 2012
 *  @brief  LinearSolver class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_LINEARSOLVER_HH_
#define callow_LINEARSOLVER_HH_

#include "callow/callow_config.hh"
#include "callow/utils/CallowDefinitions.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/preconditioner/Preconditioner.hh"
#include "utilities/SP.hh"
#include <cstdio>
#include <string>
#include <vector>

namespace callow
{

/**
 *  @class LinearSolver
 *  @brief Base class for iterative linear solvers
 *
 *  We solve problems of the form
 *  @f[
 *      \mathbf{A}x = b
 *  @f]
 *  via iterative methods.  A system is "solved" when the
 *  norm of the residual is small enough or some maximum
 *  iteration count is reached.  The residual is defined
 *  @f[
 *      \mathbf{A}x - b
 *  @f]
 *  and its norm is
 *  @f[
 *    r =  || \mathbf{A}x - b ||
 *  @f]
 *  By default, the L2 norm is used, though L1 and Linf
 *  are also recorded can can be used.  This represents
 *  the absolute norm.  Sometimes it makes more sense to
 *  check the residual with respect to the initial norm,
 *  in which case the relative norm is
 *  @f[
 *      r_n / r_0 = || \mathbf{A}x^{n} - b || / || \mathbf{A}x^{0} - b ||
 *  @f]
 *  The iteration terminates when
 *  @f[
 *      r_n < \mathrm{max} ( \tau_{\mathrm{rel}} r_0, \tau_{\mathrm{abs}} )
 *  @f]
 *
 *  The solvers currently implemented are
 *    - Richardson
 *    - Jacobi
 *    - Gauss-Seidel
 *    - GMRES(m)
 *  along with Jacobi and ILU0 preconditioners.  If PETSc is enabled,
 *  all of its solvers are potentially available.
 *
 *  Note, some linear solvers require that the matrix provides L, U, and
 *  D operations.  Here, we simply require those solvers to have a \ref
 *  Matrix operator or subclasses so that the elements can be accessed
 *  directly.
 *
 */

class LinearSolver
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum pcside
  {
    NONE, LEFT, RIGHT
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<LinearSolver>           SP_solver;
  typedef typename MatrixBase::SP_matrix               SP_matrix;
  typedef typename Preconditioner::SP_preconditioner   SP_preconditioner;
  typedef typename Vector::SP_vector                   SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  LinearSolver(const double atol,
               const double rtol,
               const int    maxit = 100,
               std::string  name = "solver");

  virtual ~LinearSolver(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  Sets the operators for the linear system to solve.
   *
   *  @param A      linear operator
   *  @param P      optional preconditioning process
   *  @param side   specifies on what side of A the preconditioner operates
   */
  virtual void set_operators(SP_matrix A,
                             SP_preconditioner P = SP_preconditioner(0),
                             const int side = LEFT);

  /**
   *  Get the preconditioner.  This allows the client to build, change, etc.
   */
  SP_preconditioner preconditioner()
  {
    return d_P;
  }

  /**
   *  @param atol   absolute tolerance (||r_n|| < atol)
   *  @param rtol   relative tolerance (||r_n|| < rtol * ||r_0||)
   *  @param maxit  maximum iterations (n < maxit)
   */
  void set_tolerances(const double atol, const double rtol, const int maxit);

  /**
   *  Print residual norms and other diagonostic information.
   *
   *  @param v  monitor via stdout
   */
  void set_monitor_level(const int v)
  {
    d_monitor_level = v;
  }

  /**
   *  Turn on monitoring of diverging iterations.
   */
  void set_monitor_diverge(const bool v)
  {
    d_monitor_diverge = v;
  }

  /**
   *  @param b  right hand side
   *  @param x  unknown vector
   */
  int solve(const Vector &b, Vector &x);

  /// return the residual norms
  std::vector<double> residual_norms()
  {
    return d_residual;
  }

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  std::string d_name;
  double d_absolute_tolerance;
  double d_relative_tolerance;
  int d_maximum_iterations;
  std::vector<double> d_residual;
  int d_number_iterations;
  SP_matrix d_A;
  /// preconditioner
  SP_preconditioner d_P;
  /// preconditioner side (one only, and left by default)
  int d_pc_side;
  int d_monitor_level;
  bool d_monitor_diverge;
  int d_norm_type;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// print out iteration and residual for initial
  virtual bool monitor_init(double r);

  // print out iteration and residual
  virtual bool monitor(int it, double r);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  virtual void solve_impl(const Vector &b, Vector &x) = 0;

private:

  int d_status;

};

} // end namespace callow

#include "LinearSolver.i.hh"

#endif /* callow_LINEARSOLVER_HH_ */