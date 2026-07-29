#include "matrix/assembler_base.h"
#include "matrix/assembler_eigen_semi_implicit_twophase.h"

#include "matrix/cg_eigen.h"
#include "matrix/solver_base.h"

// Asssembler 2D
static Register<mpm::AssemblerBase<2>,
                mpm::AssemblerEigenSemiImplicitTwoPhase<2>>
    AssemblerEigenSemiImplicitTwoPhase2d("EigenSemiImplicitTwoPhase2D");
// Asssembler 3D
static Register<mpm::AssemblerBase<3>,
                mpm::AssemblerEigenSemiImplicitTwoPhase<3>>
    AssemblerEigenSemiImplicitTwoPhase3d("EigenSemiImplicitTwoPhase3D");
// Solver 2D
static Register<mpm::SolverBase<2>, mpm::CGEigen<2>, unsigned, double>
    solvereigencg2d("EigenCG2D");
// Solver 3D
static Register<mpm::SolverBase<3>, mpm::CGEigen<3>, unsigned, double>
    solvereigencg3d("EigenCG3D");