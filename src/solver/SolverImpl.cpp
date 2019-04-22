/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#include "SolverImpl.hpp"

#include <cmath>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>

#include "bi_operators/IBoundaryIntegralOperator.hpp"
#include "cavity/Element.hpp"
#include "cavity/ICavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/Timer.hpp"

namespace pcm {
namespace solver {
namespace detail {
Eigen::MatrixXd anisotropicTEpsilon(const ICavity & cav,
                                    const IGreensFunction & gf_i,
                                    const IGreensFunction & gf_o,
                                    const IBoundaryIntegralOperator & op) {
  timer::timerON("Computing SI");
  Eigen::MatrixXd SI = op.computeS(cav, gf_i);
  timer::timerOFF("Computing SI");
  timer::timerON("Computing DI");
  Eigen::MatrixXd DI = op.computeD(cav, gf_i);
  timer::timerOFF("Computing DI");
  timer::timerON("Computing SE");
  Eigen::MatrixXd SE = op.computeS(cav, gf_o);
  timer::timerOFF("Computing SE");
  timer::timerON("Computing DE");
  Eigen::MatrixXd DE = op.computeD(cav, gf_o);
  timer::timerOFF("Computing DE");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cav.size(), cav.size());

  timer::timerON("Assemble T matrix");
  Eigen::MatrixXd T = ((2 * M_PI * Id - DE * a) * SI +
                       SE * (2 * M_PI * Id + a * DI.adjoint().eval()));
  timer::timerOFF("Assemble T matrix");

  return T;
}

Eigen::MatrixXd isotropicTEpsilon(const ICavity & cav,
                                  const IGreensFunction & gf_i,
                                  double epsilon,
                                  const IBoundaryIntegralOperator & op) {
  timer::timerON("Computing SI");
  Eigen::MatrixXd SI = op.computeS(cav, gf_i);
  timer::timerOFF("Computing SI");
  timer::timerON("Computing DI");
  Eigen::MatrixXd DI = op.computeD(cav, gf_i);
  timer::timerOFF("Computing DI");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cav.size(), cav.size());

  double fact = (epsilon + 1.0) / (epsilon - 1.0);
  timer::timerON("Assemble T matrix");
  Eigen::MatrixXd T = (2 * M_PI * fact * Id - DI * a) * SI;
  timer::timerOFF("Assemble T matrix");

  return T;
}

Eigen::MatrixXd anisotropicRinfinity(const ICavity & cav,
                                     const IGreensFunction & gf_i,
                                     const IGreensFunction & gf_o,
                                     const IBoundaryIntegralOperator & op) {
  timer::timerON("Computing SI");
  Eigen::MatrixXd SI = op.computeS(cav, gf_i);
  timer::timerOFF("Computing SI");
  timer::timerON("Computing DI");
  Eigen::MatrixXd DI = op.computeD(cav, gf_i);
  timer::timerOFF("Computing DI");
  timer::timerON("Computing SE");
  Eigen::MatrixXd SE = op.computeS(cav, gf_o);
  timer::timerOFF("Computing SE");
  timer::timerON("Computing DE");
  Eigen::MatrixXd DE = op.computeD(cav, gf_o);
  timer::timerOFF("Computing DE");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cav.size(), cav.size());

  timer::timerON("Assemble R matrix");
  Eigen::MatrixXd R =
      ((2 * M_PI * Id - DE * a) - SE * SI.llt().solve((2 * M_PI * Id - DI * a)));
  timer::timerOFF("Assemble R matrix");

  return R;
}

Eigen::MatrixXd isotropicRinfinity(const ICavity & cav,
                                   const IGreensFunction & gf_i,
                                   const IBoundaryIntegralOperator & D) {
  timer::timerON("Computing DI");
  Eigen::MatrixXd DI = D.computeD(cav, gf_i);
  timer::timerOFF("Computing DI");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cav.size(), cav.size());

  timer::timerON("Assemble R matrix");
  Eigen::MatrixXd R = (2 * M_PI * Id - DI * a);
  timer::timerOFF("Assemble R matrix");

  return R;
}

Eigen::MatrixXd anisotropicIEFMatrix(const ICavity & cav,
                                     const IGreensFunction & gf_i,
                                     const IGreensFunction & gf_o,
                                     const IBoundaryIntegralOperator & op) {
  Eigen::MatrixXd T = anisotropicTEpsilon(cav, gf_i, gf_o, op);
  Eigen::MatrixXd R = anisotropicRinfinity(cav, gf_i, gf_o, op);

  timer::timerON("Assemble T^-1R matrix");
  Eigen::MatrixXd fullPCMMatrix = T.partialPivLu().solve(R);
  timer::timerOFF("Assemble T^-1R matrix");

  return fullPCMMatrix;
}

Eigen::MatrixXd isotropicIEFMatrix(const ICavity & cav,
                                   const IGreensFunction & gf_i,
                                   double epsilon,
                                   const IBoundaryIntegralOperator & op) {
  Eigen::MatrixXd T = isotropicTEpsilon(cav, gf_i, epsilon, op);
  Eigen::MatrixXd R = isotropicRinfinity(cav, gf_i, op);

  timer::timerON("Assemble T^-1R matrix");
  Eigen::MatrixXd fullPCMMatrix = T.partialPivLu().solve(R);
  timer::timerOFF("Assemble T^-1R matrix");

  return fullPCMMatrix;
}
} // namespace detail
} // namespace solver
} // namespace pcm
