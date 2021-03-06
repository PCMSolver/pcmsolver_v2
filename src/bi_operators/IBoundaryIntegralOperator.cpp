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

#include "IBoundaryIntegralOperator.hpp"

#include <Eigen/Core>

#include "cavity/ICavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "utils/Timer.hpp"

namespace pcm {
using cavity::Element;
using utils::symmetryBlocking;

Eigen::MatrixXd IBoundaryIntegralOperator::computeS(
    const ICavity & cav,
    const IGreensFunction & gf) const {
  Eigen::MatrixXd biop = computeS_impl(cav.elements(), gf);
  Eigen::LDLT<Eigen::MatrixXd> Sldlt(biop);
  if (!Sldlt.isPositive()) {
    std::ostringstream errmsg;
    errmsg << "S matrix is not positive-definite!" << std::endl;
    errmsg << "Consider changing the average area of the cavity finite elements."
           << std::endl;
    errmsg
        << "Please report this issue: https://github.com/PCMSolver/pcmsolver/issues"
        << std::endl;
    PCMSOLVER_ERROR(errmsg.str());
  }
  // Perform symmetry blocking
  // The total size of the cavity
  int cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  if (cav.pointGroup().nrGenerators() != 0) {
    timer::timerON("Symmetry blocking");
    symmetryBlocking(biop, cavitySize, dimBlock, nrBlocks);
    timer::timerOFF("Symmetry blocking");
  }
  return biop;
}

Eigen::MatrixXd IBoundaryIntegralOperator::computeD(
    const ICavity & cav,
    const IGreensFunction & gf) const {
  Eigen::MatrixXd biop = computeD_impl(cav.elements(), gf);
  // Perform symmetry blocking
  // The total size of the cavity
  int cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  if (cav.pointGroup().nrGenerators() != 0) {
    timer::timerON("Symmetry blocking");
    symmetryBlocking(biop, cavitySize, dimBlock, nrBlocks);
    timer::timerOFF("Symmetry blocking");
  }
  return biop;
}
} // namespace pcm
