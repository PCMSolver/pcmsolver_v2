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

#pragma once

#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
class ICavity;
namespace cavity {
class Element;
} // namespace cavity
class IGreensFunction;
} // namespace pcm

namespace pcm {
class IBoundaryIntegralOperator {
public:
  virtual ~IBoundaryIntegralOperator() {}
  /*! Computes the matrix representation of the single layer operator
   *  \param[in] cav the discretized cavity
   *  \param[in] gf  a Green's function
   */
  Eigen::MatrixXd computeS(const ICavity & cav, const IGreensFunction & gf) const;
  /*! Computes the matrix representation of the double layer operator
   *  \param[in] cav the discretized cavity
   *  \param[in] gf  a Green's function
   */
  Eigen::MatrixXd computeD(const ICavity & cav, const IGreensFunction & gf) const;

private:
  /*! Computes the matrix representation of the single layer operator
   *  \param[in] elems list of finite elements of the discretized cavity
   *  \param[in] gf  a Green's function
   */
  virtual Eigen::MatrixXd computeS_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const = 0;
  /*! Computes the matrix representation of the double layer operator
   *  \param[in] elems list of finite elements of the discretized cavity
   *  \param[in] gf  a Green's function
   */
  virtual Eigen::MatrixXd computeD_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const = 0;
};
} // namespace pcm
