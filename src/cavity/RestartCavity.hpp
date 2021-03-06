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

#include <iosfwd>
#include <string>

/*! \file RestartCavity.hpp */

namespace pcm {
struct CavityData;
} // namespace pcm

#include "ICavity.hpp"

namespace pcm {
namespace cavity {
/*! \class RestartCavity
 *  \brief A class for Restart cavity.
 *  \author Roberto Di Remigio
 *  \date 2014
 */
class RestartCavity final : public ICavity {
public:
  RestartCavity(const std::string & _fname) : file(_fname) { makeCavity(); }
  virtual void makeCavity() override { loadCavity(file); }

private:
  std::string file;
  virtual std::ostream & printCavity(std::ostream & os) override;
};

ICavity * createRestartCavity(const CavityData & data);
} // namespace cavity
} // namespace pcm
