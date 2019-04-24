/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

#include "catch.hpp"

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/Collocation.hpp"
#include "cavity/TsLessCavity.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "solver/CPCMSolver.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::TsLessCavity;
using green::UniformDielectric;
using green::Vacuum;
using solver::CPCMSolver;

SCENARIO("Test solver for the C-PCM for a point charge and a TsLess cavity",
         "[solver][cpcm][cpcm_tsless-point]") {
  GIVEN("An isotropic environment modelled as a perfect conductor and a point "
        "charge") {
    double permittivity = 78.39;
    Vacuum<> gfInside;
    UniformDielectric<> gfOutside(permittivity);
    bool symm = true;
    double correction = 0.0;

    Collocation S;

    double charge = 8.0;
    double totalASC = -charge * (permittivity - 1) / (permittivity + correction);

    /*! \class CPCMSolver
     *  \test \b pointChargeTsLess tests CPCMSolver using a point charge with a
     * TsLess cavity The point charge is at the origin.
     */
    WHEN("the point charge is located at the origin") {
      Molecule point = dummy<0>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      double minDistance = 0.1;
      int derOrder = 4;
      TsLessCavity cavity =
          TsLessCavity(point, area, probeRadius, minRadius, minDistance, derOrder);
      cavity.saveCavity("tsless_point.npz");

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gfInside, gfOutside, S);

      size_t size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (size_t i = 0; i < size; ++i) {
        INFO("fake_mep(" << i << ") = " << fake_mep(i));
      }

      THEN("the apparent surface charge is") {
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        double totalFakeASC = fake_asc.sum();

        for (size_t i = 0; i < size; ++i) {
          INFO("fake_asc(" << i << ") = " << fake_asc(i));
        }

        CAPTURE(totalASC);
        CAPTURE(totalFakeASC);
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeShiftedTsLess tests CPCMSolver using a point charge with
     * a TsLess cavity The point charge is away from the origin.
     */
    AND_WHEN("the point charge is located away from the origin") {
      Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
      Molecule point = dummy<0>(2.929075493, origin);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      double minDistance = 0.1;
      int derOrder = 4;
      TsLessCavity cavity =
          TsLessCavity(point, area, probeRadius, minRadius, minDistance, derOrder);

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gfInside, gfOutside, S);

      size_t size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge, origin);
      for (size_t i = 0; i < size; ++i) {
        INFO("fake_mep(" << i << ") = " << fake_mep(i));
      }

      THEN("the surface charge is") {
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        double totalFakeASC = fake_asc.sum();

        for (size_t i = 0; i < size; ++i) {
          INFO("fake_asc(" << i << ") = " << fake_asc(i));
        }

        CAPTURE(totalASC);
        CAPTURE(totalFakeASC);
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }
  }
}
