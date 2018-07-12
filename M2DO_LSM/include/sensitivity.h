/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+lsm@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _SENSITIVITY_H
#define _SENSITIVITY_H

/*! \file Sensitivity.h
    \brief A class for calculating finite-difference boundary point sensitivities.
 */

//! Calculate the value of a function for a small displacement of a boundary point.
/*! \param point
        A reference to the boundary point.

    \return
        The value of the function.
 */
typedef std::function<double (const BoundaryPoint&)> SensitivityCallback;

/*! \brief A class for calculating finite-difference boundary point sensitivities.
*/
class Sensitivity
{
public:
    //! Constructor.
    /*! \param delta_
            The finite-difference derivative length (perturbation).
    */
    Sensitivity(double delta_ = 1e-4);

    //! Compute finite-difference sensitivity for an arbitrary function.
    /*! \param point
            A reference to the boundary point.

        \param callback
            A reference to the sensitivity callback function.

        \return
            The finite-difference sensitivity.
     */
    double computeSensitivity(BoundaryPoint&, SensitivityCallback&);

private:
    /// The finite-difference derivative length (in units of the grid spacing).
    double delta;
};

#include "../src/sensitivity.cpp"

#endif	/* _SENSITIVITY_H */
