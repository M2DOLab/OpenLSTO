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

#ifndef _HOLE_H
#define _HOLE_H

/*! \file Hole.h
    \brief A simple circular hole data type.
 */

//! Class for handling circular holes.
class Hole
{
public:
    //! Default constructor.
    Hole();

    //! Constructor.
    /*! \param x
            The x coordinate of the hole.

        \param y
            The x coordinate of the hole.

        \param r
            The radius of the hole.
     */
    Hole(double, double, double);

    //! Constructor.
    /*! \param coord_
            The x-y coordinates of the hole.

        \param r
            The radius of the hole.
     */
    Hole(Coord&, double);

    /// Coordinates.
    Coord coord;

    /// Radius.
    double r;
};

#include "../src/hole.cpp"

#endif  /* _HOLE_H */
