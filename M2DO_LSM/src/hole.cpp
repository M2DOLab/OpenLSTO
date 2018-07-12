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

/*! \file Hole.cpp
    \brief A simple circular hole data type.
 */

Hole::Hole() {}

Hole::Hole(double x, double y, double r) : r(r)
{
    coord.x = x;
    coord.y = y;
}

Hole::Hole(Coord& coord_, double r) : r(r)
{
    coord.x = coord_.x;
    coord.y = coord_.y;
}
