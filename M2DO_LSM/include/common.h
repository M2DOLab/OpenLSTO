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

#ifndef _COMMON_H
#define _COMMON_H

/*! \file Common.h
    \brief Common data types.
 */

//! Two-dimensional coordinate information.
struct Coord
{
    double x;   //!< The x coordinate.
    double y;   //!< The y coordinate.
};

//! \brief A container for storing information associated with a boundary point.
struct BoundaryPoint
{
    Coord coord;                        //!< Coordinate of the boundary point.
    Coord normal;                       //!< Inward pointing normal vector.
    double length;                      //!< Integral length of the boundary point.
    double velocity;                    //!< Normal velocity (positive acts inwards).
    double negativeLimit;               //!< Movement limit in negative direction (inwards).
    double positiveLimit;               //!< Movement limit in positive direction (outwards).
    bool isDomain;                      //!< Whether the point lies close to the domain boundary.
    bool isFixed;                       //!< Whether the point fixed.
    unsigned int nSegments;             //!< The number of boundary segments that a point belongs to.
    unsigned int segments[2];           //!< The indices of the two segments to which a point belongs.
    unsigned int neighbours[2];         //!< The indices of the neighbouring points.
    unsigned int nNeighbours;           //!< The number of neighbouring boundary points.
    std::vector<double> sensitivities;  //!< Objective and constraint sensitivities.
};

//! \brief A container for storing information associated with a boundary segment.
struct BoundarySegment
{
    unsigned int start;                 //!< Index of start point.
    unsigned int end;                   //!< Index of end point.
    unsigned int element;               //!< The element cut by the boundary segment.
    double length;                      //!< Length of the boundary segment.
    double weight;                      //!< Weighting factor for boundary segment.
};

#endif  /* _COMMON_H */
