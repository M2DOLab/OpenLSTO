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
