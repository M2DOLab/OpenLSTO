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
