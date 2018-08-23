/*! \file Sensitivity.cpp
    \brief A class for calculating finite-difference boundary point sensitivities.
 */

Sensitivity::Sensitivity(double delta_) : delta(delta_)
{
}

double Sensitivity::computeSensitivity(BoundaryPoint& point, SensitivityCallback& callback)
{
    // Store the initial boundary point coordinates.
    Coord coord = point.coord;

    // Displace the boundary point in the positive normal direction.
    point.coord.x = coord.x + delta*point.normal.x;
    point.coord.y = coord.y + delta*point.normal.y;

    // Compute the new value of the function.
    double f1 = callback(point);

    // Displace the boundary point in the negative normal direction.
    point.coord.x = coord.x - delta*point.normal.x;
    point.coord.y = coord.y - delta*point.normal.y;

    // Compute the new value of the function.
    double f2 = callback(point);

    // Compute the finite-difference derivative.
    double sens = (f1 - f2) / (2.0 * delta);

    // Divide by boundary point length (sensitivity per unit length).
    sens /= point.length;

    // Reset boundary point coordinates.
    point.coord = coord;

    // Return sensitivity.
    return sens;
}
