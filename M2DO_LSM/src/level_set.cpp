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

/*! \file LevelSet.cpp
    \brief A class for the level set function.
 */

LevelSet::LevelSet(Mesh& mesh_, double moveLimit_, unsigned int bandWidth_, bool isFixed_) :
    moveLimit(moveLimit_),
    mesh(mesh_),
    bandWidth(bandWidth_),
    isFixed(isFixed_),
    isTarget(false)
{
    int size = 0.2*mesh.nNodes;

    errno = EINVAL;
    lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");
    lsm_check(((moveLimit > 0) && (moveLimit <= 1)), "Move limit must be between 0 and 1.");

    // Resize data structures.
    signedDistance.resize(mesh.nNodes);
    velocity.resize(mesh.nNodes);
    gradient.resize(mesh.nNodes);
    narrowBand.resize(mesh.nNodes);

    // Make sure that memory is sufficient (for small test systems).
    size = std::max(25, size);
    mines.resize(size);

    // Generate a Swiss cheese structure.
    initialise();

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

LevelSet::LevelSet(Mesh& mesh_, const std::vector<Hole>& holes,
    double moveLimit_, unsigned int bandWidth_, bool isFixed_) :
    moveLimit(moveLimit_),
    mesh(mesh_),
    bandWidth(bandWidth_),
    isFixed(isFixed_),
    isTarget(false)
{
    int size = 0.2*mesh.nNodes;

    errno = EINVAL;
    lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");
    lsm_check(((moveLimit > 0) && (moveLimit <= 1)), "Move limit must be between 0 and 1.");

    // Resize data structures.
    signedDistance.resize(mesh.nNodes);
    velocity.resize(mesh.nNodes);
    gradient.resize(mesh.nNodes);
    narrowBand.resize(mesh.nNodes);

    // Make sure that memory is sufficient (for small test systems).
    size = std::max(25, size);
    mines.resize(size);

    // Initialise level set function from hole array.
    initialise(holes);

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

LevelSet::LevelSet(Mesh& mesh_, const std::vector<Coord>& points,
    double moveLimit_, unsigned int bandWidth_, bool isFixed_) :
    moveLimit(moveLimit_),
    mesh(mesh_),
    bandWidth(bandWidth_),
    isFixed(isFixed_),
    isTarget(false)
{
    int size = 0.2*mesh.nNodes;

    errno = EINVAL;
    lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");
    lsm_check(((moveLimit > 0) && (moveLimit <= 1)), "Move limit must be between 0 and 1.");

    // Resize data structures.
    signedDistance.resize(mesh.nNodes);
    velocity.resize(mesh.nNodes);
    gradient.resize(mesh.nNodes);
    narrowBand.resize(mesh.nNodes);

    // Make sure that memory is sufficient (for small test systems).
    size = std::max(25, size);
    mines.resize(size);

    // Initialise level set function from point array.
    initialise(points);

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

LevelSet::LevelSet(Mesh& mesh_, const std::vector<Hole>& initialHoles,
    const std::vector<Hole>& targetHoles, double moveLimit_, unsigned int bandWidth_, bool isFixed_) :
    moveLimit(moveLimit_),
    mesh(mesh_),
    bandWidth(bandWidth_),
    isFixed(isFixed_),
    isTarget(true)
{
    int size = 0.2*mesh.nNodes;

    errno = EINVAL;
    lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");
    lsm_check(((moveLimit > 0) && (moveLimit <= 1)), "Move limit must be between 0 and 1.");

    // Resize data structures.
    signedDistance.resize(mesh.nNodes);
    velocity.resize(mesh.nNodes);
    gradient.resize(mesh.nNodes);
    target.resize(mesh.nNodes);
    narrowBand.resize(mesh.nNodes);

    // Make sure that memory is sufficient (for small test systems).
    size = std::max(25, size);
    mines.resize(size);

    // Initialise level set function the target holes vector.
    initialise(targetHoles);
    reinitialise();

    // Copy into target signed distance function.
    target = signedDistance;

    // Initialise level set function the initial holes vector.
    initialise(initialHoles);

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

LevelSet::LevelSet(Mesh& mesh_, const std::vector<Hole>& holes,
    const std::vector<Coord>& points, double moveLimit_, unsigned int bandWidth_, bool isFixed_) :
    moveLimit(moveLimit_),
    mesh(mesh_),
    bandWidth(bandWidth_),
    isFixed(isFixed_),
    isTarget(true)
{
    int size = 0.2*mesh.nNodes;

    errno = EINVAL;
    lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");
    lsm_check(((moveLimit > 0) && (moveLimit <= 1)), "Move limit must be between 0 and 1.");

    // Resize data structures.
    signedDistance.resize(mesh.nNodes);
    velocity.resize(mesh.nNodes);
    gradient.resize(mesh.nNodes);
    target.resize(mesh.nNodes);
    narrowBand.resize(mesh.nNodes);

    // Make sure that memory is sufficient (for small test systems).
    size = std::max(25, size);
    mines.resize(size);

    // Initialise level set function from point array.
    initialise(points);
    reinitialise();

    // Copy into target signed distance function.
    target = signedDistance;

    // Initialise level set function from hole array.
    initialise(holes);

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

LevelSet::LevelSet(Mesh& mesh_, const std::vector<Coord>& initialPoints,
    const std::vector<Coord>& targetPoints, double moveLimit_, unsigned int bandWidth_, bool isFixed_) :
    moveLimit(moveLimit_),
    mesh(mesh_),
    bandWidth(bandWidth_),
    isFixed(isFixed_),
    isTarget(true)
{
    int size = 0.2*mesh.nNodes;

    errno = EINVAL;
    lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");
    lsm_check(((moveLimit > 0) && (moveLimit <= 1)), "Move limit must be between 0 and 1.");

    // Resize data structures.
    signedDistance.resize(mesh.nNodes);
    velocity.resize(mesh.nNodes);
    gradient.resize(mesh.nNodes);
    target.resize(mesh.nNodes);
    narrowBand.resize(mesh.nNodes);

    // Make sure that memory is sufficient (for small test systems).
    size = std::max(25, size);
    mines.resize(size);

    // Initialise level set function from target point array.
    initialise(targetPoints);
    reinitialise();

    // Copy into target signed distance function.
    target = signedDistance;

    // Initialise level set function from initialisation points array.
    initialise(initialPoints);

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

bool LevelSet::update(double timeStep)
{
    // Loop over all nodes in the narrow band.
    for (unsigned int i=0;i<nNarrowBand;i++)
    {
        unsigned int node = narrowBand[i];
        signedDistance[node] -= timeStep * gradient[node] * velocity[node];

        // If node is on domain boundary.
        if (mesh.nodes[node].isDomain)
        {
            // Enforce boundary condition.
            if (signedDistance[node] > 0)
                signedDistance[node] = 0;
        }

        // Reset the number of boundary points.
        mesh.nodes[node].nBoundaryPoints = 0;
    }

    // Check mine nodes.
    for (unsigned int i=0;i<nMines;i++)
    {
        // Boundary is within one grid spacing of the mine.
        if (std::abs(signedDistance[mines[i]]) < 1.0)
        {
            // Reinitialise the signed distance function.
            reinitialise();

            return true;
        }
    }

    return false;
}

void LevelSet::killNodes(const std::vector<Coord>& points)
{
    // Kill nodes in a rectangular region.
    if (points.size() == 2)
    {
        // Loop over all nodes.
        for (unsigned int i=0;i<mesh.nNodes;i++)
        {
            // Point is inside the rectangle.
            if (mesh.nodes[i].coord.x > points[0].x &&
                mesh.nodes[i].coord.y > points[0].y &&
                mesh.nodes[i].coord.x < points[1].x &&
                mesh.nodes[i].coord.y < points[1].y)
            {
                signedDistance[i] = -1e-6;
                mesh.nodes[i].isFixed = true;
            }
        }
    }

    // Mask off a piece-wise linear shape.
    else
    {
        // Loop over all nodes.
        for (unsigned int i=0;i<mesh.nNodes;i++)
        {
            // Point is inside the polygon.
            if (isInsidePolygon(mesh.nodes[i].coord, points))
            {
                signedDistance[i] = -1e-6;
                mesh.nodes[i].isFixed = true;
            }
        }
    }
}

void LevelSet::fixNodes(const std::vector<Coord>& points)
{
    // Loop over all nodes.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Point is inside the rectangle.
        if (mesh.nodes[i].coord.x > points[0].x &&
            mesh.nodes[i].coord.y > points[0].y &&
            mesh.nodes[i].coord.x < points[1].x &&
            mesh.nodes[i].coord.y < points[1].y)
        {
            mesh.nodes[i].isFixed = true;
        }
    }
}

void LevelSet::createLevelSetBoundary(const std::vector<Coord>& points)
{
    // Loop over all nodes.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Point is inside the rectangle.
        if (mesh.nodes[i].coord.x > points[0].x &&
            mesh.nodes[i].coord.y > points[0].y &&
            mesh.nodes[i].coord.x < points[1].x &&
            mesh.nodes[i].coord.y < points[1].y)
        {
            signedDistance[i] = 0;
        }
    }
}

void LevelSet::computeVelocities(const std::vector<BoundaryPoint>& boundaryPoints)
{
    // Initialise velocity (map boundary points to boundary nodes).
    initialiseVelocities(boundaryPoints);

    // Initialise fast marching method object.
    FastMarchingMethod fmm(mesh, false);

    // Reinitialise the signed distance function.
    fmm.march(signedDistance, velocity);
}

double LevelSet::computeVelocities(std::vector<BoundaryPoint>& boundaryPoints,
    double& timeStep, const double temperature, MersenneTwister& rng)
{
    // Square root of two times temperature, sqrt(2T).
    double sqrt2T = sqrt(2.0 * temperature);

    // Time step scale factor.
    double scale = 1.0;

    // Check that noise won't lead to severe CFL violation.
    if ((sqrt2T * sqrt(timeStep)) > (0.5 * moveLimit))
    {
        // Calculate time step scale factor.
        scale = (8.0 * timeStep * temperature) / (moveLimit * moveLimit);

        // Scale the time step.
        timeStep /= scale;
    }

    // Calculate noise prefactor.
    double noise = sqrt2T / sqrt(timeStep);

    // Add random noise to velocity of each boundary point.
    for (unsigned int i=0;i<boundaryPoints.size();i++)
        boundaryPoints[i].velocity += noise * rng.normal(0, 1);

    // Perform velocity extension.
    computeVelocities(boundaryPoints);

    return scale;
}

void LevelSet::computeGradients()
{
    // Compute gradient of the signed distance function using upwind finite difference.
    // This function assumes that velocities have already been calculated.

    // Reset gradients.
    std::fill(gradient.begin(), gradient.end(), 0.0);

    // std::cout << "\n\ngradient.size() = " << gradient.size() << std::endl ;
    // std::cout << "\n\nnNarrowBand = " << nNarrowBand << std::endl ;

    // Loop over all nodes in the narrow band region.
    for (unsigned int i=0;i<nNarrowBand;i++)
    {
        // Compute the nodal gradient.
        unsigned int node = narrowBand[i];
        // std::cout << "\nx" << std::endl ; 
        gradient[node] = computeGradient(node);
        // std::cout << "\ny" << std::endl ; 
    }
}

void LevelSet::reinitialise()
{
    // Initialise fast marching method object.
    FastMarchingMethod fmm(mesh, false);

    // Reinitialise the signed distance function.
    fmm.march(signedDistance);

    // Reinitialise the narrow band.
    initialiseNarrowBand();
}

void LevelSet::initialise()
{
    // Generate a swiss cheese arrangement of holes.
    // The holes have a default radius of 5.

    // Number of holes in x and y directions.
    unsigned int nx = std::round((double) mesh.width / 30);
    unsigned int ny = std::round((double) mesh.height / 30);

    // Number of holes.
    unsigned int n1 = (nx * ny);                // outer grid
    unsigned int n2 = ((nx - 1) * (ny - 1));    // inner grid
    unsigned int nHoles = n1 + n2;

    // Initialise a vector of holes.
    std::vector<Hole> holes(nHoles);

    // Hole separations.
    double dx, dy;

    // Check that mesh is large enough.
    lsm_check(((nx > 2) && (ny > 2)), "Mesh is too small for Swiss cheese initialisation.");

    // Initialise hole separations.
    dx = ((double) mesh.width / (2 * nx));
    dy = ((double) mesh.height / (2 * ny));

    // Calculate hole coordinates. (outer grid)
    for (unsigned int i=0;i<n1;i++)
    {
        // Work out x and y indices for hole.
        unsigned x = i % nx;
        unsigned int y = int(i / nx);

        // Set hole coordinates and radius.
        holes[i].coord.x = dx + (2 * x * dx);
        holes[i].coord.y = dy + (2 * y * dy);
        holes[i].r = 5;
    }

    // Calculate hole coordinates. (inner grid)
    for (unsigned int i=0;i<n2;i++)
    {
        // Work out x and y indices for hole.
        unsigned x = i % (nx - 1);
        unsigned int y = int(i / (nx -1));

        // Set hole coordinates and radius.
        holes[i + n1].coord.x = 2 * (dx + (x * dx));
        holes[i + n1].coord.y = 2 * (dy + (y * dy));
        holes[i + n1].r = 5;
    }

    // Now pass the holes to the initialise method.
    initialise(holes);

    return;

error:
    exit(EXIT_FAILURE);
}

void LevelSet::initialise(const std::vector<Hole>& holes)
{
    // First initialise the signed distance based on domain boundary.
    closestDomainBoundary();

    /* Now test signed distance against the surface of each hole.
       Update signed distance function when distance to hole surface
       is less than the current value. Since this is only done once, we
       use the simplest implementation possible.
     */

    // Loop over all nodes.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Loop over all holes.
        for (unsigned int j=0;j<holes.size();j++)
        {
            // Work out x and y distance of the node from the hole centre.
            double dx = holes[j].coord.x - mesh.nodes[i].coord.x;
            double dy = holes[j].coord.y - mesh.nodes[i].coord.y;

            // Work out distance (Pythag).
            double dist = sqrt(dx*dx + dy*dy);

            // Signed distance from the hole surface.
            dist -= holes[j].r;

            // If distance is less than current value, then update.
            if (dist < signedDistance[i])
                signedDistance[i] = dist;
        }
    }
}

void LevelSet::initialise(const std::vector<Coord>& points)
{
    // First initialise the signed distance based on domain boundary.
    closestDomainBoundary();

    /* The points define a piece-wise linear interface, i.e. they are
       assumed to be ordered (clockwise) with each pair of points
       forming a segment of the boundary. The first and last point
       in the vector should be identical, i.e. we have a closed loop
       (for an n-gon there would be n+1 points).

       The signed distance at each node is tested against each segment.

       N.B. We currently only support single closed interface.
     */

    // Loop over all nodes.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Loop over all interface segments.
        for (unsigned int j=0;j<points.size()-1;j++)
        {
            // Compute the minimum distance to the line segment j --> j+1.
            double dist = pointToLineDistance(points[j], points[j+1], mesh.nodes[i].coord);

            // If distance is less than current value, then update.
            if (dist < signedDistance[i])
                signedDistance[i] = dist;
        }

        // Invert the signed distance function if the point lies inside the polygon.
        if (isInsidePolygon(mesh.nodes[i].coord, points))
            signedDistance[i] *= -1;
    }
}

void LevelSet::closestDomainBoundary()
{
    // Initial LSF is distance from closest domain boundary.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Closest edge in x.
        unsigned int minX = std::min(mesh.nodes[i].coord.x, mesh.width - mesh.nodes[i].coord.x);

        // Closest edge in y.
        unsigned int minY = std::min(mesh.nodes[i].coord.y, mesh.height - mesh.nodes[i].coord.y);

        // Signed distance is the minimum of minX and minY;
        signedDistance[i] = double(std::min(minX, minY));
    }
}

void LevelSet::initialiseNarrowBand()
{
    unsigned int mineWidth = bandWidth - 1;

    // Reset the number of nodes in the narrow band.
    nNarrowBand = 0;

    // Reset the number of mines.
    nMines = 0;

    // Loop over all nodes.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Flag node as inactive.
        mesh.nodes[i].isActive = false;

        /* Check that the node isn't in a masked region. If it's not, then check
           whether it lies on the domain boundary, and if it does then check that
           the boundary isn't fixed.
         */
        if (!mesh.nodes[i].isFixed && (!mesh.nodes[i].isDomain || !isFixed))
        {
            // Absolute value of the signed distance function.
            double absoluteSignedDistance = std::abs(signedDistance[i]);

            // Node lies inside band.
            if (absoluteSignedDistance < bandWidth)
            {
                // Flag node as active.
                mesh.nodes[i].isActive = true;

                // Update narrow band array.
                narrowBand[nNarrowBand] = i;

                // Increment number of nodes.
                nNarrowBand++;

                // Node lines at edge of band.
                if (absoluteSignedDistance > mineWidth)
                {
                    // Node is a mine.
                    mesh.nodes[i].isMine = true;

                    // Update mine array.
                    mines[nMines] = i;

                    // Increment mine count.
                    nMines++;

                    // TODO:
                    // Check when number of mines exceeds array size!
                }
            }
        }
    }
}

void LevelSet::initialiseVelocities(const std::vector<BoundaryPoint>& boundaryPoints)
{
    // Map boundary point velocities to nodes of the level set domain
    // using inverse squared distance interpolation.

    // Whether the velocity at a node has been set.
    bool isSet[mesh.nNodes];

    // Weighting factor for each node.
    double weight[mesh.nNodes];

    // Initialise arrays.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        isSet[i] = false;
        weight[i] = 0;
        velocity[i] = 0;
    }

    // Loop over all boundary points.
    for (unsigned int i=0;i<boundaryPoints.size();i++)
    {
        // Find the closest node.
        unsigned int node = mesh.getClosestNode(boundaryPoints[i].coord);

        // Distance from the boundary point to the node.
        double dx = mesh.nodes[node].coord.x - boundaryPoints[i].coord.x;
        double dy = mesh.nodes[node].coord.y - boundaryPoints[i].coord.y;

        // Squared distance.
        double rSqd = dx*dx + dy*dy;

        // If boundary point lies exactly on the node, then set velocity
        // to that of the boundary point.
        if (rSqd < 1e-6)
        {
            velocity[node] = boundaryPoints[i].velocity;
            weight[node] = 1.0;
            isSet[node] = true;
        }
        else
        {
            // Update velocity estimate if not already set.
            if (!isSet[node])
            {
                velocity[node] += boundaryPoints[i].velocity / rSqd;
                weight[node] += 1.0 / rSqd;
            }
        }

        // Loop over all neighbours of the node.
        for (unsigned int j=0;j<4;j++)
        {
            // Index of the neighbouring node.
            unsigned int neighbour = mesh.nodes[node].neighbours[j];

            // Make sure neighbour is in bounds.
            if (neighbour < mesh.nNodes)
            {
                // Distance from the boundary point to the node.
                double dx = mesh.nodes[neighbour].coord.x - boundaryPoints[i].coord.x;
                double dy = mesh.nodes[neighbour].coord.y - boundaryPoints[i].coord.y;

                // Squared distance.
                double rSqd = dx*dx + dy*dy;

                // If boundary point lies exactly on the node, then set velocity
                // to that of the boundary point.
                if (rSqd < 1e-6)
                {
                    velocity[neighbour] = boundaryPoints[i].velocity;
                    weight[neighbour] = 1.0;
                    isSet[neighbour] = true;
                }
                else if (rSqd <= 1.0)
                {
                    // Update velocity estimate if not already set.
                    if (!isSet[neighbour])
                    {
                        velocity[neighbour] += boundaryPoints[i].velocity / rSqd;
                        weight[neighbour] += 1.0 / rSqd;
                    }
                }
            }
        }
    }

    // Compute interpolated velocity.
    for (unsigned int i=0;i<nNarrowBand;i++)
    {
        unsigned int node = narrowBand[i];
        if (velocity[node]) velocity[node] /= weight[node];
    }
}

double LevelSet::computeGradient(const unsigned int node)
{
    // Nodal coordinates.
    unsigned int x = mesh.nodes[node].coord.x;
    unsigned int y = mesh.nodes[node].coord.y;

    // Nodal signed distance.
    double lsf = signedDistance[node];

    // Whether gradient has been computed.
    bool isGradient = false;

    // Zero the gradient.
    double grad = 0;

    // Node is on the left edge.
    if (x == 0)
    {
        // Node is at bottom left corner.
        if (y == 0)
        {
            // If signed distance at nodes to right and above is the same, then use
            // the diagonal node for computing the gradient.
            if ((std::abs(signedDistance[mesh.xyToIndex[x+1][y]] - lsf) < 1e-6) &&
                (std::abs(signedDistance[mesh.xyToIndex[x][y+1]] - lsf) < 1e-6))
            {
                // Calculate signed distance to diagonal node.
                grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x+1][y+1]]);
                grad *= sqrt(2.0);
                isGradient = true;
            }
        }

        // Node is at top left corner.
        else if (y == mesh.height)
        {
            // If signed distance at nodes to right and below is the same, then use
            // the diagonal node for computing the gradient.
            if ((std::abs(signedDistance[mesh.xyToIndex[x+1][y]] - lsf) < 1e-6) &&
                (std::abs(signedDistance[mesh.xyToIndex[x][y-1]] - lsf) < 1e-6))
            {
                // Calculate signed distance to diagonal node.
                grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x+1][y-1]]);
                grad *= sqrt(2.0);
                isGradient = true;
            }
        }
    }

    // Node is on the right edge.
    else if (x == mesh.width)
    {
        // Node is at bottom right corner.
        if (y == 0)
        {
            // If signed distance at nodes to left and above is the same, then use
            // the diagonal node for computing the gradient.
            if ((std::abs(signedDistance[mesh.xyToIndex[x-1][y]] - lsf) < 1e-6) &&
                (std::abs(signedDistance[mesh.xyToIndex[x][y+1]] - lsf) < 1e-6))
            {
                // Calculate signed distance to diagonal node.
                grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x-1][y+1]]);
                grad *= sqrt(2.0);
                isGradient = true;
            }
        }

        // Node is at top right corner.
        else if (y == mesh.height)
        {
            // If signed distance at nodes to left and below is the same, then use
            // the diagonal node for computing the gradient.
            if ((std::abs(signedDistance[mesh.xyToIndex[x-1][y]] - lsf) < 1e-6) &&
                (std::abs(signedDistance[mesh.xyToIndex[x][y-1]] - lsf) < 1e-6))
            {
                // Calculate signed distance to diagonal node.
                grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x-1][y-1]]);
                grad *= sqrt(2.0);
                isGradient = true;
            }
        }
    }

    // Gradient hasn't already been calculated.
    if (!isGradient)
    {
        // Stencil values for the WENO approximation.
        double v1, v2, v3, v4, v5;

        // Upwind direction.
        int sign = velocity[node] < 0 ? -1 : 1;

        // Derivatives to right.

        // Node on left-hand edge.
        if (x == 0)
        {
            v1 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
            v2 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
            v3 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

            // Approximate derivatives outside of domain.
            v4 = v3;
            v5 = v3;
        }

        // One node to right of left-hand edge.
        else if (x == 1)
        {
            v1 = signedDistance[mesh.xyToIndex[4][y]] - signedDistance[mesh.xyToIndex[3][y]];
            v2 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
            v3 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
            v4 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

            // Approximate derivatives outside of domain.
            v5 = v4;
        }

        // Node on right-hand edge.
        else if (x == mesh.width)
        {
            v5 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x-1][y]];

            // Approximate derivatives outside of domain.
            v3 = v4;
            v2 = v4;
            v1 = v4;
        }

        // One node to left of right-hand edge.
        else if (x == (mesh.width - 1))
        {
            v5 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x-1][y]];
            v3 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];

            // Approximate derivatives outside of domain.
            v2 = v3;
            v1 = v3;
        }

        // Two nodes to left of right-hand edge.
        else if (x == (mesh.width - 2))
        {
            v5 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x-1][y]];
            v3 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];
            v2 = signedDistance[mesh.xyToIndex[x+2][y]] - signedDistance[mesh.xyToIndex[x+1][y]];

            // Approximate derivatives outside of domain.
            v1 = v2;
        }

        // Node lies in bulk.
        else
        {
            v1 = signedDistance[mesh.xyToIndex[x+3][y]] - signedDistance[mesh.xyToIndex[x+2][y]];
            v2 = signedDistance[mesh.xyToIndex[x+2][y]] - signedDistance[mesh.xyToIndex[x+1][y]];
            v3 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x-1][y]];
            v5 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
        }

        double gradRight = sign * gradHJWENO(v1, v2, v3, v4, v5);

        // Derivatives to left.

        // Node on right-hand edge.
        if (x == mesh.width)
        {
            v1 = signedDistance[mesh.xyToIndex[x-2][y]] - signedDistance[mesh.xyToIndex[x-3][y]];
            v2 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
            v3 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x-1][y]];

            // Approximate derivatives outside of domain.
            v4 = v3;
            v5 = v3;
        }

        // One node to left of right-hand edge.
        else if (x == (mesh.width-1))
        {
            v1 = signedDistance[mesh.xyToIndex[x-2][y]] - signedDistance[mesh.xyToIndex[x-3][y]];
            v2 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
            v3 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x-1][y]];
            v4 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];

            // Approximate derivatives outside of domain.
            v5 = v4;
        }

        // Node on left-hand edge.
        else if (x == 0)
        {
            v5 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
            v4 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

            // Approximate derivatives outside of domain.
            v3 = v4;
            v2 = v4;
            v1 = v4;
        }

        // One node to right of left-hand edge.
        else if (x == 1)
        {
            v5 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
            v4 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
            v3 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

            // Approximate derivatives outside of domain.
            v2 = v3;
            v1 = v3;
        }

        // Two nodes to right of left-hand edge.
        else if (x == 2)
        {
            v5 = signedDistance[mesh.xyToIndex[4][y]] - signedDistance[mesh.xyToIndex[3][y]];
            v4 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
            v3 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
            v2 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

            // Approximate derivatives outside of domain.
            v1 = v2;
        }

        // Node lies in bulk.
        else
        {
            v1 = signedDistance[mesh.xyToIndex[x-2][y]] - signedDistance[mesh.xyToIndex[x-3][y]];
            v2 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
            v3 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x-1][y]];
            v4 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];
            v5 = signedDistance[mesh.xyToIndex[x+2][y]] - signedDistance[mesh.xyToIndex[x+1][y]];
        }

        double gradLeft = sign * gradHJWENO(v1, v2, v3, v4, v5);

        // Upward derivatives.

        // Node on bottom edge.
        if (y == 0)
        {
            v1 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
            v2 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
            v3 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

            // Approximate derivatives outside of domain.
            v4 = v3;
            v5 = v3;
        }

        // One node above bottom edge.
        else if (y == 1)
        {
            v1 = signedDistance[mesh.xyToIndex[x][4]] - signedDistance[mesh.xyToIndex[x][3]];
            v2 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
            v3 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
            v4 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

            // Approximate derivatives outside of domain.
            v5 = v4;
        }

        // Node is on top edge.
        else if (y == mesh.height)
        {
            v5 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x][y-1]];

            // Approximate derivatives outside of domain.
            v3 = v4;
            v2 = v4;
            v1 = v4;
        }

        // One node below top edge.
        else if (y == (mesh.height - 1))
        {
            v5 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x][y-1]];
            v3 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];

            // Approximate derivatives outside of domain.
            v2 = v3;
            v1 = v3;
        }

        // Two nodes below top edge.
        else if (y == (mesh.height - 2))
        {
            v5 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x][y-1]];
            v3 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];
            v2 = signedDistance[mesh.xyToIndex[x][y+2]] - signedDistance[mesh.xyToIndex[x][y+1]];

            // Approximate derivatives outside of domain.
            v1 = v2;
        }

        // Node lies in bulk.
        else
        {
            v1 = signedDistance[mesh.xyToIndex[x][y+3]] - signedDistance[mesh.xyToIndex[x][y+2]];
            v2 = signedDistance[mesh.xyToIndex[x][y+2]] - signedDistance[mesh.xyToIndex[x][y+1]];
            v3 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];
            v4 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x][y-1]];
            v5 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
        }

        double gradUp = sign * gradHJWENO(v1, v2, v3, v4, v5);

        // Downward derivative.

        // Node on top edge.
        if (y == mesh.height)
        {
            v1 = signedDistance[mesh.xyToIndex[x][y-2]] - signedDistance[mesh.xyToIndex[x][y-3]];
            v2 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
            v3 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x][y-1]];

            // Approximate derivatives outside of domain.
            v4 = v3;
            v5 = v3;
        }

        // One node below top edge.
        else if (y == (mesh.height - 1))
        {
            v1 = signedDistance[mesh.xyToIndex[x][y-2]] - signedDistance[mesh.xyToIndex[x][y-3]];
            v2 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
            v3 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x][y-1]];
            v4 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];

            // Approximate derivatives outside of domain.
            v5 = v4;
        }

        // Node lies on bottom edge
        else if (y == 0)
        {
            v5 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
            v4 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

            // Approximate derivatives outside of domain.
            v3 = v4;
            v2 = v4;
            v1 = v4;
        }

        // One node above bottom edge.
        else if (y == 1)
        {
            v5 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
            v4 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
            v3 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

            // Approximate derivatives outside of domain.
            v2 = v3;
            v1 = v3;
        }

        // Two nodes above bottom edge.
        else if (y == 2)
        {
            v5 = signedDistance[mesh.xyToIndex[x][4]] - signedDistance[mesh.xyToIndex[x][3]];
            v4 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
            v3 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
            v2 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

            // Approximate derivatives outside of domain.
            v1 = v2;
        }

        // Node lies in bulk.
        else
        {
            v1 = signedDistance[mesh.xyToIndex[x][y-2]] - signedDistance[mesh.xyToIndex[x][y-3]];
            v2 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
            v3 = signedDistance[mesh.xyToIndex[x][y]]   - signedDistance[mesh.xyToIndex[x][y-1]];
            v4 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];
            v5 = signedDistance[mesh.xyToIndex[x][y+2]] - signedDistance[mesh.xyToIndex[x][y+1]];
        }

        double gradDown = sign * gradHJWENO(v1, v2, v3, v4, v5);

        // Compute gradient using upwind scheme.

        if (gradDown > 0)   grad += gradDown * gradDown;
        if (gradLeft > 0)   grad += gradLeft * gradLeft;
        if (gradUp < 0)     grad += gradUp * gradUp;
        if (gradRight < 0)  grad += gradRight * gradRight;

        grad = sqrt(grad);
    }

    // Return gradient.
    return grad;
}

double LevelSet::gradHJWENO(double v1, double v2, double v3, double v4, double v5)
{
    // Approximate the gradient using the 5th order Hamilton-Jacobi WENO approximation.
    // Taken from pages 34-35 of "Level Set Methods and Dynamic Implicit Surfaces".
    // See: http://web.stanford.edu/class/cs237c/Lecture16.pdf

    double oneQuarter        = 1.0  / 4.0;
    double thirteenTwelths   = 13.0 / 12.0;
    double eps               = 1e-6;

    // Estimate the smoothness of each stencil.

    double s1 = thirteenTwelths * (v1 - 2*v2 + v3)*(v1 - 2*v2 + v3)
              + oneQuarter * (v1 - 4*v2 + 3*v3)*(v1 - 4*v2 + 3*v3);

    double s2 = thirteenTwelths * (v2 - 2*v3 + v4)*(v2 - 2*v3 + v4)
              + oneQuarter * (v2 - v4)*(v2 - v4);

    double s3 = thirteenTwelths * (v3 - 2*v4 + v5)*(v3 - 2*v4 + v5)
              + oneQuarter * (3*v3 - 4*v4 + v5)*(3*v3 - 4*v4 + v5);

    // Compute the alpha values for each stencil.

    double alpha1 = 0.1 / ((s1 + eps)*(s1 + eps));
    double alpha2 = 0.6 / ((s2 + eps)*(s2 + eps));
    double alpha3 = 0.3 / ((s3 + eps)*(s3 + eps));

    // Calculate the normalised weights.

    double totalWeight = alpha1 + alpha2 + alpha3;

    double w1 = alpha1 / totalWeight;
    double w2 = alpha2 / totalWeight;
    double w3 = alpha3 / totalWeight;

    // Sum the three stencil components.
    double grad = w1 * (2*v1 - 7*v2 + 11*v3)
                + w2 * (5*v3 - v2 + 2*v4)
                + w3 * (2*v3 + 5*v4 - v5);

    grad *= (1.0 / 6.0);

    return grad;
}

double LevelSet::pointToLineDistance(const Coord& vertex1, const Coord& vertex2, const Coord& point) const
{
    // Separation components between vertices.
    double dx = vertex2.x - vertex1.x;
    double dy = vertex2.y - vertex1.y;

    // Squared separation.
    double rSqd = dx*dx + dy*dy;

    // Return separation between point and either vertex.
    if (rSqd < 1e-6)
    {
        dx = point.x - vertex1.x;
        dy = point.y - vertex1.y;

        return sqrt(dx*dx + dy*dy);
    }

    /* Consider the line extending the segment, parameterized as v + t (w - v).
       We find the projection of the point p onto the line. It falls where

         t = [(p-v) . (w-v)] / |w-v|^2

       We clamp t from [0,1] to handle points outside the segment vw, i.e. if
       t < 0 then we use the distance from p to v, and if t > 1 we use the
       distance from p to w.

       (Where v = vertex1, w = vertex2, and p = point)
     */
    else
    {
        double t = ((point.x - vertex1.x) * dx + (point.y - vertex1.y) * dy ) / rSqd;

        t = std::max(0.0, std::min(1.0, t));

        // Project onto the line.
        double x = vertex1.x + t * dx;
        double y = vertex1.y + t * dy;

        // Compute distance from the point.
        dx = x - point.x;
        dy = y - point.y;

        return sqrt(dx*dx + dy*dy);
    }
}

bool LevelSet::isInsidePolygon(const Coord& point, const std::vector<Coord>& vertices) const
{
    /* Test whether a point lies inside a polygon.

       The polygon is defined by a vector of vertices. These must be ordered
       and closed, i.e. there are n vertices with vertices[n] = vertices[0].

       The test calculates the winding number of the point and returns false
       only when it is zero.

       Adapted from: http://geomalgorithms.com/a03-_inclusion.html
     */

    int windingNumber = 0;

    // Loop through all vertices.
    for (unsigned int i=0;i<vertices.size()-1;i++)
    {
        // Point lies above the vertex.
        if (vertices[i].y <= point.y)
        {
            // Point lies below the next vertex.
            if (vertices[i+1].y > point.y)

                // Point is to the left of the edge.
                if (isLeftOfLine(vertices[i], vertices[i+1], point) > 0)
                    windingNumber++;
        }
        else
        {
            // Point lies above the vertex.
            if (vertices[i+1].y <= point.y)

                // Point lies to the right of the edge.
                if (isLeftOfLine(vertices[i], vertices[i+1], point) < 0)
                    windingNumber--;
        }
    }

    if (windingNumber == 0) return false;
    else return true;
}

int LevelSet::isLeftOfLine(const Coord& vertex1, const Coord& vertex2, const Coord& point) const
{
    return ((vertex2.x - vertex1.x) * (point.y - vertex1.y)
        - (point.x -  vertex1.x) * (vertex2.y - vertex1.y));
}