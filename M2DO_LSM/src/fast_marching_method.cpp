// Adapted from Scikit-FMM: https://github.com/scikit-fmm/scikit-fmm

/*! \file FastMarchingMethod.cpp
    \brief An implementation of the Fast Marching Method.
 */

FastMarchingMethod::FastMarchingMethod(const Mesh& mesh_, bool isTest_) :
    mesh(mesh_),
    isTest(isTest_),
    outOfBounds(mesh.nNodes)
{
    heap = nullptr;

    // Resize data structures.
    heapPtr.resize(mesh.nNodes);
    nodeStatus.resize(mesh.nNodes);
    signedDistanceCopy.resize(mesh.nNodes);
}

FastMarchingMethod::~FastMarchingMethod()
{
    delete heap;
}

void FastMarchingMethod::march(std::vector<double>& signedDistance_)
{
    signedDistance = &signedDistance_;
    isVelocity = false;

    // Initialise the set of frozen boundary nodes.
    initialiseFrozen();

    // Initialise the heap data structure.
    initialiseHeap();

    // Initialise the set of trial nodes adjacent to the boundary.
    initialiseTrial();

    // Find the fast marching solution.
    solve();
}

void FastMarchingMethod::march(std::vector<double>& signedDistance_, std::vector<double>& velocity_)
{
    /* Extend boundary velocities to all nodes within the narrow band region.

       Note that this method assumes that boundary point velocities have
       already been mapped to the level set nodes using inverse squared
       distance interpolation, or similar.
     */

    signedDistance = &signedDistance_;
    velocity = &velocity_;
    isVelocity = true;

    // Initialise the set of frozen boundary nodes.
    initialiseFrozen();

    // Initialise the heap data structure.
    initialiseHeap();

    // Initialise the set of trial nodes adjacent to the boundary.
    initialiseTrial();

    // Find the fast marching solution.
    solve();

    // Restore the original signed distance function. Only update velocities.
    (*signedDistance) = signedDistanceCopy;
}

void FastMarchingMethod::initialiseFrozen()
{
    // The number of frozen nodes.
    unsigned int nFrozen = 0;

    // First find all zero values of the level set.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Store a copy of the level set.
        signedDistanceCopy[i] = (*signedDistance)[i];

        // Make sure node isn't masked.
        if (nodeStatus[i] != FMM_NodeStatus::MASKED)
        {
            // Zero contour passes through node.
            if (signedDistanceCopy[i] == 0)
            {
                // Mark node as frozen.
                nodeStatus[i] = FMM_NodeStatus::FROZEN;

                // Increment number of frozen nodes.
                nFrozen++;
            }
        }
    }

    // Now check whether the neighbours of each node (in any direction)
    // are on opposite sides of the zero contour.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Whether level set changes sign between a node and its neighbour.
        bool isBorder = false;

        // Only consider nodes that haven't yet been frozen.
        if (nodeStatus[i] == FMM_NodeStatus::NONE)
        {
            // Initialise distance array.
            double dist[2] = {0, 0};

            // Loop over all neighbours.
            for (unsigned int j=0;j<4;j++)
            {
                // Neighbours are ordered: left, right, down, up.

                // Get index of neighbour.
                unsigned int neighbour = mesh.nodes[i].neighbours[j];

                // Make sure neighbour lies inside domain boundary.
                if (neighbour != outOfBounds)
                {
                    // Level set changes sign along direction.
                    if ((signedDistanceCopy[i] * signedDistanceCopy[neighbour]) < 0)
                    {
                        isBorder = true;

                        // Calculate the distance to the zero contour (linear interpolation).
                        double d = signedDistanceCopy[i] / (signedDistanceCopy[i] - signedDistanceCopy[neighbour]);

                        // Set dimension (neighbours 0 and 1 are x dimension, 2 and 3 are y).
                        unsigned int dim = (j < 2) ? 0 : 1;

                        // Check if distance is less than current value.
                        if (dist[dim] == 0 || dist[dim] > d)
                            dist[dim] = d;
                    }
                }
            }

            // Node and neighbour span the zero contour.
            if (isBorder)
            {
                double distSum = 0;

                // Calculate perpendicular distance to boundary (Pythag.)

                for (unsigned int j=0;j<2;j++)
                {
                    if (dist[j] > 0)
                        distSum += 1.0 / (dist[j] * dist[j]);
                }

                // Update signed distance.
                if (signedDistanceCopy[i] < 0) (*signedDistance)[i] = -sqrt(1.0 / distSum);
                else (*signedDistance)[i] = sqrt(1.0 / distSum);

                // Flag node as frozen.
                nodeStatus[i] = FMM_NodeStatus::FROZEN;

                // Increment number of frozen nodes.
                nFrozen++;
            }
        }
    }

    errno = 0;
    lsm_check(nFrozen > 0, "There are no frozen nodes!");

    return;

error:
    exit(EXIT_FAILURE);
}

void FastMarchingMethod::initialiseHeap()
{
    // Initialise the maximum heap size.
    unsigned int maxHeapSize = 0;

    // Loop over all nodes.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Find nodes without a status (far field).
        if (nodeStatus[i] == FMM_NodeStatus::NONE)
            maxHeapSize++;
    }

    heap = new Heap(maxHeapSize, isTest);
}

void FastMarchingMethod::initialiseTrial()
{
    // For each node, check whether it has a frozen neighbour.
    // If so, calculate the distance from the zero contour and insert into heap.

    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Node hasn't yet been given a status.
        if (nodeStatus[i] == FMM_NodeStatus::NONE)
        {
            // Loop over all nearest neighbour nodes.
            for (unsigned int j=0;j<4;j++)
            {
                unsigned int neighbour = mesh.nodes[i].neighbours[j];

                // Neighbour lies within the domain boundary.
                if (neighbour != outOfBounds)
                {
                    // Neighbour is frozen.
                    if (nodeStatus[neighbour] & FMM_NodeStatus::FROZEN)
                    {
                        // Check that node status hasn't been updated.
                        if (nodeStatus[i] == FMM_NodeStatus::NONE)
                        {
                            if (isVelocity)
                            {
                                // Node lies inside the narrow band region.
                                if (mesh.nodes[i].isActive)
                                {
                                    // Flag node as in trial band.
                                    nodeStatus[i] = FMM_NodeStatus::TRIAL;

                                    // Get distance from zero contour.
                                    (*signedDistance)[i] = updateNode(i);

                                    // Add to heap.
                                    heapPtr[i] = heap->push(i, std::abs((*signedDistance)[i]));
                                }
                            }
                            else
                            {
                                // Flag node as in trial band.
                                nodeStatus[i] = FMM_NodeStatus::TRIAL;

                                // Get distance from zero contour.
                                (*signedDistance)[i] = updateNode(i);

                                // Add to heap.
                                heapPtr[i] = heap->push(i, std::abs((*signedDistance)[i]));
                            }
                        }
                    }
                }
            }
        }
    }
}

void FastMarchingMethod::solve()
{
    /* This is the fast marching method main loop. The order
       of operations is as follows...

        (1) Choose from the trial set the node with the smallest
            absolute distance from the zero contour and freeze it.

        (2) For each neighbour of the frozen node, calculate a new
            distance based on frozen elements.

             - Mark each neighbour as a trial node and insert into
               the heap.

             - If the neighbour is already in the heap, update the
               distance value.
     */

    // Number of nodes to freeze.
    unsigned int nFrozen;

    // Indices of nodes that are to be frozen.
    std::vector<unsigned int> toFreeze(mesh.nNodes);

    while (!heap->empty())
    {
        unsigned int addr;
        double value;

        // Zero number of frozen nodes.
        nFrozen = 0;

        // Pop top entry off heap.
        heap->pop(addr, value);

        // Mark node as frozen.
        nodeStatus[addr] = FMM_NodeStatus::FROZEN;

        // Set final velocity.
        if (isVelocity) finaliseVelocity(addr);

        // Increment number of frozen nodes.
        toFreeze[nFrozen] = addr;
        nFrozen++;

        // Inititialise loop termination flag.
        bool isDone = false;

        while (!isDone)
        {
            if (!heap->empty() && (value == heap->peek()))
            {
                unsigned int l_addr;
                double l_value;

                // Pop top entry off heap.
                heap->pop(l_addr, l_value);

                // Mark node as frozen.
                nodeStatus[l_addr] = FMM_NodeStatus::FROZEN;

                // Set final velocity.
                if (isVelocity) finaliseVelocity(l_addr);

                // Increment number of frozen nodes.
                toFreeze[nFrozen] = l_addr;
                nFrozen++;
            }
            else isDone = true;
        }

        // Loop over all frozen nodes.
        for (unsigned int i=0;i<nFrozen;i++)
        {
            // Get node address.
            unsigned int addr = toFreeze[i];

            // Loop over all neighbours of frozen node.
            for (unsigned int j=0;j<4;j++)
            {
                // Get address of neighbour.
                unsigned int naddr = mesh.nodes[addr].neighbours[j];

                // Neighbour lies within domain boundary.
                if (naddr != outOfBounds)
                {
                    // Neighbour hasn't been frozen.
                    if (nodeStatus[naddr] != FMM_NodeStatus::FROZEN)
                    {
                        // Calculate and store udpdated distance estimate.
                        double d = updateNode(naddr);
                        (*signedDistance)[naddr] = d;

                        // Neighbour is in trial band.
                        if (nodeStatus[naddr] & FMM_NodeStatus::TRIAL)
                        {
                            // Update value in heap.
                            heap->set(heapPtr[naddr], std::abs(d));
                        }
                        // Neighbour has no status (far field).
                        else if (nodeStatus[naddr] == FMM_NodeStatus::NONE)
                        {
                            if (isVelocity)
                            {
                                // Node lies inside the narrow band region.
                                if (mesh.nodes[naddr].isActive)
                                {
                                    // Mark node as in trial band.
                                    nodeStatus[naddr] = FMM_NodeStatus::TRIAL;

                                    // Push onto heap.
                                    heapPtr[naddr] = heap->push(naddr, std::abs(d));
                                }
                            }
                            else
                            {
                                // Mark node as in trial band.
                                nodeStatus[naddr] = FMM_NodeStatus::TRIAL;

                                // Push onto heap.
                                heapPtr[naddr] = heap->push(naddr, std::abs(d));
                            }
                        }

                        // Now update the far field point in the second order stencil.
                        // "jump" over a frozen node if needed.

                        // Address of second nearest neighbour in the same direction.
                        naddr = mesh.nodes[naddr].neighbours[j];

                        // Next nearest neighbour lies within domain boundary.
                        if (naddr != outOfBounds)
                        {
                            // Neighbour is in the trial band.
                            if (nodeStatus[naddr] & FMM_NodeStatus::TRIAL)
                            {
                                // Calculate and store udpdated distance estimate.
                                double d = updateNode(naddr);
                                (*signedDistance)[naddr] = d;

                                // Update value in heap.
                                heap->set(heapPtr[naddr], std::abs(d));
                            }
                        }
                    }
                }
            }
        }
    }
}

double FastMarchingMethod::updateNode(unsigned int node)
{
    // Reused constants.
    const double aa = 9.0/4.0;
    const double oneThird = 1.0/3.0;

    // Quadratic coefficients.
    double a, b, c;

    // Zero coefficients.
    a = b = c = 0;

    double dist1, dist2;

    // Loop over all dimensions.
    for (unsigned int i=0;i<2;i++)
    {
        // Initialise distances.
        dist1 = maxDouble;
        dist2 = maxDouble;

        // Loop over all directions.
        for (unsigned int j=0;j<2;j++)
        {
            // Work out index of neighbour.
            unsigned int index = 2*i + j;

            // First neighbour.
            unsigned int n1 = mesh.nodes[node].neighbours[index];

            // Neighbour is within the domain boundary.
            if (n1 != outOfBounds)
            {
                // Neighbour is frozen.
                if (nodeStatus[n1] & FMM_NodeStatus::FROZEN)
                {
                    // Make sure neighbour is closer to the zero contour (upwind).
                    if (std::abs((*signedDistance)[n1]) < std::abs(dist1))
                    {
                        // Store distance.
                        dist1 = (*signedDistance)[n1];

                        // Second neighbour in same direction.
                        unsigned int n2 = mesh.nodes[n1].neighbours[index];

                        // Neighbour is within the domain boundary.
                        if (n2 != outOfBounds)
                        {
                            // Neighbour is frozen.
                            if (nodeStatus[n2] & FMM_NodeStatus::FROZEN)
                            {
                                // Make sure neighbour is closer to the zero contour (upwind).
                                if (std::abs((*signedDistance)[n2]) <= std::abs(dist1))
                                {
                                    // Store distance.
                                    dist2 = (*signedDistance)[n2];
                                }
                            }
                        }
                    }
                }
            }
        }

        // Second order finite difference.
        if (dist2 < maxDouble)
        {
            double tp = oneThird*(4*dist1 - dist2);

            a += aa;
            b -= 2*aa*tp;
            c += aa*tp*tp;
        }
        // First order finite difference.
        else if (dist1 < maxDouble)
        {
            a += 1;
            b -= 2*dist1;
            c += dist1*dist1;
        }
    }

    // Update coefficient.
    c -= 1;

    return solveQuadratic(node, a, b, c);
}

void FastMarchingMethod::finaliseVelocity(unsigned int node)
{
    // Set the velocity of this node, i.e.
    // find v_ext, where grad v_ext . grad phi = 0

    // Technically, we don't need to calculate this extension
    // velocity until the point is frozen.

    // Initialise distance array.
    double dist[2] = {0, 0};

    // Initialise front (zero contour) distance array.
    double frontDist[2] = {0, 0};

    // Initialise velocity array.
    double vel[2] = {0, 0};

    // Loop over all neighbours of the node.
    for (unsigned int i=0;i<4;i++)
    {
        // Set dimension (neighbours 0 and 1 are x dimension, 2 and 3 are y).
        unsigned int dim = (i < 2) ? 0 : 1;

        // Get index of neighbour.
        unsigned int neighbour = mesh.nodes[node].neighbours[i];

        // Neighbour is within domain boundary.
        if (neighbour != outOfBounds)
        {
            // Neighbour is frozen.
            if (nodeStatus[neighbour] & FMM_NodeStatus::FROZEN)
            {
                // Absolute signed distance of the neighbouring node.
                double d = std::abs((*signedDistance)[neighbour]);

                // Check whether the neighbour is closer to the zero contour.
                if ((frontDist[dim] == 0) || (frontDist[dim] > d))
                {
                    // Store updated distance to the front.
                    frontDist[dim] = d;

                    // Calculate the distance to the front in this direction.
                    d = (*signedDistance)[node] - (*signedDistance)[neighbour];

                    // Store absolute distance and velocity.
                    dist[dim] = std::abs(d);
                    vel[dim] = (*velocity)[neighbour];
                }
            }
        }
    }

    double numerator = 0;
    double denominator = 0;

    for (unsigned int i=0;i<2;i++)
    {
        numerator += dist[i] * vel[i];
        denominator += dist[i];
    }

    errno = 0;
    lsm_check(denominator != 0, "Divide by zero error.");

    (*velocity)[node] = numerator / denominator;

    return;

error:
    exit(EXIT_FAILURE);
}

double FastMarchingMethod::solveQuadratic(unsigned int node, const double& a, const double& b, const double& c) const
{
    // Initialise roots.
    double r0, r1;

    // Calculate discriminant.
    double discrim = b*b - 4*a*c;

    // Solve for roots.
    if (discrim > 0)
    {
        r0 = (-b + sqrt(discrim)) / (2.0 * a);
        r1 = (-b - sqrt(discrim)) / (2.0 * a);
    }
    // Use previous estimate.
    else return (*signedDistance)[node];

    if (signedDistanceCopy[node] > doubleEpsilon) return r0;
    else return r1;
}
