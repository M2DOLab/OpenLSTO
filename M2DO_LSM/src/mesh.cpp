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

/*! \file Mesh.cpp
    \brief A class for the level-set domain fixed-grid mesh.
 */

Mesh::Mesh(unsigned int width_,
           unsigned int height_,
           bool isPeriodic_) :

           width(width_),
           height(height_),
           nElements(width*height),
           nNodes((1+width)*(1+height)),
           isPeriodic(isPeriodic_)
{
    // Resize element and node data structures.
    elements.resize(nElements);
    nodes.resize(nNodes);

    // Resize 2D to 1D mapping vector.
    xyToIndex.resize(width+1);
    for (unsigned int i=0;i<width+1;i++)
        xyToIndex[i].resize(height+1);

    // Calculate node nearest neighbours.
    initialiseNodes();

    // Initialise elements (and node to element connectivity).
    initialiseElements();
}

unsigned int Mesh::getClosestNode(const Coord& point) const
{
    // Get element index.
    unsigned int element = getElement(point);

    // Work out distance relative to element centre.
    double dx = point.x - elements[element].coord.x;
    double dy = point.y - elements[element].coord.y;

    // Point lies in left half.
    if (dx < 0)
    {
        // Lower left quadrant.
        if (dy < 0) return elements[element].nodes[0];

        // Upper left quadrant.
        else return elements[element].nodes[3];
    }

    // Point lies in right half.
    else
    {
        // Lower right quadrant.
        if (dy < 0) return elements[element].nodes[1];

        // Upper right quadrant.
        else return elements[element].nodes[2];
    }
}

unsigned int Mesh::getClosestNode(double x, double y) const
{
    // Get element index.
    unsigned int element = getElement(x, y);

    // Work out distance relative to element centre.
    double dx = x - elements[element].coord.x;
    double dy = y - elements[element].coord.y;

    // Point lies in left half.
    if (dx < 0)
    {
        // Lower left quadrant.
        if (dy < 0) return elements[element].nodes[0];

        // Upper left quadrant.
        else return elements[element].nodes[3];
    }

    // Point lies in right half.
    else
    {
        // Lower right quadrant.
        if (dy < 0) return elements[element].nodes[1];

        // Upper right quadrant.
        else return elements[element].nodes[2];
    }
}

unsigned int Mesh::getElement(const Coord& point) const
{
    // Subtract a small value to ensure that we round down.
    double x = point.x - 1e-6;
    double y = point.y - 1e-6;

    // Enforce lower bound.
    if (x < 0) x = 0;
    if (y < 0) y = 0;

    // Work out x and y element indices (cells are unit width).
    unsigned int elementX = std::floor(x);
    unsigned int elementY = std::floor(y);

    // Return global element index.
    return (elementY*width + elementX);
}

unsigned int Mesh::getElement(double x, double y) const
{
    // Subtract a small value to ensure that we round down.
    x -= 1e-6;
    y -= 1e-6;

    // Enforce lower bound.
    if (x < 0) x = 0;
    if (y < 0) y = 0;

    // Work out x and y element indices (cells are unit width).
    unsigned int elementX = std::floor(x);
    unsigned int elementY = std::floor(y);

    // Return global element index.
    return (elementY*width + elementX);
}

void Mesh::createMeshBoundary(const std::vector<Coord>& points)
{
    // Loop over all nodes.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Point is inside the rectangle.
        if (nodes[i].coord.x > points[0].x &&
            nodes[i].coord.y > points[0].y &&
            nodes[i].coord.x < points[1].x &&
            nodes[i].coord.y < points[1].y)
        {
            // Make mesh node a new domain boundary.
            nodes[i].isDomain = true;
        }
    }
}

void Mesh::initialiseNodes()
{
    // Coordinates of the node.
    unsigned int x, y;

    // Loop over all nodes.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Mark node as in bulk.
        nodes[i].isDomain = false;

        // Mark node as unmasked.
        nodes[i].isMasked = false;

        // Mark node as not fixed.
        nodes[i].isFixed = false;

        // Zero number of connected elements.
        nodes[i].nElements = 0;

        // Zero number of boundary points.
        nodes[i].nBoundaryPoints = 0;

        // Work out node coordinates.
        x = i % (width + 1);
        y = int(i / (width + 1));

        // Node lies on the domain boundary.
        if ((x == 0) || (x == width) || (y == 0) || (y == height))
            nodes[i].isDomain = true;

        // Set node coordinates.
        nodes[i].coord.x = x;
        nodes[i].coord.y = y;

        // Add to 2D mapping vector.
        xyToIndex[x][y] = i;

        // Determine nearest neighbours.
        initialiseNeighbours(i, x, y);
    }
}

void Mesh::initialiseElements()
{
    // Coordinates of the element.
    unsigned int x, y;

    // Number of nodes along width of mesh (number of elements plus one).
    unsigned int w = width + 1;

    // Loop over all elements.
    for (unsigned int i=0;i<nElements;i++)
    {
        // Work out element coordinates.
        x = i % width;
        y = int(i / width);

        // Store coordinates of elemente centre.
        elements[i].coord.x = x + 0.5;
        elements[i].coord.y = y + 0.5;

        // Store connectivity (element --> node)

        // Node on bottom left corner of element.
        elements[i].nodes[0] = x + (y * w);

        // Node on bottom right corner of element.
        elements[i].nodes[1] = x + 1 + (y * w);

        // Node on top right corner of element.
        elements[i].nodes[2] = x + 1 + ((y + 1) * w);

        // Node on top right corner of element.
        elements[i].nodes[3] = x + ((y + 1) * w);

        // Fill reverse connectivity arrays (node --> element)
        for (unsigned int j=0;j<4;j++)
        {
            unsigned int node = elements[i].nodes[j];
            nodes[node].elements[nodes[node].nElements] = i;
            nodes[node].nElements++;
        }
    }
}

void Mesh::initialiseNeighbours(unsigned int node, unsigned int x, unsigned int y)
{
    // Number of nodes along width and height of mesh (number of elements plus one).
    unsigned int w = width + 1;
    unsigned int h = height + 1;

    // Neighbours to left and right.
    nodes[node].neighbours[0] = (x - 1 + w) % w + (y * w);
    nodes[node].neighbours[1] = (x + 1 + w) % w + (y * w);

    // Neighbours below and above.
    nodes[node].neighbours[2] = x + (w * ((y - 1 + h) % h));
    nodes[node].neighbours[3] = x + (w * ((y + 1 + h) % h));

    // The mesh isn't periodic, flag out of bounds neighbours.
    if (!isPeriodic)
    {
        // Node is on first or last row.
        if (x == 0) nodes[node].neighbours[0] = nNodes;
        else if (x == width) nodes[node].neighbours[1] = nNodes;

        // Node is on first or last column.
        if (y == 0) nodes[node].neighbours[2] = nNodes;
        else if (y == height) nodes[node].neighbours[3] = nNodes;
    }
}

