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
  
/*! \file Boundary.cpp
    \brief A class for the discretised boundary.
 */

#ifndef _BOUNDARY_HOLE_H
#define _BOUNDARY_HOLE_H
    

namespace M2DO_LSM {

    class Boundary_hole
    {
    public:
        /// Vector of boundary points.
        std::vector<BoundaryPoint> points;

        /// Vector of boundary segments.
        std::vector<BoundarySegment> segments;

        /// The number of boundary points.
        unsigned int nPoints;

        /// The number of boundary segments.
        unsigned int nSegments;

        /// The total length of the boundary.
        double length;

        /// The total area fraction of the mesh.
        double area;

        LevelSet levelSet;
        std::vector<double>* signedDistance;


        //! Constructor.
        /*! \param levelSet_
                A reference to the level set object.
         */

        Boundary_hole(LevelSet levelSet_, std::vector<double>* signedDistance_) : levelSet(levelSet_), signedDistance(signedDistance_){}
        //Boundary_hole(LevelSet levelSet){}

        //! Use linear interpolation to compute the discretised boundary
        /*! \param isTarget
                Whether to discretise the target signed distance function.
         */
        double computeAreaFractions(bool isTarget = false) {

            // Clear and reserve vector memory.
            points.clear();
            segments.clear();
            points.reserve(levelSet.mesh.nNodes);
            segments.reserve(levelSet.mesh.nNodes);

            // Get level set mesh elements and nodes;
            std::vector<M2DO_LSM::Element> elements_temp(levelSet.mesh.nElements);
            for (int i = 0; i<levelSet.mesh.nElements; i++) {
                elements_temp[i] = levelSet.mesh.elements[i];
            }

            std::vector<M2DO_LSM::Node> nodes_temp(levelSet.mesh.nNodes);
            for (int i = 0; i<levelSet.mesh.nNodes; i++) {
                nodes_temp[i] = levelSet.mesh.nodes[i];
            }

            // std::vector<M2DO_LSM::Element>* elements_temp;
            // std::vector<M2DO_LSM::Node>* nodes_temp;
            // elements_temp = &levelSet.mesh.elements;
            // nodes_temp = &levelSet.mesh.nodes;

            // Reset the number of points and segments.
            nPoints = nSegments = 0;

            // Zero the boundary length.
            length = 0;

            // // Initialise a pointer to the signed distance vector.
            // std::vector<double>* signedDistance;

            // // Point to the target signed distance.
            // if (isTarget) signedDistance = &levelSet.target;

            // // Point to the current signed distance function.
            // else signedDistance = &levelSet.signedDistance;

            // Compute the status of nodes and elements in level-set mesh.
            computeMeshStatus(signedDistance, elements_temp, nodes_temp);

            // Calculate the material area fraction in each element.
            //computeAreaFractions(elements_temp, nodes_temp);
            
        //}

        //! Calculate the material area fraction in each element.
        /*! \return
                The total element area fraction.
         */
        //double computeAreaFractions(vector<M2DO_LSM::Element>& elements_temp, vector<M2DO_LSM::Node>& nodes_temp) {
            // Zero the total area fraction.
            area = 0;

            for (unsigned int i=0;i<levelSet.mesh.nElements;i++)
            {
                // Element is inside structure.
                if (elements_temp[i].status & ElementStatus::INSIDE)
                    elements_temp[i].area = 1.0;

                // Element is outside structure.
                else if (elements_temp[i].status & ElementStatus::OUTSIDE)
                    elements_temp[i].area = 0.0;

                // Element is cut by the boundary.
                else elements_temp[i].area = cutArea(elements_temp[i],nodes_temp);

                // Add the area to the running total.
                area += elements_temp[i].area;
            }

            return area;
        }

        //! Determine the status of the elements and nodes of the level set mesh.
        /*! \param signedDistance
                A pointer to the signed distance function vector.
         */
        void computeMeshStatus(const std::vector<double>* signedDistance, vector<M2DO_LSM::Element>& elements_temp, vector<M2DO_LSM::Node>& nodes_temp) const {
            // Calculate node status.
            for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
            {
                // Reset the number of boundary points associated with the element.
                nodes_temp[i].nBoundaryPoints = 0;

                // Flag node as being on the boundary if the signed distance is within
                // a small tolerance of the zero contour. This avoids problems with
                // rounding errors when generating the discretised boundary.
                if (std::abs((*signedDistance)[i]) < 1e-6)
                {
                    nodes_temp[i].status = NodeStatus::BOUNDARY;
                }
                else if ((*signedDistance)[i] < 0)
                {
                    nodes_temp[i].status = NodeStatus::OUTSIDE;
                }
                else nodes_temp[i].status = NodeStatus::INSIDE;
            }

            // Calculate element status.
            for (unsigned int i=0;i<levelSet.mesh.nElements;i++)
            {
                // Tally counters for the element's node statistics.
                unsigned int tallyInside = 0;
                unsigned int tallyOutside = 0;

                // Reset the number of boundary segments associated with the element.
                elements_temp[i].nBoundarySegments = 0;

                // Loop over each node of the element.
                for (unsigned int j=0;j<4;j++)
                {
                    unsigned int node = elements_temp[i].nodes[j];

                    if (nodes_temp[node].status & NodeStatus::INSIDE) tallyInside++;
                    else if (nodes_temp[node].status & NodeStatus::OUTSIDE) tallyOutside++;
                }

                // No nodes are outside: element is inside the structure.
                if (tallyOutside == 0) elements_temp[i].status = ElementStatus::INSIDE;

                // No nodes are inside: element is outside the structure.
                else if (tallyInside == 0) elements_temp[i].status = ElementStatus::OUTSIDE;

                // Otherwise no status.
                else elements_temp[i].status = ElementStatus::NONE;
            }
        }

        //! Calculate the material area for an element cut by the boundary.
        /*! \param element
                A reference to the element.

            \return
                The area fraction.
         */
        double cutArea(const Element& element, vector<M2DO_LSM::Node>& nodes_temp) {
            // Number of polygon vertices.
            unsigned int nVertices = 0;

            // Polygon vertices (maximum of six).
            std::vector<Coord> vertices(6);

            // Whether we're searching for nodes that are inside or outside the boundary.
            NodeStatus::NodeStatus status;

            if (element.status & ElementStatus::CENTRE_OUTSIDE) status = NodeStatus::OUTSIDE;
            else status = NodeStatus::INSIDE;

            // Check all nodes of the element.
            for (unsigned int i=0;i<4;i++)
            {
                // Node index;
                unsigned int node = element.nodes[i];

                // Node matches status.
                if (nodes_temp[node].status & status)
                {
                    // Add coordinates to vertex array.
                    vertices[nVertices].x = nodes_temp[node].coord.x;
                    vertices[nVertices].y = nodes_temp[node].coord.y;

                    // Increment number of vertices.
                    nVertices++;
                }

                // Node is on the boundary.
                else if (nodes_temp[node].status & NodeStatus::BOUNDARY)
                {
                    // Next node.
                    unsigned int n1 = (i == 3) ? 0 : (i + 1);
                    n1 = element.nodes[n1];

                    // Previous node.
                    unsigned int n2 = (i == 0) ? 3 : (i - 1);
                    n2 = element.nodes[n2];

                    // Check that node isn't part of a boundary segment, i.e. both of its
                    // neighbours are inside the structure.
                    if ((nodes_temp[n1].status & NodeStatus::INSIDE) &&
                        (nodes_temp[n2].status & NodeStatus::INSIDE))
                    {
                        // Add coordinates to vertex array.
                        vertices[nVertices].x = nodes_temp[node].coord.x;
                        vertices[nVertices].y = nodes_temp[node].coord.y;

                        // Increment number of vertices.
                        nVertices++;
                    }
                }
            }

            // Add boundary segment start and end points.
            for (unsigned int i=0;i<element.nBoundarySegments;i++)
            {
                // Segment index.
                unsigned int segment = element.boundarySegments[i];

                // Add start point coordinates to vertices array.
                vertices[nVertices].x = points[segments[segment].start].coord.x;
                vertices[nVertices].y = points[segments[segment].start].coord.y;

                // Increment number of vertices.
                nVertices++;

                // Add end point coordinates to vertices array.
                vertices[nVertices].x = points[segments[segment].end].coord.x;
                vertices[nVertices].y = points[segments[segment].end].coord.y;

                // Increment number of vertices.
                nVertices++;
            }

            // Return area of the polygon.
            if (element.status & ElementStatus::CENTRE_OUTSIDE)
                return (1.0 - polygonArea(vertices, nVertices, element.coord));
            else
                return polygonArea(vertices, nVertices, element.coord);
        }

        //! Whether a point is clockwise of another. The origin point is 12 o'clock.
        /*! \param point1
                The coordinates of the first point.

            \param point2
                The coordinates of the second point.

            \param centre
                The coordinates of the element centre.

            \return
                Whether the first point is clockwise of the second.
         */
        bool isClockwise(const Coord& point1, const Coord& point2, const Coord& centre) const {
            if ((point1.x - centre.x) >= 0 && (point2.x - centre.x) < 0)
                return false;

            if ((point1.x - centre.x) < 0 && (point2.x - centre.x) >= 0)
                return true;

            if ((point1.x - centre.x) == 0 && (point2.x - centre.x) == 0)
            {
                if ((point1.y - centre.y) >= 0 || (point2.y - centre.y) >= 0)
                    return (point1.y > point2.y) ? false : true;

                return (point2.y > point1.y) ? false : true;
            }

            // Compute the cross product of the vectors (centre --> point1) x (centre --> point2).
            double det = (point1.x - centre.x) * (point2.y - centre.y)
                       - (point2.x - centre.x) * (point1.y - centre.y);

            if (det < 0) return false;
            else return true;

            // Points are on the same line from the centre, check which point is
            // closer to the centre.

            double d1 = (point1.x - centre.x) * (point1.x - centre.x)
                      + (point1.y - centre.y) * (point1.y - centre.y);

            double d2 = (point2.x - centre.x) * (point2.x - centre.x)
                      + (point2.y - centre.y) * (point2.y - centre.y);

            return (d1 > d2) ? false : true;
        }

        //! Return the area of a polygon.
        /*! \param vertices
                A clockwise ordered vector of polygon vertices.

            \param nVertices
                The number of vertices.

            \param centre
                The coordinates of the element centre.

            \return
                The area of the polygon.
         */
        double polygonArea(std::vector<Coord>& vertices, const unsigned int& nVertices, const Coord& centre) const {
            double area = 0;

            // Sort vertices in anticlockwise order.
            std::sort(vertices.begin(), vertices.begin() + nVertices, std::bind(&Boundary_hole::isClockwise,
                this, std::placeholders::_1, std::placeholders::_2, centre));

            // Loop over all vertices.
            for (unsigned int i=0;i<nVertices;i++)
            {
                // Next point around (looping back to beginning).
                unsigned int j = (i == (nVertices - 1)) ? 0 : (i + 1);

                area += vertices[i].x * vertices[j].y;
                area -= vertices[j].x * vertices[i].y;
            }

            area *= 0.5;

            return (std::abs(area));
        }
    };

}

#endif



