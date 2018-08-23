#ifndef _LEVELSET_H
#define _LEVELSET_H

/*! \file LevelSet.h
    \brief A class for the level set function.
 */

/*! \brief A class for the level set function.

    The level set is represented as a signed distance function from the zero
    countour. Positive values are inside the structure, negative values are
    outside.

    The class provides methods for initialising and reinitialising the signed
    distance function. The default initialisation uses a set of linearly spaced
    circular holes to create the classic "Swiss cheese" starting configuration.
    Alternatively, the user can pass the constructor a vector of holes or point
    coordinates that will be used to define the initial structure. A vector
    of coordinates defining the interface of a target structure can be used
    for shape matching simulations.

    Reinitialisation of the signed distance function is performed using an
    implementation of the Fast Marching Method. A second-order stencil is used
    for the upwind finite difference scheme where possible.

    Functionality is also provided for tracking nodes that are part of the
    narrow band region around the zero contour, as well as mine nodes at
    the edge of the narrow band.
 */
class LevelSet
{
public:
    //! Constructor.
    /*! \param mesh_
            A reference to the fixed-grid mesh.

        \param moveLimit_
            The CFL limit (in units of the mesh grid spacing).

        \param bandWidth_
            The width of the narrow band region.

        \param isFixed_
            Whether the domain boundary is fixed.
     */
    LevelSet(Mesh&, double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixed_ = false);

    //! Constructor.
    /*! \param mesh_
            A reference to the fixed-grid mesh.

        \param holes
            A vector of holes.

        \param moveLimit_
            The CFL limit (in units of the mesh grid spacing).

        \param bandWidth_
            The width of the narrow band region.

        \param isFixed_
            Whether the domain boundary is fixed.
     */
    LevelSet(Mesh&, const std::vector<Hole>&, double moveLimit_ = 0.5,
        unsigned int bandWidth_ = 6, bool isFixed_ = false);

    //! Constructor.
    /*! \param mesh_
            A reference to the fixed-grid mesh.

        \param points
            A vector of point coordinates (clockwise ordered and closed).

        \param moveLimit_
            The CFL limit (in units of the mesh grid spacing).

        \param bandWidth_
            The width of the narrow band region.

        \param isFixed_
            Whether the domain boundary is fixed.
     */
    LevelSet(Mesh&, const std::vector<Coord>&, double moveLimit_ = 0.5,
        unsigned int bandWidth_ = 6, bool isFixed_ = false);

    //! Constructor.
    /*! \param mesh_
            A reference to the fixed-grid mesh.

        \param initialHoles
            A vector of holes for the initial interface.

        \param targetHoles
            A vector of holes for the target interface.

        \param moveLimit_
            The CFL limit (in units of the mesh grid spacing).

        \param bandWidth_
            The width of the narrow band region.

        \param isFixed_
            Whether the domain boundary is fixed.
     */
    LevelSet(Mesh&, const std::vector<Hole>&, const std::vector<Hole>&,
        double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixed_ = false);

    //! Constructor.
    /*! \param mesh_
            A reference to the fixed-grid mesh.

        \param holes
            A vector of holes for the initial interface.

        \param points
            A vector of point coordinates for the target interface (clockwise ordered and closed).

        \param moveLimit_
            The CFL limit (in units of the mesh grid spacing).

        \param bandWidth_
            The width of the narrow band region.

        \param isFixed_
            Whether the domain boundary is fixed.
     */
    LevelSet(Mesh&, const std::vector<Hole>&, const std::vector<Coord>&,
        double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixed_ = false);

    //! Constructor.
    /*! \param mesh_
            A reference to the fixed-grid mesh.

        \param initialPoints
            A vector of point coordinates for the initial interface (clockwise ordered and closed).

        \param targetPoints
            A vector of point coordinates for the target interface (clockwise ordered and closed).

        \param moveLimit_
            The CFL limit (in units of the mesh grid spacing).

        \param bandWidth_
            The width of the narrow band region.

        \param isFixed_
            Whether the domain boundary is fixed.
     */
    LevelSet(Mesh&, const std::vector<Coord>&, const std::vector<Coord>&,
        double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixed_ = false);

    //! Update the level set function.
    /*! \param timeStep
            The time step.

        \return
            Whether the signed distance was reinitialised.
     */
    bool update(double);
    
    //! Kill nodes in a region of the domain.
    /*! param points
            A reference to a vector of points.
     */
    void killNodes(const std::vector<Coord>&);

    //! Fix nodes in a region of the domain.
    /*! param points
            A reference to a vector of points.
     */
    void fixNodes(const std::vector<Coord>&);

    //! Create Level Set boundary in a region of the domain.
    /*! param points
            A reference to a vector of points.
     */
    void createLevelSetBoundary(const std::vector<Coord>&);

    //! Reinitialise the level set to a signed distance function.
    void reinitialise();

    //! Extend boundary point velocities to the level set nodes.
    /*! \param boundaryPoints
            A reference to a vector of boundary points.
     */
    void computeVelocities(const std::vector<BoundaryPoint>&);

    //! Extend boundary point velocities to the level set nodes.
    /*! \param boundaryPoints
            A reference to a vector of boundary points.

        \param timeStep
            The time step for the level set update.

        \param temperature
            The temperature of the thermal bath.

        \param rng
            A reference to the random number generator.

        \return
            The time step scaling factor.
     */
    double computeVelocities(std::vector<BoundaryPoint>&, double&, const double, MersenneTwister&);

    //! Compute the modulus of the gradient of the signed distance function.
    void computeGradients();

    std::vector<double> signedDistance;     //!< The nodal signed distance function (level set).
    std::vector<double> velocity;           //!< The nodal signed normal velocity.
    std::vector<double> gradient;           //!< The nodal gradient of the level set function (modulus).
    std::vector<double> target;             //!< Signed distance target (for shape matching).
    std::vector<unsigned int> narrowBand;   //!< Indices of nodes in the narrow band.
    std::vector<unsigned int> mines;        //!< Indices of nodes at the edge of the narrow band.
    unsigned int nNarrowBand;               //!< The number of nodes in narrow band.
    unsigned int nMines;                    //!< The number of mine nodes.
    const double moveLimit;                 //!< The boundary movement limit (CFL condition).

    Mesh& mesh;                             //!< A reference to the fixed-grid mesh.

private:
    unsigned int bandWidth;                 //!< The width of the narrow band region.
    bool isFixed;                           //!< Whether the domain boundary is fixed.
    bool isTarget;                          //!< Whether a target interface is defined.

    //! Default initialisation of the level set function (Swiss cheese configuration).
    void initialise();

    //! Initialise the level set from a vector of user-defined holes.
    /* \param holes
            A vector of holes.
     */
    void initialise(const std::vector<Hole>&);

    //! Initialise the level set from a vector of user-defined points.
    /* \param points
            A vector of point coordinates.
     */
    void initialise(const std::vector<Coord>&);

    //! Helper function for initialise methods.
    //! Initialises the level set function as the distance to the closest domain boundary.
    void closestDomainBoundary();

    //! Initialise the narrow band region.
    void initialiseNarrowBand();

    //! Initialise velocities for boundary nodes.
    /*! \param boundaryPoints
            A reference to a vector of boundary points.
     */
    void initialiseVelocities(const std::vector<BoundaryPoint>&);

    //! Compute the modulus of the gradient of the signed distance function at a node.
    /*! \param node
            The node index.

        \return
            The gradient at the node.
     */
    double computeGradient(const unsigned int);

    //! Compute Hamilton-Jacobi WENO gradient approximation.
    /*! \param v1
            The value of the function at the first stencil point.

        \param v2
            The value of the function at the second stencil point.

        \param v3
            The value of the function at the third stencil point.

        \param v4
            The value of the function at the fourth stencil point.

        \param v5
            The value of the function at the fifth stencil point.

        \return
            The smoothed function (gradient).
     */
    double gradHJWENO(double, double, double, double, double);

    //! Compute the minimum distance between a point and a line segment.
    /*! \param vertex1
            The coordinate of the first vertex.

        \param vertex2
            The coordinate of the second vertex.

        \param point
            The coordinates of the point of interest.

        \return
            The minimum distance to the line segment.
     */
    double pointToLineDistance(const Coord&, const Coord&, const Coord&) const;

    //! Test if a point lies inside a polygon.
    /*! \param point
            The coordinates of the point.

        \param vertices
            The vertices of the polygon (closed and ordered).

        \return
            Whether the point lies inside the polygon.
     */
    bool isInsidePolygon(const Coord&, const std::vector<Coord>&) const;

    //! Test if a point lies left, on, or right of an infinite line.
    /*! \param vertex1
            The coordinate of the first vertex.

        \param vertex2
            The coordinate of the second vertex.

        \param point
            The coordinates of the point of interest.

        \return
            > 0 if the point is left of the line through vertex1 and vertex2,
            = 0 if the point is on the line,
            < 0 if right of the line.
     */
    int isLeftOfLine(const Coord&, const Coord&, const Coord&) const;
};

#include "../src/level_set.cpp"

#endif  /* _LEVELSET_H */
