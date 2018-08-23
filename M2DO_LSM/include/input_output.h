#ifndef _INPUTOUTPUT_H
#define _INPUTOUTPUT_H

/*! \file InputOutput.h
    \brief A class for reading and writing data.
 */

//! A class for reading and writing data.
class InputOutput
{
public:
    //! Constructor.
    InputOutput();

    //! Save the level set function as a ParaView VTK file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param levelSet
            A reference to the level set object.

        \param isVelocity
            Whether to write velocity information to file (optional).

        \param isGradient
            Whether to write gradient information to file (optional).

        \param outputDirectory
            The output directory path (optional).
     */
    void saveLevelSetVTK(const unsigned int&, const LevelSet&, bool isVelocity = false,
        bool isGradient = false, const std::string& outputDirectory = "") const;

    //! Save the level set function as a ParaView VTK file.
    /*! \param fileName
            The name of the data file.

        \param levelSet
            A reference to the level set object.

        \param isVelocity
            Whether to write velocity information to file (optional).

        \param isGradient
            Whether to write gradient information to file (optional).
     */
    void saveLevelSetVTK(const std::ostringstream&, const LevelSet&,
        bool isVelocity = false, bool isGradient = false) const;

    //! Save the level set function as a plain text file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param levelSet
            A reference to the level set object.

        \param outputDirectory
            The output directory path (optional).

        \param isXY
            Whether to also output the nodal x/y coordinates (optional).
     */
    void saveLevelSetTXT(const unsigned int&, const LevelSet&,
        const std::string& outputDirectory = "", bool isXY = false) const;

    //! Save the level set function as a plain text file.
    /*! \param fileName
            The name of the data file.

        \param levelSet
            A reference to the level set object.

        \param isXY
            Whether to also output the nodal x/y coordinates (optional).
     */
    void saveLevelSetTXT(const std::ostringstream&, const LevelSet&, bool isXY = false) const;

    //! Save boundary points as a plain text file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param boundary
            A reference to the boundary object.

        \param outputDirectory
            The output directory path (optional).
     */
    void saveBoundaryPointsTXT(const unsigned int&, const Boundary&, const std::string& outputDirectory = "") const;

    //! Save boundary points as a plain text file.
    /*! \param fileName
            The name of the data file.

        \param boundary
            A reference to the boundary object.
     */
    void saveBoundaryPointsTXT(const std::ostringstream&, const Boundary&) const;

    //! Save boundary segments as a plain text file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param boundary
            A reference to the boundary object.

        \param outputDirectory
            The output directory path (optional).
     */
    void saveBoundarySegmentsTXT(const unsigned int&,
        const Boundary&, const std::string& outputDirectory = "") const;

    //! Save boundary segments as a plain text file.
    /*! \param fileName
            The name of the data file.

        \param boundary
            A reference to the boundary object.
     */
    void saveBoundarySegmentsTXT(const std::ostringstream&, const Boundary&) const;

    //! Save the element area fractions as a ParaView VTK file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param mesh
            A reference to the level set mesh.

        \param outputDirectory
            The output directory path (optional).
     */
    void saveAreaFractionsVTK(const unsigned int&, const Mesh&, const std::string& outputDirectory = "") const;

    //! Save the element area fractions as a ParaView VTK file.
    /*! \param fileName
            The name of the data file.

        \param mesh
            A reference to the level set mesh.
     */
    void saveAreaFractionsVTK(const std::ostringstream&, const Mesh&) const;

    //! Save element area fractions as a plain text file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param mesh
            A reference to the level set mesh.

        \param outputDirectory
            The output directory path (optional).

        \param isXY
            Whether to also output the element x/y coordinates (optional).
     */
    void saveAreaFractionsTXT(const unsigned int&, const Mesh&,
        const std::string& outputDirectory = "", bool isXY = false) const;

    //! Save element area fractions as a plain text file.
    /*! \param fileName
            The name of the data file.

        \param mesh
            A reference to the level set mesh.

        \param isXY
            Whether to also output the element x/y coordinates (optional).
     */
    void saveAreaFractionsTXT(const std::ostringstream&, const Mesh&, bool isXY = false) const;

    // Plot boundary sensitivities
    bool BoundaryVTK (const std::ostringstream& FieldLabel, std::vector < std::vector <double> >& BoundaryPoints, std::vector < std::vector <double> >& Sensitivities, std::vector < std::vector <unsigned int> >& Neighbors);

    // Write optimisation history in .txt file.
    void WriteOptimisationHistoryTXT(std::vector<double>, std::vector<std::vector<double> >);
};

#include "../src/input_output.cpp"

#endif  /* _INPUTOUTPUT_H */
