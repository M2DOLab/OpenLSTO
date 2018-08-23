/*! \file InputOutput.cpp
    \brief A class for reading and writing data.
 */

InputOutput::InputOutput() {}

void InputOutput::saveLevelSetVTK(const unsigned int& datapoint, const LevelSet& levelSet,
    bool isVelocity, bool isGradient, const std::string& outputDirectory) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "level-set_" << num.str() << ".vtk";

    saveLevelSetVTK(fileName, levelSet);
}

void InputOutput::saveLevelSetVTK(const std::ostringstream& fileName,
    const LevelSet& levelSet, bool isVelocity, bool isGradient) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        lsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Set up ParaView header information.
    fprintf(pFile, "# vtk DataFile Version 3.0\n");
    fprintf(pFile, "Para0\n");
    fprintf(pFile, "ASCII\n");
    fprintf(pFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(pFile, "DIMENSIONS %d %d %d\n", 1 + levelSet.mesh.width, 1 + levelSet.mesh.height, 1);
    fprintf(pFile, "X_COORDINATES %d int\n", 1 + levelSet.mesh.width);
    for (unsigned int i=0;i<=levelSet.mesh.width;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nY_COORDINATES %d int\n", 1 + levelSet.mesh.height);
    for (unsigned int i=0;i<=levelSet.mesh.height;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nZ_COORDINATES 1 int\n0\n\n");
    fprintf(pFile, "POINT_DATA %d\n", levelSet.mesh.nNodes);

    // Write the nodal signed distance to file.
    fprintf(pFile, "SCALARS distance float 1\n");
    fprintf(pFile, "LOOKUP_TABLE default\n");
    for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
        fprintf(pFile, "%lf\n", levelSet.signedDistance[i]);

    // Write the nodal velocity to file.
    if (isVelocity)
    {
        fprintf(pFile, "SCALARS velocity float 1\n");
        fprintf(pFile, "LOOKUP_TABLE default\n");
        for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
            fprintf(pFile, "%lf\n", levelSet.velocity[i]);
    }

    // Write the nodal gradient to file.
    if (isGradient)
    {
        fprintf(pFile, "SCALARS gradient float 1\n");
        fprintf(pFile, "LOOKUP_TABLE default\n");
        for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
            fprintf(pFile, "%lf\n", levelSet.gradient[i]);
    }

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveLevelSetTXT(const unsigned int& datapoint,
    const LevelSet& levelSet, const std::string& outputDirectory, bool isXY) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "level-set_" << num.str() << ".txt";

    saveLevelSetTXT(fileName, levelSet, isXY);
}

void InputOutput::saveLevelSetTXT(const std::ostringstream& fileName,
    const LevelSet& levelSet, bool isXY) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        lsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Write the nodal signed distance to file.
    for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
    {
        if (isXY) fprintf(pFile, "%lf %lf ", levelSet.mesh.nodes[i].coord.x, levelSet.mesh.nodes[i].coord.y);
        fprintf(pFile, "%lf %lf %lf\n", levelSet.signedDistance[i], levelSet.velocity[i], levelSet.gradient[i]);
    }

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveBoundaryPointsTXT(const unsigned int& datapoint,
    const Boundary& boundary, const std::string& outputDirectory) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "boundary-points_" << num.str() << ".txt";

    saveBoundaryPointsTXT(fileName, boundary);
}

void InputOutput::saveBoundaryPointsTXT(const std::ostringstream& fileName, const Boundary& boundary) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        lsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Write the boundary points to file.
    for (unsigned int i=0;i<boundary.nPoints;i++)
        fprintf(pFile, "%lf %lf %lf\n",
            boundary.points[i].coord.x, boundary.points[i].coord.y, boundary.points[i].length);

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveBoundarySegmentsTXT(const unsigned int& datapoint,
    const Boundary& boundary, const std::string& outputDirectory) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "boundary-segments_" << num.str() << ".txt";

    saveBoundarySegmentsTXT(fileName, boundary);
}

void InputOutput::saveBoundarySegmentsTXT(const std::ostringstream& fileName, const Boundary& boundary) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        lsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Write the boundary points to file.
    for (unsigned int i=0;i<boundary.nSegments;i++)
    {
        // Start and end points of the segment.
        unsigned int start = boundary.segments[i].start;
        unsigned int end = boundary.segments[i].end;

        // Coordinates.
        double x, y;

        // First point.
        x = boundary.points[start].coord.x;
        y = boundary.points[start].coord.y;

        // Write boundary point to file.
        fprintf(pFile, "%lf %lf\n", x, y);

        // Second point.
        x = boundary.points[end].coord.x;
        y = boundary.points[end].coord.y;

        // Write boundary point to file.
        fprintf(pFile, "%lf %lf\n\n", x, y);
    }

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveAreaFractionsVTK(const unsigned int& datapoint,
    const Mesh& mesh, const std::string& outputDirectory) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "area_" << num.str() << ".vtk";

    saveAreaFractionsVTK(fileName, mesh);
}

void InputOutput::saveAreaFractionsVTK(const std::ostringstream& fileName, const Mesh& mesh) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        lsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Set up ParaView header information.
    fprintf(pFile, "# vtk DataFile Version 3.0\n");
    fprintf(pFile, "Para0\n");
    fprintf(pFile, "ASCII\n");
    fprintf(pFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(pFile, "DIMENSIONS %d %d %d\n", 1 + mesh.width, 1 + mesh.height, 1);
    fprintf(pFile, "X_COORDINATES %d int\n", 1 + mesh.width);
    for (unsigned int i=0;i<=mesh.width;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nY_COORDINATES %d int\n", 1 + mesh.height);
    for (unsigned int i=0;i<=mesh.height;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nZ_COORDINATES 1 int\n0\n\n");

    // Write the element area fractions to file.
    fprintf(pFile, "CELL_DATA %d\n", mesh.nElements);
    fprintf(pFile, "SCALARS area float 1\n");
    fprintf(pFile, "LOOKUP_TABLE default\n");
    for (unsigned int i=0;i<mesh.nElements;i++)
        fprintf(pFile, "%lf\n", mesh.elements[i].area);

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveAreaFractionsTXT(const unsigned int& datapoint,
    const Mesh& mesh, const std::string& outputDirectory, bool isXY) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "area_" << num.str() << ".txt";

    saveAreaFractionsTXT(fileName, mesh, isXY);
}

void InputOutput::saveAreaFractionsTXT(const std::ostringstream& fileName, const Mesh& mesh, bool isXY) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        lsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Write the element area fractions to file.
    for (unsigned int i=0;i<mesh.nElements;i++)
    {
        if (isXY) fprintf(pFile, "%lf %lf ", mesh.elements[i].coord.x, mesh.elements[i].coord.y);
        fprintf(pFile, "%lf\n", mesh.elements[i].area);
    }

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

bool InputOutput::BoundaryVTK (const std::ostringstream& FieldLabel, std::vector < std::vector <double> >& BoundaryPoints, std::vector < std::vector <double> >& Sensitivities, std::vector < std::vector <unsigned int> >& Neighbors)
{
    //std::stringstream FileName;
    //FileName = *FieldLabel + ".vtk";
    std::ofstream File;
    File.open (FieldLabel.str().c_str());
    if (!File)
    {
        // std::cout << "NULL File pointer while printing node data into " << FileName << "\n";
        return false;
    }
    unsigned int nPoints = BoundaryPoints.size ();
    unsigned int nSegments = Neighbors.size ();
    File << "# vtk DataFile Version 3.0\n";
    File << "Para0\n";
    File << "ASCII\n";
    File << "DATASET UNSTRUCTURED_GRID\n";
    File << "POINTS\t" << nPoints << "\tdouble\n";
    unsigned int SpaceDimension = BoundaryPoints [0].size ();
    for (unsigned int i = 0; i < nPoints; i = i + 1)
    {
        for (unsigned int j = 0; j < SpaceDimension; j = j + 1)
        {
            File << BoundaryPoints [i] [j] << "\t";
        }
        File << "0\n";
    }
    File << "CELLS\t" << nSegments << "\t" << nSegments * 3 << "\n";
    for (unsigned int i = 0; i < nSegments; i = i + 1)
        File << 2 << "\t" << Neighbors [i] [0] << "\t" << Neighbors [i] [1] << "\n";
    File << "CELL_TYPES\t" << nSegments << "\n";
    for (unsigned int i = 0; i < nSegments; i = i + 1)
        File << "3\n";
    File << "POINT_DATA\t" << nPoints << "\n";
    for (unsigned int i = 0; i < Sensitivities.size (); i = i + 1)
    {
        File << "SCALARS\tSensitivity" << std::to_string (i + 1) << "\tdouble\t1\n";
        File << "LOOKUP_TABLE DEFAULT\n";
        for (unsigned int j = 0; j < nPoints; j = j + 1)
        {
            File << Sensitivities [i] [j] << "\n";
        }
        File << "\n";
    }
    File.close ();
    return true;
}

// Functionality to write sensitivities for the current load case stored.
void InputOutput::WriteOptimisationHistoryTXT(std::vector<double> objectives, std::vector<std::vector<double> > constraints) {

    // Defining file name variables.
    std::ostringstream fileName, num;

    // Creating file name.
    fileName.str("");
    fileName << "Output/optimisation_history.txt";

    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    // Total number of iterations (including iteration 0).
    int n_it = objectives.size();

    // Total number of constraints.
    int n_cons = constraints.size();

    for (int i = 0; i < n_it; i++)
    {
        // Printing objective and first constraint values.
        fprintf(pFile, "%lf \t %lf", objectives[i], constraints[0][i]);

        // Printing other constraints
        for (int j = 1; j < n_cons; j++)
        {
            // Printing each stress value.
            fprintf(pFile, "%lf \t", constraints[j][i]);
        }
        fprintf(pFile, "\n");
    }

    fclose(pFile);

    return;
}

