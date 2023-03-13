#include <SFML/Window/Keyboard.hpp>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <chrono>
#include "cartesian_grid_of_speed.hpp"
#include "vortex.hpp"
#include "cloud_of_points.hpp"
#include "runge_kutta.hpp"
#include "screen.hpp"

#include <mpi.h>

auto readConfigFile(std::ifstream &input)
{
    using point = Simulation::Vortices::point;

    int isMobile;
    std::size_t nbVortices;
    Numeric::CartesianGridOfSpeed cartesianGrid;
    Geometry::CloudOfPoints cloudOfPoints;
    constexpr std::size_t maxBuffer = 8192;
    char buffer[maxBuffer];
    std::string sbuffer;
    std::stringstream ibuffer;


    // Lit la première ligne de commentaire :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lecture de la grille cartésienne
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    double xleft, ybot, h;
    std::size_t nx, ny;
    ibuffer >> xleft >> ybot >> nx >> ny >> h;
    cartesianGrid = Numeric::CartesianGridOfSpeed({nx, ny}, point{xleft, ybot}, h);


    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit mode de génération des particules
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    int modeGeneration;
    ibuffer >> modeGeneration;


    if (modeGeneration == 0) // Génération sur toute la grille
    {
        std::size_t nbPoints;
        ibuffer >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {cartesianGrid.getLeftBottomVertex(), cartesianGrid.getRightTopVertex()});
    }
    else
    {
        std::size_t nbPoints;
        double xl, xr, yb, yt;
        ibuffer >> xl >> yb >> xr >> yt >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {point{xl, yb}, point{xr, yt}});
    }


    // Lit le nombre de vortex :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit le nombre de vortex
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    try
    {
        ibuffer >> nbVortices;
    }
    catch (std::ios_base::failure &err)
    {
        std::cout << "Error " << err.what() << " found" << std::endl;
        std::cout << "Read line : " << sbuffer << std::endl;
        throw err;
    }
    Simulation::Vortices vortices(nbVortices, {cartesianGrid.getLeftBottomVertex(),
                                               cartesianGrid.getRightTopVertex()});


    input.getline(buffer, maxBuffer); // Relit un commentaire
    for (std::size_t iVortex = 0; iVortex < nbVortices; ++iVortex)
    {
        input.getline(buffer, maxBuffer);
        double x, y, force;
        std::string sbuffer(buffer, maxBuffer);
        std::stringstream ibuffer(sbuffer);
        ibuffer >> x >> y >> force;
        vortices.setVortex(iVortex, point{x, y}, force);
    }
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit le mode de déplacement des vortex
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    ibuffer >> isMobile;
    return std::make_tuple(vortices, isMobile, cartesianGrid, cloudOfPoints);
}

int main(int nargs, char *argv[])
{
    MPI_Init(&nargs, &argv);

    int nProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

    int process;
    MPI_Comm_rank(MPI_COMM_WORLD, &process);

    double dt;
    // if (process == 0)
    // {

    std::size_t resx = 800, resy = 800; // square window
    if (nargs > 3)
    {
        resx = std::stoull(argv[2]);
        resy = std::stoull(argv[3]);
    }

    auto vortices = std::get<0>(config);
    auto isMobile = std::get<1>(config);
    auto grid = std::get<2>(config);
    auto cloud = std::get<3>(config);

    std::cout << "######## Vortex simultor ########" << std::endl
              << std::endl;
    std::cout << "Press P for play animation " << std::endl;
    std::cout << "Press S to stop animation" << std::endl;
    std::cout << "Press right cursor to advance step by step in time" << std::endl;
    std::cout << "Press down cursor to halve the time step" << std::endl;
    std::cout << "Press up cursor to double the time step" << std::endl;

    grid.updateVelocityField(vortices);

    Graphisme::Screen myScreen({resx, resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()});
    bool animate = false;
    double dt = 0.1;

            if (process == 0)
            {


                        MPI_Send(&dt, 1, MPI_DOUBLE, 1, 101, MPI_COMM_WORLD);
                        // MPI_Send(&dt, 1, MPI_DOUBLE, 1, 101, MPI_COMM_WORLD);
    MPI_Finalize();
            if (process == 1)
            {
                MPI_Recv(&dt, 1, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Recv(&dt, 1, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD, &status);
                std::cout << "print dt: " << dt << std::endl;
            }
    MPI_Finalize();

    return EXIT_SUCCESS;
}