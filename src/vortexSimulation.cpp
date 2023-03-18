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
    // initializing MPI
    int numProcesses, process;

    MPI_Init(&nargs, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses); //
    MPI_Comm_rank(MPI_COMM_WORLD, &process);      //

    MPI_Status status;   // for MPI_Iprobe
    MPI_Request request; // for MPI_Isend

    // both processes will create the following variables
    char const *filename;
    if (nargs == 1)
    {
        std::cout << "Usage : vortexsimulator <nom fichier configuration>" << std::endl;
        return EXIT_FAILURE;
    }

    // initial variables
    filename = argv[1];
    std::ifstream fich(filename);
    auto config = readConfigFile(fich);
    fich.close();

    std::size_t resx = 1080, resy = 1080;
    if (nargs > 3)
    {
        resx = std::stoull(argv[2]);
        resy = std::stoull(argv[3]);
    }

    // initializing variables from file
    auto vortices = std::get<0>(config);
    auto isMobile = std::get<1>(config);
    auto grid = std::get<2>(config);
    auto cloud = std::get<3>(config);

    // initialize grid with vortices
    grid.updateVelocityField(vortices);

    bool animate = true; // continue calculation?
    double dt = 0.1;     // time step

    int isReceiving = false; // flag, are there any message to receive?
    int isFinished = false;  // flag, window was closed?

    // program will be divided in 2 process
    //  0:  update and management of screen
    //  1:  execute calculus
    if (process == 0)
    {
        // user instructions
        Graphisme::Screen myScreen({resx, resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()});
        std::cout << "######## Vortex simultor ########" << std::endl
                  << std::endl;
        std::cout << "Press P for play animation " << std::endl;
        std::cout << "Press S to stop animation" << std::endl;
        std::cout << "Press right cursor to advance step by step in time" << std::endl;
        std::cout << "Press down cursor  to halve  the time step" << std::endl;
        std::cout << "Press up cursor    to double the time step" << std::endl;

        // initializing variable to hold keyboard's commands
        auto start = std::chrono::system_clock::now();

        // comunication between processes is time consuming
        // processes will only exchange a variable to reduce time lost with comunication
        char command = '_';
        int FPS = 0;
        int frameCount = 0;
        while (myScreen.isOpen())
        {
            // on inspecte tous les évènements de la fenêtre qui ont été émis depuis la précédente itération
            sf::Event event;
            // read screen for key pressed to evaluate
            while (myScreen.pollEvent(event) && command != 'E' && command != 'K')
            {
                if (event.type == sf::Event::Resized)
                    myScreen.resize(event); // event resize screen
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::P))
                    command = 'P'; // play animation
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
                    command = 'S'; // stop animation
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
                {
                    command = 'U'; // +speed animation
                    dt *= 2;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
                {
                    command = 'D'; // -speed animatin
                    dt /= 2;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
                    command = 'A'; // advance
                if (event.type == sf::Event::Closed || sf::Keyboard::isKeyPressed(sf::Keyboard::E))
                    command = 'E'; // close window and terminate other processes

                if (command != '_')
                    MPI_Send(&command, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
            }

            // scaning for message to receive
            MPI_Iprobe(1, MPI_ANY_TAG, MPI_COMM_WORLD, &isReceiving, &status);
            if (isReceiving)
            {
                // scaning for end process
                MPI_Iprobe(1, 0, MPI_COMM_WORLD, &isFinished, &status);
                if (isFinished)
                {
                    MPI_Recv(&command, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    if (command == 'K')
                    {
                        myScreen.close();
                        MPI_Finalize();
                        return 0;
                    }
                }

                // receiving data
                MPI_Recv(cloud.data(), 2 * cloud.numberOfPoints(), MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(grid.data(), 2 * grid.getSizeGrid(), MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(vortices.data(), 3 * vortices.numberOfVortices(), MPI_DOUBLE, 1, 3, MPI_COMM_WORLD, &status);
                frameCount++;

                // updating screen
                myScreen.clear(sf::Color::Black);

                std::string strDt = std::string("Time step : ") + std::to_string(dt);
                myScreen.drawText(strDt, Geometry::Point<double>{50, double(myScreen.getGeometry().second - 96)}); // not time consuming
                myScreen.displayVelocityField(grid, vortices);

                myScreen.displayParticles(grid, vortices, cloud); // time consuming

                auto end = std::chrono::system_clock::now();
                std::chrono::duration<double> diff = end - start;
                if (diff.count() >= 1.0)
                {
                    FPS = frameCount;
                    start = end;
                    frameCount = 0;
                }

                std::string str_fps = std::string("FPS : ") + std::to_string(FPS);
                myScreen.drawText(str_fps, Geometry::Point<double>{300, double(myScreen.getGeometry().second - 96)}); // not time consuming
                myScreen.display();                                                                                   // not time consuming
            }
        }
    }

    else
    {
        while (true)
        {
            bool advance = false;
            char command = '_';

            // scaning for message to receive
            MPI_Iprobe(0, 0, MPI_COMM_WORLD, &isReceiving, &status);
            if (isReceiving)
            {
                // reading keyboard pressed
                MPI_Recv(&command, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
                switch (command)
                {
                case 'E':
                    command = 'K';
                    MPI_Send(&command, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD); // process terminated
                    MPI_Finalize();
                    return EXIT_SUCCESS;
                case 'P':
                    animate = true;
                    break;
                case 'S':
                    animate = false;
                    break;
                case 'U':
                    dt *= 2;
                    break;
                case 'D':
                    dt /= 2;
                    break;
                case 'A':
                    advance = true;
                    break;
                default:
                    break;
                }
            }

            if (animate | advance)
            {
                // mobile decide between fix or dynamical execution
                if (isMobile)
                    cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, cloud);
                else
                    cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);

                // grid size is multipled by 2 because each entry of the 2D vector has x and y values
                // number of vortices is multipled by 3 because each vortice has v1, v2 and v3
                // MPI_Isend(cloud.data(), 2 * cloud.numberOfPoints(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
                // MPI_Isend(grid.data(), 2 * grid.getSizeGrid(), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &request);
                // MPI_Isend(vortices.data(), 3 * vortices.numberOfVortices(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &request);
                MPI_Send(cloud.data(), 2 * cloud.numberOfPoints(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                MPI_Send(grid.data(), 2 * grid.getSizeGrid(), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
                MPI_Send(vortices.data(), 3 * vortices.numberOfVortices(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();



    return EXIT_SUCCESS;
}



// mpirun -np 2 ./vortexSimulation.exe ./data/simpleSimulation.dat