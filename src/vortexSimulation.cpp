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

typedef struct
{
    bool advance, animate;
    double dt;
} eventVars;

eventVars getEvents(Graphisme::Screen myScreen, bool animate, bool advance, float dt)
{
    eventVars vars;
    vars.advance = advance;
    vars.animate = animate;
    vars.dt = dt;

    sf::Event event;
    while (myScreen.pollEvent(event))
    {
        // évènement "fermeture demandée" : on ferme la fenêtre
        if (event.type == sf::Event::Closed)
            myScreen.close();
        if (event.type == sf::Event::Resized)
        {
            // on met à jour la vue, avec la nouvelle taille de la fenêtre
            myScreen.resize(event);
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::P))
            vars.animate = true;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
            vars.animate = false;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
            vars.dt *= 2;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
            vars.dt /= 2;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
            vars.advance = true;
    }

    return vars;
};

int main(int nargs, char *argv[])
{
    MPI_Init(&nargs, &argv);

    int nProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

    int process;
    MPI_Comm_rank(MPI_COMM_WORLD, &process);

    // MPI_Request request;

    // MPI_Status status;
    // if (process == 0)
    // {
    //     MPI_Send(&token, 1, MPI_INT, 1, 101, globComm);
    //     MPI_Recv(&token, 1, MPI_INT, nProcesses - 1, 101, globComm, &status);
    // }
    // else
    // {
    //     MPI_Recv(&token, 1, MPI_INT, process - 1, 101, globComm, &status);
    //     MPI_Send(&token, 1, MPI_INT, (process + 1) % nProcesses, 101, globComm);
    // }
    // MPI_Recv(
    // void* data,
    // int count,
    // MPI_Datatype datatype,
    // int source,
    // int tag,
    // MPI_Comm communicator,
    // MPI_Status* status)

    // int MPI_Irecv(
    //           void* buf,
    //           int count,
    //           MPI_Datatype datatype,
    //           int source,
    //           int tag,
    //           MPI_Comm comm,
    //           MPI_Request *request)

    // MPI_Send(
    // void* data,
    // int count,
    // MPI_Datatype datatype,
    // int destination,
    // int tag,
    // MPI_Comm communicator)

    // int MPI_Isend(
    //           const void* buf,
    //           int count,
    //           MPI_Datatype datatype,
    //           int dest,
    //           int tag,
    //           MPI_Comm comm,
    //           MPI_Request *request)




    bool isActive;




    if (process == 0)
    {
        char const *filename;
        if (nargs == 1)
        {
            std::cout << "Usage : vortexsimulator <nom fichier configuration>" << std::endl;
            return EXIT_FAILURE;
        }

        filename = argv[1];
        std::ifstream fich(filename);
        auto config = readConfigFile(fich);
        fich.close();

        std::size_t resx = 1080, resy = 1080; // square window
        if (nargs > 3)
        {
            resx = std::stoull(argv[2]);
            resy = std::stoull(argv[3]);
        }

        auto vortices = std::get<0>(config);
        auto isMobile = std::get<1>(config);
        auto grid = std::get<2>(config);
        auto cloud = std::get<3>(config);
        double dt;

        std::cout << "######## Vortex simultor ########" << std::endl
                  << std::endl;
        std::cout << "Press P for play animation " << std::endl;
        std::cout << "Press S to stop animation" << std::endl;
        std::cout << "Press right cursor to advance step by step in time" << std::endl;
        std::cout << "Press down cursor to halve the time step" << std::endl;
        std::cout << "Press up cursor to double the time step" << std::endl;

        grid.updateVelocityField(vortices);

        Graphisme::Screen myScreen({resx, resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()});
        dt = 0.1;
        bool animate = false;

        // myScreen({resx, resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()});
        // dt = 0.1;
        // animate = false;

        while(isActive)
        // while (myScreen.isOpen())
        {
        // if (process == 0)
        // {
            auto start = std::chrono::system_clock::now();
            bool advance = false;

            // advance = false;

            // évènements
            // on inspecte tous les évènements de la fenêtre qui ont été émis depuis la précédente itération
            sf::Event event;
            while (myScreen.pollEvent(event))
            {
                // évènement "fermeture demandée" : on ferme la fenêtre
                if (event.type == sf::Event::Closed)
                    myScreen.close();
                    // isActive = myScreen.isOpen();

                if (event.type == sf::Event::Resized)
                {
                    // on met à jour la vue, avec la nouvelle taille de la fenêtre
                    myScreen.resize(event);
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::P))
                    animate = true;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
                    animate = false;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
                    dt *= 2;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
                    dt /= 2;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
                    advance = true;
            }


            // calcule
            if (animate | advance)
            {
                if (isMobile)
                {
                    MPI_Send(&dt,       1, MPI_DOUBLE, 1, 100, MPI_COMM_WORLD);
                    // MPI_Send(&grid,     1, MPI_DOUBLE, 1, 101, MPI_COMM_WORLD);
                    // MPI_Send(&vortices, 1, MPI_DOUBLE, 1, 102, MPI_COMM_WORLD);
                    // MPI_Send(&cloud,    1, MPI_DOUBLE, 1, 103, MPI_COMM_WORLD);
                    cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, cloud);
                    // MPI_Recv(&cloud,    1, MPI_DOUBLE, 1, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                }
                else
                {
                    // MPI_Send(&dt, 1, MPI_DOUBLE, 1, 101, MPI_COMM_WORLD);
                    cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);
                }
            }

            // affichage
            myScreen.clear(sf::Color::Black);

            std::string strDt = std::string("time step : ") + std::to_string(dt);
            myScreen.drawText(strDt, Geometry::Point<double>{50, double(myScreen.getGeometry().second - 96)});
            myScreen.displayVelocityField(grid, vortices);
            myScreen.displayParticles(grid, vortices, cloud);

            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = end - start;

            std::string str_fps = std::string("FPS : ") + std::to_string(1. / diff.count());
            myScreen.drawText(str_fps, Geometry::Point<double>{300, double(myScreen.getGeometry().second - 96)});

            myScreen.display();

            isActive = myScreen.isOpen();
            MPI_Send(&isActive, 1, MPI_CXX_BOOL, 1, 99, MPI_COMM_WORLD);
            // MPI_Isend(&isActive, 1, MPI_CXX_BOOL, 1, 99, MPI_COMM_WORLD, &request);

        }
    }


    if (process == 1)
    {
        isActive = true;

        // char const *filename;

        // filename = argv[1];
        // std::ifstream fich(filename);
        // auto config = readConfigFile(fich);
        // fich.close();

        // auto vortices = std::get<0>(config);
        // auto isMobile = std::get<1>(config);
        // auto grid = std::get<2>(config);
        // auto cloud = std::get<3>(config);
        Simulation::Vortices  vortices;
        bool isMobile;
        Numeric::CartesianGridOfSpeed grid;
        Geometry::CloudOfPoints cloud;

        double dt;


        while(isActive){
            // if (isMobile)
            // {
                MPI_Recv(&dt,       1, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Recv(&grid,     1, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Recv(&vortices, 1, MPI_DOUBLE, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Recv(&cloud,    1, MPI_DOUBLE, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Send(&cloud,    1, MPI_DOUBLE, 0, 104, MPI_COMM_WORLD);

                cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, cloud);

                std::cout << "print dt: " << dt << std::endl;

                MPI_Recv(&isActive, 1, MPI_CXX_BOOL, 0, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Irecv(&isActive, 1, MPI_CXX_BOOL, 0, 99, MPI_COMM_WORLD, &request);
                // std::cout << "isActive: " << isActive << std::endl;

            // }
        }

    }

    MPI_Finalize();



    return EXIT_SUCCESS;
}



// mpirun -np 2 ./vortexSimulation.exe ./data/simpleSimulation.dat