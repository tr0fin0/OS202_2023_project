# OS202_project

à la suite de ce document on va présenter une analyse des données rencontrées pendant l'exécution des étapes suivantes:
- **partie 1**: séparation;
- **partie 2**: parallélisation en mémoire distribuée;
- **partie 3**: parallélisation en mémoire distribuée et partagée des calculs;
- **partie 4**: analyse approche Eulérienne-Lagrangienne;


## partie 1
on considère comme reference les données suivantes:
```bash

```
on considère que le code séquentielle a la performance suivante:


## partie 2


## partie 3


## partie 4
### eulérien
le méthode Eulérien décrit une fluide à partir des champs de vélocité et de pression.

pour faire le calcule d'un champs décrit pour cette hypothèse le champ de vitesse calcule par sous-domaine où chaque processus MPI serait responsable pour une sous-ensemble du champs de vitesse.

avec cette division du champs les processus ont besoin d'echanger les informations de bornés pour que les autres puissent faire les calcules nécessaires en demandant une importante gestion de communication.

cette méthode est plus efficace pour décrire une ensemble grande de données comme des océans et la circulation atmospheric.


### lagrangienne
le méthode Lagrangienne décrit une fluide à partir des particules individuelles.

pour faire le calcule les particulles sont distribuées entre différents processus à partir de la région où ils sont dans le champs en demandant une échange entre les processus pour qu'une particule change de processus au fur et à mesure de l'exécution. 

cette méthode est plus efficace pour étudier les trajectories individuales de chaque particule si l'ensemble de données n'est pas très significative.

### proposition
on considère la structure suivante pour décrire la mise en oeuvre de l'approche mixte Eulérien-Lagrangienne:
```c++
#include <mpi.h>
#include <vector>

// Function to compute the velocity field in a subdomain using an Eulerian approach
void computeVelocityField(std::vector<double>& velocityField) {
    // ...
}

// Function to advect particles through the velocity field using a Lagrangian approach
void advectParticles(std::vector<double>& particles, const std::vector<double>& velocityField) {
    // ...
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set up the domain decomposition
    int numSubdomains = size;
    std::vector<int> subdomainStartIndices(numSubdomains);
    std::vector<int> subdomainSizes(numSubdomains);
    // Compute the start index and size of each subdomain
    // ...

    // Compute the velocity field in each subdomain using an Eulerian approach
    std::vector<double> velocityFields[numSubdomains];
    for (int i = 0; i < numSubdomains; i++) {
        velocityFields[i].resize(subdomainSizes[i]);
        computeVelocityField(velocityFields[i]);
    }

    // Initialize the particles in each subdomain
    std::vector<double> particles[numSubdomains];
    for (int i = 0; i < numSubdomains; i++) {
        particles[i].resize(subdomainSizes[i]);
        // Initialize the particles
        // ...
    }

    // Advect the particles through the velocity field using a Lagrangian approach
    for (int timestep = 0; timestep < numTimesteps; timestep++) {
        for (int i = 0; i < numSubdomains; i++) {
            advectParticles(particles[i], velocityFields[i]);
        }
        // Communicate particles between subdomains
        // ...
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
```
quand l'ensemble de très grande dimension le temps de communication entre les processus devient important et l’efficacité de l'algorithme sera moins important et par consequence l'effort de réaliser la parallélisation ne serait pas compense à cause de:
- **load balancing**: le lagrangienne dépende de la trajectoire de chaque points mais l’eulérien dépende du champs de vélocité sur une région, pendant le derroulement de l'exécution du code il pourra avoir une desequilibre entre la quantité de calcule necessaire à faire pour chaque processus;
- **communication overhead**: si le champs de vélocité est très heterogene ou il y a une grande numéro de particules les exchanges entre les processus peut causer une temps significatif d'attend ou une blocage quand les processus ont besoin d'echanger;
- **implementation complexity**: le mise en oeuvre du code peut le rendre très difficile à comprendre et à le modifier car il y aura beaucoup de communication entre processus et des restritions à suivre;
