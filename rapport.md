# OS202_project

à la suite de ce document on va présenter une analyse des données rencontrées pendant l'exécution des étapes suivantes:
- **partie 1**: séparation;
- **partie 2**: parallélisation en mémoire distribuée;
- **partie 3**: parallélisation en mémoire distribuée et partagée des calculs;
- **partie 4**: analyse approche Eulérienne-Lagrangienne;
pour faire la comparasion entre des différents l'ensemble de données utilisé était le:

on considère comme reference les données suivantes pour mon ordniateur:
```bash
Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         39 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  20
  On-line CPU(s) list:   0-19
Vendor ID:               GenuineIntel
  Model name:            12th Gen Intel(R) Core(TM) i7-12700H
    CPU family:          6
    Model:               154
    Thread(s) per core:  2
    Core(s) per socket:  14
    Socket(s):           1
    Stepping:            3
    CPU max MHz:         4700.0000
    CPU min MHz:         400.0000


Caches (sum of all):     
  L1d:                   544 KiB (14 instances)
  L1i:                   704 KiB (14 instances)
  L2:                    11.5 MiB (8 instances)
  L3:                    24 MiB (1 instance)

```
on utilisera le FPS comme mesure de la vélocité du code.


## partie 1, 2 et 3
avec une fenêtre de taille `540 x 540 pixels` avec un `dt = 0.1` et on considère les performance suivantes:
| code | main | part1 | part2 | part3 |
|------|-----:|------:|------:|------:|
| oneVortexSimulation | 8 FPS| 7 FPS| 12 FPS| - FPS|
| cornerTest | 62 FPS| 59 FPS| 92 FPS| - FPS|
| simpleSimulation | 21 FPS| 18 FPS| 28 FPS| - FPS|
| manyVortices | 11 FPS| 10 FPS| 14 FPS| - FPS|

avec une fenêtre de taille `1080 x 1080 pixels` avec un `dt = 0.1` et on considère les performance suivantes:
| code | main | part1 | part2 | part3 |
|------|-----:|------:|------:|------:|
| oneVortexSimulation | 8 FPS| 7 FPS| 12 FPS| - FPS|
| cornerTest | 60 FPS| 58 FPS| 82 FPS| - FPS|
| simpleSimulation | 21 FPS| 19 FPS| 26 FPS| - FPS|
| manyVortices | 11 FPS| 10 FPS| 14 FPS| - FPS|

avec une fenêtre de taille `2160 x 2160 pixels` avec un `dt = 0.1` et on considère les performance suivantes:
| code | main | part1 | part2 | part3 |
|------|-----:|------:|------:|------:|
| oneVortexSimulation | 8 FPS| 7 FPS| 12 FPS| - FPS|
| cornerTest | 58 FPS| 55 FPS| 82 FPS| - FPS|
| simpleSimulation | 20 FPS| 19 FPS| 27 FPS| - FPS|
| manyVortices | 10 FPS| 10 FPS| 14 FPS| - FPS|

on note que la part1 n'a pas ajouté des performances. l'adition des processus de communication entre les processus a fait que le FPS reste très pareil.

même que l'exécution soit parallélise, gagne de performance, les exchanges entre les processus fait qui, à la fin, l'exécution soit pareil.


par contre, on note que la part2 a augmente beaucoup l'efficacité du algorithme. les calculus parallelises rendent le code beaucoup plus rapide comment on peut voir avec le gagne de FPS.

on note que la taille de la simulation n'a pas apporté trop de changement de temps d'exécution.



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
