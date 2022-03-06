#include "library.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define RNG_UNIFORM() (rand()/(double)RAND_MAX)

#define RNG_UNIFORM_INT(s) (rand()%s)
typedef struct Particles
{
    int position;
    int velocity; //
    int col;
}Particles;

typedef Particles  *part;

int c1=2;
int c2=2;

int villes[5][5]={
        {0 ,10 ,50 ,200 ,300 },
        {10 ,0 ,90 ,160 ,180 },
        {50 ,90 ,0 ,240 ,60 },
        {200 ,160 ,240 ,0 ,70 },
        {300 ,180 ,60 ,70 ,0 }
};
int particules[40][5];
int pBest[40][5];
int GBest[5];
int cost[40];

bool exist(int T[], int val){
    for(int i =0;i<5;i++){
        if(T[i]==val)
        {
            return true;
        }
    }
    return false;
}
int *InitialiseParticles(int p){
    int value;
    int v=0;
    for(int i = 0;i<5;i++){
        value=rand()%5+1;
        if (exist(&particules[p][i],value))
        {
            particules[p][i]=value;
            pBest[p][i]=value;
            v++;
        }else{
            while(exist(particules[p],value))
            {
                value= rand()%5+1;
            }
            particules[p][i]=value;
            pBest[p][i]=value;
        }
    }
    return *particules;
}
int villesDistance(int t1, int t2){
    int dist =0;
    dist = villes[t1-1][t2-1];
    return dist;
}
int particleCost(int p){
    int coste=0;
    for(int i=0;i<4;i++){
        coste+= villesDistance(particules[p][i],particules[p][i+1]);
    }
    coste+= villesDistance(particules[p][6], particules[p][0]);

    return coste;
}
int costPBest(int p){
    int coste=0;
    for(int i=0;i<4;i++){
        coste+= villesDistance(pBest[p][i],pBest[p][i+1]);
    }
    coste+= villesDistance(pBest[p][6], pBest[p][0]);

    return coste;
}
int costGBest(){
    int coste=0;
    for(int i=0;i<4;i++){
        coste+= villesDistance(GBest[i],GBest[i+1]);
    }
    coste+= villesDistance(GBest[6], GBest[0]);

    return coste;
}
part *tabP;
part *tmpP;
part *Somme_PV;
part *vitesse_calc;

part calDif(int **T,int **M,int i,int j){
    int v=0;
    for(int i=0;i<5;i++){
        if(T[i][j]!=M[i][j])
        {
            tabP[v]=(part)(malloc(sizeof (Particles)));
            tabP[v]->position=i;
            tabP[v]->velocity=M[j][i];
            tabP[0]->col=v+1;
            v++;
        }
    }
    return *tabP;
}
part MultipliciteCoe(float c, part parti, part pParticles) {
    int random = RNG_UNIFORM();
    random = random * c;

    float *tab;
    for(int i =0;i<parti[0].col;i++){
        tab[i]=*(float *)malloc(sizeof(float));
        tab[i]=RNG_UNIFORM();
    }
    for(int i=0;i<parti[0].col;i++){
        if(tab[i]<random){
            tmpP[i]=(part)malloc(sizeof(Particles));
            tmpP[i]->position=parti[i].position;
            tmpP[i]->velocity=parti[i].velocity;
        }
    }
    return *tmpP;
}
part sommeVitesse(part *tab1, part *tab2)
{
    for (int j = 0; j < tab1[0]->col; j++)
    {
        Somme_PV[j] = (part)malloc(sizeof(Particles));
        Somme_PV[j]->position = tab1[j]->position;
        Somme_PV[j]->velocity = tab1[j]->velocity;
    }
    for (int j = tab1[0]->col; j < tab1[0]->col + tab2[0]->col; j++)
    {
        Somme_PV[j] = (part)malloc(sizeof(Particles));
        Somme_PV[j]->position = tab2[j]->position;
        Somme_PV[j]->velocity = tab2[j]->velocity;
    }
    return Somme_PV;
}

int *sommeVitessePos(part *tab, int p)
{
    int *t;

    for (int i = 0; i < 5; i++)
    {
        t[i] = particules[p][i];
    }
    for (int i = 0; i < tab[0]->col; i++)
    {
        t[(int)tab[i]->position] =(int) (tab[(int)tab[i]->position] + (int)tab[i]->velocity);
    }

    return t;
}
void pso()
{
    int t = 0;
    int fitnesseParticule = 0;
    for (int i = 0; i < 40; i++)
    {
        InitialiseParticles(i);
    }

    do
    {
        for (int j = 0; j < 40; j++)
        {
            if (j == 0)
            {
                for (int i = 0; i < 6; i++)
                {
                    GBest[i] = particules[j][i];
                }
            }
            else if (costGBest() > particleCost((int)particules[j]))
            {
                for (int i = 0; i < 6; i++)
                {
                    GBest[i] = particules[j][i];
                }
            }

            vitesse_calc = sommeVitesse((part)MultipliciteCoe(c1, calDif(particules[0], particules[1], 0, 1), NULL),(part)MultipliciteCoe(c1, calDif(particules[j], particules[j], j, j),
                                                                                                                                          MultipliciteCoe(c2, calDif(particules[j], GBest[j], j, j), NULL)
                                        )
            );

            for (int i = 0; i < 6; i++)
            {
                particules[j][i] = sommeVitessePos(vitesse_calc, j)[i];
            }
            fitnesseParticule = particleCost(j);
            if (particleCost(j) < costPBest(j))
            {
                for (int i = 0; i < 6; i++)
                {
                    pBest[j + 1][i] = particules[j][i];
                }
            }
        }
        t++;
    } while (t < 60);
}
