#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct vitesse{
	int posit;
	int val;
}Vitesse;

Vitesse *tabvit;

int get_rand(){
	int ran;
	ran = (rand()%6)+1;
	return ran;
}

float get_rand_entre_0_1(){
	int ran;
	ran = (rand()%10);
	return ran;
}

int *get_tab(int **tab_in, int i){
	int *tab;
	int j;
	tab = malloc(sizeof(int)*6);
	for(j=0;j<6;j++){
		tab[j] = tab_in[i][j];
	}
	return tab;
}

int Distance(int i,int j){
	int dist[6][6] = {{0,1,2,3,4,5,6},{1,0,70,60,300,140,100},{2,70,0,90,260,100,50},{3,60,90,0,240,60,170},{4,300,260,240,0,100,200},{5,140,100,60,100,0,100},{6,100,50,170,200,100,0}};
	return dist[0][i]+dist[j][0];
}

bool verifier(int element, int **tab, int linge, int colonne){
	if(colonne==0){
		if(tab[linge][colonne]==element){
			return false;
		}
	}
	else{
		int i;
		for(i=0;i<colonne;i++){
			if(tab[linge][i]==element){
				return false;
			}
		}
	}
	return true;
}

int **initialisation_solutions(){
    int **arr=malloc(sizeof(int *)*60);
    int i,j;
    int rand = get_rand();
    for(i=0;i<40;i++)
    {
        arr[i]=malloc(sizeof(int)*6);
        for(j=0;j<6;j++)
        {
        	while(verifier(rand=get_rand(),arr,i,j)){
        		arr[i][j]=rand;
			}
        }
    }
    return arr;
}

Vitesse *calcule_de_difference(int *tab1,int *tab2){
	tabvit=malloc(sizeof(tabvit));
	int pos=0;
	int i,j;
	for(i=0;i<sizeof(tab1);i++)
    {
    	if(tab1[i]!=tab2[i]){
    		tabvit[i].val =tab2[i];
    		tabvit[i].posit=pos+1;
    		pos++;
		}
		else{
			pos++;
		}
    }
    return tabvit;
}

Vitesse *somme_vitesse(Vitesse *vit1, Vitesse *vit2){
	Vitesse *som = malloc(sizeof(som));
	int i;
	
	for(i=0;i<sizeof(vit1);i++)
    {
    	som[i].val = vit1[i].val;
    	som[i].posit = vit1[i].posit;
    }
    
    for(i=sizeof(vit1);i<sizeof(vit1)+sizeof(vit2);i++)
    {
    	som[i].val = vit2[i-sizeof(vit1)].val;
    	som[i].posit = vit2[i-sizeof(vit1)].posit;
    }
    
    return som;
}

Vitesse *mul_coe_par_vit(Vitesse *vit,float co){
	Vitesse *resu = malloc(sizeof(resu));
	int value_ali[sizeof(vit)];
	int i;
	for(i=0;i<sizeof(vit);i++)
    {
    	value_ali[i] = get_rand_entre_0_1()/10;
    }
    for(i=0;i<sizeof(vit);i++)
    {
    	if(value_ali[i]<co){
    		resu[i].val=vit[i].val;
    		resu[i].posit=vit[i].posit;
		}
    }
    return resu;
}

int fitnesse_particule(int **liste, int i){
	int fit;
	int j;
	for(j=0;j<5;j++){
		fit = fit + Distance(liste[i][j],liste[i][j+1]);
	}
	return fit;
}

int *min(int **liste){
	int *resul = malloc(sizeof(resul));
	int i;
	int min[2];
	for(i=0;i<40;i++){
		resul[i] = fitnesse_particule(liste, i);
	}
	
	//recherche le miminimale avec leur position
	min[0] = resul[0];
	min[1] = 0;
	for(i=1;i<40;i++){
		if(resul[i]<min){
			min[0] = resul[i];
			min[1] = i;
		}
	}
	return *min;
}

int *mise_a_jour_vitesse_position1(Vitesse *vitesse, int pos_act, int best_pos, int best_pos_glob, int **liste){
	Vitesse *vit = malloc(sizeof(vit));
	Vitesse *vit1 = malloc(sizeof(vit1));
	Vitesse *vit2 = malloc(sizeof(vit2));
	Vitesse *vit3 = malloc(sizeof(vit3));
	Vitesse *vit4 = malloc(sizeof(vit4));
	Vitesse *vit5 = malloc(sizeof(vit5));
	
	vit = mul_coe_par_vit(vitesse ,get_rand_entre_0_1()/10);
	vit1 = calcule_de_difference(liste[best_pos],liste[pos_act]);
	vit2 = mul_coe_par_vit(vit1 ,get_rand_entre_0_1()/10);
	vit3 = calcule_de_difference(liste[best_pos_glob],liste[pos_act]);
	vit4 = mul_coe_par_vit(vit3 ,get_rand_entre_0_1()/10);
	
	vit = somme_vitesse(vit, somme_vitesse(vit2 , vit4));
	
	return vit;
}

int **mise_a_jour_vitesse_position2(Vitesse *vitesse, int **liste,int position){
	int **resultat = malloc(sizeof(resultat));
	int i;
	
	for(i=0;i<6;i++){
		if(vitesse[i].posit==i){
			resultat[position][i]=vitesse[i].val;
		}
	}
	
	return resultat;
}

int main()
{
	int **solu_ini;
	int **mettre_a_jour;
	int *minipos;
	Vitesse *vitesse;
	int mini, best_pos;
	int fit_part,fit_best_part;
	int t;
	int i;
	
	minipos = malloc(sizeof(minipos));
	solu_ini = initialisation_solutions();
	mettre_a_jour = malloc(sizeof(mettre_a_jour));
	vitesse = malloc(sizeof(vitesse));
	
	t=0;
	do{
		
		minipos = min(solu_ini);
		mini = minipos[0];
		best_pos = minipos[1];
		
		for(i=0;i<40;i++){
			
			vitesse = mise_a_jour_vitesse_position1(vitesse, i, best_pos, best_pos, solu_ini);
			mettre_a_jour = mise_a_jour_vitesse_position2(vitesse, solu_ini, i);
			fit_part = fitnesse_particule(solu_ini, i);
			fit_best_part = fitnesse_particule(solu_ini, best_pos);
			
			if(fit_part<fit_best_part){
				solu_ini = mettre_a_jour;
			}
			
		}
		
		t=t+1;
		
	}while(t!=60);
	
    return 0;
}
