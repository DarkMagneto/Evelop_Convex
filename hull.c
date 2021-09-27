// HAOUCHE Achour, TP1 B

#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>

#define BUFSIZE 1024

//____________________________ Partie 1 : Structures pour la géométrie ___________________________

struct vec {
  double x;
  double y;
};

// ----------------  Question 1.1 ----------------
/*
 * Fonction dot qui calcule le produit scalaire
 * @param const struct vec *v1
 * @param const struct vec *v2
 * @return double produit scalaire
 */
double dot(const struct vec *v1, const struct vec *v2){
  return (v1->x * v2->x) + (v1->y * v2->y);
}

// ----------------  Question 1.2 ----------------
/*
 * Fonction cross qui calcule le produit vectoriel 2D de deux vecteurs (P1,P2) et (P1,P3)
 * @param const struct vec *p1
 * @param const struct vec *p2
 * @param const struct vec *p3
 * @return double produit vectoriel
 */
double cross(const struct vec *p1,const struct vec *p2, const struct vec *p3){
  return ((p2->x - p1->x) * (p3->y - p1->y)) - ((p2->y - p1->y) * (p3->x - p1->x));
}

// ----------------  Question 1.3 ----------------
/*
 * Fonction qui dit si la suite de point P1,P2,P3 constitue un tournant à gauche
 * Si le produit vectoriel est plus grand que zero alors P3 est a droite
 * @param const struct vec *p1
 * @param const struct vec *p2
 * @param const struct vec *p3
 * @return bool tourne a gauche
 */
bool is_left_turn(const struct vec *p1,const struct vec *p2, const struct vec *p3){
  return cross(p1,p2,p3) > 0;
}

//_______________________________ Partie 2 : Ensemble de points _____________________________________

struct vecset {
  struct vec *data;
  size_t size;
  size_t capacity;
};

// ----------------  Question 2.1 ----------------
/*
 * Fonction qui créer un ensemble de points vide
 * @param const struct vecset *self
 */
void vecset_create(struct vecset *self){
  self -> data = malloc(10 * sizeof(struct vec));
  self -> capacity = 10;
  self -> size = 0;
}

// ----------------  Question 2.2 ----------------
/*
 * Fonction qui détruit un ensemble de points
 * @param const struct vecset *self
 */

void vecset_destroy(struct vecset *self){
  free(self -> data);
  self -> capacity = 0;
  self -> size = 0;
}

/*
 * Augmente la capacite du vecset self
 * @param struct vecset *self
 */
void vecset_grow(struct vecset *self){
  size_t capacity = self->capacity*2;
  struct vec *data = calloc(capacity,sizeof(struct vec));
  memcpy(data,self->data,self->size*sizeof(struct vec));
  free(self->data);
  self->data = data;
  self->capacity = capacity;
}

// ----------------  Question 2.3 ----------------
/*
 * Fonction qui ajoute un point à un ensemble de points
 * @param struct vecset *self
 * @param struct vec p
 */
void vecset_add(struct vecset *self, struct vec p){
  if (self->size == self->capacity) {
    vecset_grow(self);
  }
  self->data[self->size] = p;
  self->size++;
}

typedef int (*comp_func_t)(const struct vec *p1,const struct vec *p2, const void *ctx);

/*
 * Fonction qui compare les abcsisses de deux points p1 et p2
 * @param const struct vec *p1
 * @param const struct vec *p2
 * @param const void *ctx
 */
int vec_compare(const struct vec *p1, const struct vec *p2, const void *ctx){
  if (p1->x < p2->x) return -1;
  if (p1->x > p2->x) return 1;
  return 0;
}

/*
 * Fonction qui compare les abcsisses et ordonnées de deux points 
 * @param const struct vec *p1
 * @param const struct vec *p2
 * @param const void *ctx
 */
int vec_position_compare(const struct vec *p1, const struct vec *p2, const void *ctx){
  if (p1->y < p2->y) return -1;
  if (p1->y > p2->y) return 1;
  vec_compare(p1, p2, ctx);
  return 0;
}

// ----------------  Question 2.4 ----------------
/*
 * Fonction qui renvoie le maximum d’un ensemble de points suivant une fonction de comparaison donnée
 * @param struct vecset *self
 * @param comp_func_t func
 * @param const void *ctx
 * @return const struct vec
 */ 
const struct vec *vecset_max(const struct vecset *self, comp_func_t func, const void *ctx){
  struct vec *max = &self -> data[0];
  for (size_t i = 1; i < self -> size; ++i) {
    if (func(max, &self -> data[i], &ctx) < 0) {
      max = &self -> data[i];
    }
  }
  return max;
}

// ----------------  Question 2.5 ----------------
/*
 * Fonction qui renvoie le minimum d’un ensemble de points suivant une fonction de comparaison donnée
 * @param struct vecset *self
 * @param comp_func_t func
 * @param const void *ctx
 * @return const struct vec
 */
const struct vec *vecset_min(const struct vecset *self, comp_func_t func, const void *ctx){
  struct vec *min = &self -> data[0];
  for (size_t i = 1; i < self -> size; ++i) {
    if (func(min, &self -> data[i], &ctx) > 0) {
      min = &self -> data[i];
    }
  }
  return min;
}

// ----------------  Question 2.6 ----------------
/*
 * function swap beetwin data[i] and data[j]
 * @param struct vecset *self
 * @param size_t i
 * @param size_t j
 *
 */
void vecset_swap(struct vecset *self, size_t i, size_t j){
  struct vec swap = self->data[i];
  self->data[i] = self->data[j];
  self->data[j] = swap;
}

/*
 * Fonction vecset_partition
 * @param struct vecset *self
 * @param comp_func_t func
 * @param ptrdiff_t i
 * @param ptrdiff_t j
 * @param const void *ctx
 */
ptrdiff_t vecset_partition(struct vecset *self, ptrdiff_t i, ptrdiff_t j, comp_func_t func, const void *ctx){
  ptrdiff_t  pivot_index = i;
  const struct vec pivot = self->data[pivot_index];
  vecset_swap(self, pivot_index, j);
  ptrdiff_t  l = i;

  for (ptrdiff_t  k = i; k < j; ++k) {
    if (func(&self->data[k], &pivot, ctx) < 0) {
      vecset_swap(self, k, l);
      l++;
    }
  }
  vecset_swap(self, l, j);
  return l;
}

/*
 * Fonction qui fait un appel recurssif pour partitionée le tableau
 * @param struct vecset *self
 * @param comp_func_t func
 * @param ptrdiff_t i
 * @param ptrdiff_t j
 * @param const void *ctx
 */
void vecset_quick_sort_partial(struct vecset *self, ptrdiff_t i, ptrdiff_t j, comp_func_t func, const void *ctx){
  if (i < j){
    ptrdiff_t p = vecset_partition(self, i, j, func, ctx);
    vecset_quick_sort_partial(self, i, p - 1, func, ctx);
    vecset_quick_sort_partial(self, p + 1, j, func, ctx);
  }
}

/*
 * Fonction qui trie l’ensemble de points suivant la fonction de comparaison donnée
 * @param struct vecset *self
 * @param comp_func_t func
 * @param const void *ctx
 */
void vecset_sort(struct vecset *self, comp_func_t func,const void *ctx){
  ptrdiff_t size = self->size;
  vecset_quick_sort_partial(self, 0, size - 1, func, ctx);
}

// ----------------  Question 2.7 ----------------
/*
 * Fonction similaire vecset_add
 * @param struct vecset *self
 * @param struct vec p
 */
void vecset_push(struct vecset *self, struct vec p){
  vecset_add(self, p);
}

// ----------------  Question 2.8 ----------------
/*
 * Fonction *vecset_pop
 * @param const struct vecset *self
 */
void vecset_pop(struct vecset *self){
  struct vec * data = calloc(self->capacity, sizeof(struct vec));
  memcpy(data, self->data, (self->size - 1) * sizeof(struct vec));
  free(self->data);
  self->data = data;
  self->size--;
}

// ----------------  Question 2.9 ----------------
/*
 * Fonction *vecset_top qui renvoie le premier élément de la pile
 * @param const struct vecset *self
 * @return const struct 
 */
const struct vec *vecset_top(const struct vecset *self){
  return &self -> data[self -> size - 1];
}

// ----------------  Question 2.10 ----------------
/*
 * Fonction *vecset_second
 * @param const struct vecset *self
 * @return const struct
 */
const struct vec *vecset_second(const struct vecset *self){
  return &self -> data[self -> size - 2];
}


//____________________________________ Partie 3 : Marche de Jarvis ______________________________


/*
 * Fonction de comparaison entre deux points pour obtenir celui le plus a gauche des abcsisses
 * @param const struct vecset *self
 * @param const void *ctx
 * @return int
 */
const struct vec *vecset_leftmost_point(const struct vecset *self, const void *ctx){
  return vecset_min(self, &vec_compare, NULL);
}

// ----------------  Question 3.1 ----------------
/* Fonction jarvis_march
 * @param const struct vecset *in
 * @param struct vecset *out
 */
void jarvis_march(const struct vecset *in, struct vecset *out){
  
  const struct vec *vec_first = vecset_leftmost_point(in, NULL);
  const struct vec *vec_current = vec_first;
  struct vec * vec_next;
  int point_in_s = 0;

  do {
    vecset_add(out, *vec_current);

    //On utilise rand() pour avoir un index de tableau in.data au hasard qui change a chaque debut de boucle 
    //et qui ne peut pas etre repris une 10eme fois lorsqu'il a ete stocké dans le tableau out.data
    do {
      point_in_s = rand() % (in -> size - 1);
    } while (&in -> data[point_in_s] == vec_current);
    vec_next = &in -> data[point_in_s];

    for (int i = 0; i < in -> size; ++i) { 
      if (is_left_turn(vec_current, &in -> data[i], vec_next)) {
        vec_next = &in->data[i];
      }
    }
    vec_current = vec_next;
  } while (vec_first != vec_current);
} 


//____________________________________ Partie 4 : Parcours de Graham _____________________________


/*
 * Fonction vec_angle_calcule calcul l'angle
 * @param const struct vec *p1
 * @param const struct vec *p2
 * @param const void *ctx
 */
float vec_angle_calcule(const struct vec *p1, const struct vec *p2, const void *ctx){
  //atan2 : L'arc tangente du quotient formé par les deux arguments, 
  //c'est-à-dire l'angle, exprimé en radian entre l'axe des abscisses et la droite passant par l'origin (0,0) et le point de coordonnées (x,y).
  const struct vec *p3 = ctx;
  /*  
                      |
                      |
                      p3        p3
                      |
                      |
                      p1              p1
                      |
                      |
                      |
                      |
  ____________________|___x_____p3____p1_________
                      |
                      |
                      x   x
                      |
                      |
                      |
                      |
                      |
                      |
                      |
                      |
  */
  //(p1->x - p3->x), (p1->y - p3->y) nous donne un angle en valeur absolu egale a l'angle p1,0,P3
  //(p2->x - p3->x), (p2->y - p3->y) nous donne un angle en valeur absolu egale a l'angle p2,0,P3
  //Don il y a une premiere arc tangente de l'angle p1,0,P3 et p2,0,P3,
  //on soustrait le premier angle du 2eme pour avoir la valeur de l'angle p1,p3,p2
  float angle = atan2((p1->x - p3->x), (p1->y - p3->y)) - atan2((p2->x - p3->x), (p2->y - p3->y));
  return angle;
}

/*
 * Fonction vec_angle_compare 
 * @param const struct vec *p1
 * @param const struct vec *p2
 * @param const void *ctx*
 * @return int
 */
//si l'angle est positif alors le point est a gauche
//si l'angle est negatif le point est a droite
//si les angles sont egaux a zero alors nos deux vecteurs sont superposés donc on fait une comparaison sur les positions des points des vecteurs
int vec_angle_compare(const struct vec *p1, const struct vec *p2, const void *ctx){
  if(vec_angle_calcule(p1,p2,ctx) > 0)return 1;  
  if(vec_angle_calcule(p1,p2,ctx) < 0)return -1;   
  return vec_position_compare(p1, p2, ctx);
}

/*
 * Fonction copy du vecset car on ne peut pas modifier le vecset in dans le main car c'est un const
 * @param const struct vecset *in
 * @param const struct vecset *out
 * @param const void *ctx
 */
void vecset_copy(const struct vecset * in, struct vecset *out, const void *ctx){
  size_t pos_current = 0;
  do{
    vecset_add(out, in->data[pos_current]);
    pos_current++;
  }while(pos_current < in->size);
}

// ----------------  Question 4.1 ----------------
/* Fonction graham_scan
 * @param const struct vecset *in
 * @param struct vecset *out
 */
void graham_scan(const struct vecset * in, struct vecset * out){
  const struct vec *vec_low = vecset_min(in, &vec_position_compare, NULL);
  //On créer une copy de in que l'on pourra modifier
  struct vecset *copy = malloc(sizeof(struct vecset));
  vecset_create(copy);
  vecset_copy(in,copy,vec_low);

  vecset_sort(copy, &vec_angle_compare, vec_low);

  const struct vec *vec_first = vecset_top(out);
  const struct vec *vec_top;
  const struct vec *vec_second;

  vecset_push(out, *vec_low);
  vecset_push(out, *vec_first);

  for (size_t i = 0; i < copy->size; ++i) {
    vec_top = vecset_top(out);
    vec_second = vecset_second(out);

    while ( out->size >= 2 && is_left_turn(vec_second, vec_top, &copy->data[i])) {
      vecset_pop(out);
      //On réeinitialise vec_top et vec second car out change donc c'est deux derniers aussi 
      //si on ne le fait pas, on garde les premieres valeurs qui deviennent fausses
      vec_top = vecset_top(out);    
      vec_second = vecset_second(out);
    }     
    vecset_push(out, copy->data[i]);
  }
  vecset_destroy(copy);
  free(copy);
}


//____________________________________ Partie 5 : Enveloppe Rapide _____________________________

/*
 * Fonction vecset_distance_between_point_vector
 * @param const struct vec *X
 * @param const struct vec *Y
 * @param double x
 * @param double y
 * @param const void *ctx
 * @return double
 * formule de calcule : sqrt(pow2(p2.x-p1.x) + pow2(p2.y-p1.y))
 */
double vecset_distance_between_point_vector(double x, double y, const struct vec *X, const struct vec *Y, const void *ctx){
  //coefficient directeur cf de la droite(XY)
  double cf = (X -> y - Y -> y) / (X -> x - Y -> x);
  //on cherche l'ordonnée oo à l'origine en utilisant les coordonnées du point X ou du point Y 
  //qui vérifient l'équation y = cf*x + oo et on en déduit oo en résolvant cette équation. 
  double oo = X -> y - cf * X -> x;
  return fabs(cf * x - y + oo) / sqrt(pow(cf, 2) + 1);
}

/*
 * Fonction qui trouve le point le plus loin d'un segment
 * @param struct vecset *in
 * @param const struct vec *X
 * @param const struct vec *Y
 * @param const void *ctx
 * @return struct vec fartherest_point le point le plus loin
 */
struct vec *vecset_fartherest_point(struct vecset *in, const struct vec *X, const struct vec *Y, const void *ctx){
  struct vec *fartherest_point = &in->data[0];
  double distance_i_moins_1 = 0;

  for (size_t i = 0; i < in->size; i++) {
    if (vecset_distance_between_point_vector(in->data[i].x, in->data[i].y, X, Y, NULL) > distance_i_moins_1) {
      fartherest_point = &in->data[i];
      distance_i_moins_1 = vecset_distance_between_point_vector(in->data[i].x, in->data[i].y, X, Y, NULL);
    }
  }
  return fartherest_point;
}

/*
 * Fonction qui dit si l'ensemble de vecteur est vide ou non
 * @param const struct vec *in
 * @return bool
 */
bool vecset_is_empty(const struct vecset *in){
  return in -> size == 0 ;
}

/*
 * Fonction compare deux points pour obtenir celui le plus a droite des abscisses
 * @param const struct vecset *self
 * @param const void *ctx
 * @return const struct vec
 */
const struct vec *vecset_rightmost_point(const struct vecset *self, const void *ctx){
  return vecset_max(self, &vec_compare, NULL);
}

/*
 * Fonction inverse de is_left_turn 
 * @param const struct vec *p1
 * @param const struct vec *p2
 * @param const struct vec *p3
 * @return bool
 */
bool is_right_turn(const struct vec *p1, const struct vec *p2, const struct vec *p3){
  return cross(p1,p2,p3) < 0;
}


// ----------------  Question 5.1 ----------------

/*
 * Fonction findhull
 * @param struct vecset *in
 * @param const struct vec *X
 * @param const struct vec *Y
 * @return struct vecset R
 */
struct vecset *findhull(struct vecset *in, const struct vec *X, const struct vec *Y){
 
  struct vecset *vecset_vide = malloc(sizeof(struct vecset));
  vecset_create(vecset_vide);

  if (vecset_is_empty(in)) {
    return vecset_vide;
  }

  const struct vec *farthest_point = vecset_fartherest_point(in, X, Y, NULL);
  struct vecset *S1 = malloc(sizeof(struct vecset));
  vecset_create(S1);
  struct vecset *S2 = malloc(sizeof(struct vecset));
  vecset_create(S2);

  for (size_t i = 0; i < in->size; i++) {
    const struct vec cdi = in->data[i];
    if (is_left_turn(X, farthest_point, &cdi)) {
      vecset_add(S1, in->data[i]);
    }
    if (is_right_turn(Y, farthest_point, &cdi)) {
      vecset_add(S2, in->data[i]);
    }
  }

  struct vecset *R1 = findhull(S1, X, farthest_point);
  struct vecset *R2 = findhull(S2, farthest_point, Y);

  struct vecset *R = malloc(sizeof(struct vecset));
  vecset_create(R);

  for (size_t i = 0; i < R1->size; i++) {
    vecset_add(R, R1->data[i]);
  }

  vecset_add(R, *farthest_point);

  for (size_t i = 0; i < R2->size; i++) {
    vecset_add(R, R2->data[i]);
  }
  vecset_destroy(vecset_vide);
  free(vecset_vide);

  vecset_destroy(S1);
  free(S1);
  vecset_destroy(S2);
  free(S2);
  vecset_destroy(R1);
  free(R1);
  vecset_destroy(R2);
  free(R2);


 return R;
}

/*
 * Fonction quickhull
 * @param const struct vecset *in
 * @param struct vecset *out
 */
void quickhull(const struct vecset *in, struct vecset *out){

  const struct vec A = *vecset_leftmost_point(in, NULL);
  const struct vec B = *vecset_rightmost_point(in, NULL);

  struct vecset *S1 = malloc(sizeof(struct vecset));
  vecset_create(S1);
  struct vecset *S2 = malloc(sizeof(struct vecset));
  vecset_create(S2);

  for (size_t i = 0; i < in->size; i++) {
    //On utilise un const struct vec pour la valeur in.data[i] car expected ‘const struct vec *’
    //but argument is of type ‘struct vec’ 
    const struct vec data_i = in->data[i];
    if (is_left_turn(&A, &B, &data_i)) {
      vecset_push(S1, in->data[i]);
    }
    if(is_left_turn(&B, &A, &data_i)){
      vecset_push(S2, in->data[i]);
    }
  }

  struct vecset *R1 = findhull(S1, &A, &B);
  struct vecset *R2 = findhull(S2, &B, &A);

  vecset_add(out, A);

  for (size_t i = 0; i < R1->size; i++) {
    vecset_add(out, R1->data[i]);
  }   
  vecset_add(out, B);

  for (size_t i = 0; i < R2->size; i++) {
    vecset_add(out, R2->data[i]);
  }

  vecset_destroy(S1);
  free(S1);
  vecset_destroy(S2);
  free(S2);
  vecset_destroy(R1);
  free(R1);
  vecset_destroy(R2);
  free(R2);
}


//____________________________________ Partie 6 : Pilote _____________________________

int main(){

  //Partie initialisation
  setbuf(stdout, NULL); 
  struct vecset *self = malloc(sizeof(struct vecset));
  vecset_create(self);
  struct vecset *output = malloc(sizeof(struct vecset));
  vecset_create(output);

  //Partie recuperation stdin dans buffer
  char buffer[BUFSIZE];
  fgets(buffer, BUFSIZE, stdin);
  size_t count = strtol(buffer, NULL, 10);

  //Partie copie des coordonnées genérées par hull-generator dans vecset self
  for (size_t i = 0; i < count; ++i) {
    struct vec p;
    fgets(buffer, BUFSIZE, stdin);
    char *endptr = buffer;
    p.x = strtod(endptr, &endptr);
    p.y = strtod(endptr, &endptr);
    vecset_add(self,p);
  }

  //Partie Enveloppe Convexe
  //jarvis_march(self,output);
  //graham_scan(self,output);
  quickhull(self,output);


  //Partie Affichage
  //%zu pour afficher les size_t size
  printf("%zu\n", output -> size);
  //%f pour afficher les double x et y
  for(size_t i = 0; i < output -> size;i++){
    printf("%f %f\n", output -> data[i].x, output -> data[i].y);  
  }
               
  //Partie netoyage et destruction
  vecset_destroy(self);
  free(self);
  self=NULL;

  vecset_destroy(output);
  free(output);
  output=NULL; 
  
  return 0;
}
