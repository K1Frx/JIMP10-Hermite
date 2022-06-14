#include "makespl.h"
#include "gaus/piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

double Hermit(int n, double x) {
	if (n == 0)
		return 1;
	else if(n == 1)
		return 2*x;
	else
		return 2*x* Hermit(n-1,x) - 2*(n-1)* Hermit(n-2,x);
}

/* Funkcja podstawowa, która zwraca wielomian Hermite'a
Wykorzystany wzór:
H.n(x) = 2x * H.n-1(x) - 2n * H.n-2(x)

Z założeniami że:
H.0(x) = 1
H.1(x) = 2x

Gdzie: 
n to stopień wielomianu
x argument funkcji 

Źródło: https://pl.wikipedia.org/wiki/Wielomiany_Hermite’a*/


double ddx( int n, double x) {
	if (n == 0)
		return 0;
	else if (n == 1)
		return 2;
	else
		return 2* Hermit(n-1,x) + 2*x* ddx(n-1,x) - 2*(n-1)* ddx(n-2,x);
}
/* Pierwsza pochodna
Wzory na pochodne "zaczęrpnięte" z informacji znalezionych w internecie
Chyba gdzieś na stronie AGH, ale zgubiłem źródło i pewności nie mam
Wzory sprawdziłem licząc je na kartce dla n E {1,2,3,4,5} i x = 1 (bo najszybciej)
Wyniki programu pokrywają się z wynikami na kartce, więc wzory wydają się poprawne
Nie do końca rozumiem matematycznie dlaczego ten wzór działa, ale działa, więc z niego skorzystam*/

double ddx2( int n, double x) {
	if(n == 0)
		return 0;
	else if (n == 1)
		return 0;
	else
		return 4* ddx(n-1,x) + 2*x* ddx2(n-1,x) - 2*(n-1)* ddx2(n-2,x);
}
/* Druga pochodna */

double ddx3( int n, double x) {
	if(n == 0)
		return 0;
	if(n == 1)
		return 0;
	else
		return 6* ddx2(n-1,x) + 2*x* ddx3(n-1,x) - 2*(n-1)* ddx3(n-2,x);
}
/* Trzecia pochodna */

/*int main(){
    for(int i=1; i<6; i++){
        printf("Dla n = %d\n",i);
        printf("Hermit: %f\n",Hermit(i,1));
        printf("1 Pochodna: %f\n", ddx(i,1));
        printf("2 Pochodna: %f\n", ddx(i,1));
        printf("3 Pochodna: %f\n\n", ddx(i,1));
    }
    return 0;
}

Ten main pomagał mi, kiedy testowałem czy wzory zostały dobrze przepisane
Zostawiam go w komentarzu, nie wiem w sumie po co, ale nie szkodzi*/


void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
    char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

    for (j = 0; j < nb; j++) {
		    for (i = 0; i < nb; i++)
			    for (k = 0; k < pts->n; k++)
				    add_to_entry_matrix(eqs, j, i, Hermit(i, x[k]) * Hermit(j, x[k]));

		    for (k = 0; k < pts->n; k++)
			    add_to_entry_matrix(eqs, j, nb, y[k] * Hermit(j, x[k]));
    }

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * Hermit(k, xx);
				spl->f1[i] += ck * ddx(k, xx);
				spl->f2[i] += ck * ddx2(k, xx);
				spl->f3[i] += ck * ddx3(k, xx);
			}
		}
	}
}


