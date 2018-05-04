#include "boundary_val.h"

void boundaryvalues(int imax, int jmax, double **U, double **V) {
	
	for (int j = 1; j <= jmax; j++) {
		// left boundary - no-slip condition
		U[0][j] = 0;
		// right boundary - no-slip condition
		U[imax][j] = 0;
		// left boundary - no-slip condition - no v component on right wall of the cell - we need to enforce the condition v = 0 at the boundary by interpolation
		V[0][j] = -V[1][j];
		// right boundary - no-slip condition - no v component on right wall of the cell - we need to enforce the condition v = 0 at the boundary by interpolation
		V[imax+1][j] = -V[imax][j];
	}

	for (int i = 1; i <= imax; i++) {
		// lower boundary - no-slip condition
		V[i][0] = 0;
		// upper boundary - moving wall condition - set v = 0
		V[i][jmax] = 0;
		// lower boundary - no-slip condition - no u component on the horizontal wall of the cell - we need to enforce the condition u = 0 at the boundary by interpolation
		U[i][0] = -U[i][1];
		// upper boundary - moving wall condition - no u component on the horizontal wall of the cell - we need to enforce the condition u = 1 at the boundary by interpolation
		U[i][jmax+1] = 2.0 - U[i][jmax];
	}
}

//eof
