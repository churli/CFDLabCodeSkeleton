#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "boundary_val.h"
#include "uvp.h"
#include <stdio.h>


int main(int argn, char** args){
	int** Flag = read_pgm("P1.pgm");

	for(int i = 0; i < 53; i++){
		for(int j = 0; j < 53; j++){
			printf("%d", Flag[i][j]);
		}
		printf("\n");
	}

	free_imatrix(Flag, 0, 52, 0, 52);

}