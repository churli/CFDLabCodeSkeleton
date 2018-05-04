#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "boundary_val.h"
#include "uvp.h"
#include "logger.h"


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

int main(int argn, char** args){

	char* szFileName = "cavity100.dat";
	char* szProblem = "visualization";
	double Re;                /* reynolds number   */
    double UI;                /* velocity x-direction */
    double VI;                /* velocity y-direction */
    double PI;                /* pressure */
    double GX;                /* gravitation x-direction */
    double GY;                /* gravitation y-direction */
    double t_end;             /* end time */
    double xlength;           /* length of the domain x-dir.*/
    double ylength;           /* length of the domain y-dir.*/
    double dt;                /* time step */
    double dx;                /* length of a cell x-dir. */
    double dy;                /* length of a cell y-dir. */
    int  imax;                /* number of cells x-direction*/
    int  jmax;                /* number of cells y-direction*/
    double alpha;             /* uppwind differencing factor*/
    double omg;               /* relaxation factor */
    double tau;               /* safety factor for time step*/
	int  itermax;			  /* max. number of iterations  */
    double eps;               /* accuracy bound for pressure*/
    double dt_value;          /* time for output */
	int n = 0;				  /* timestep iteration counter */
	double res = 10;		  /* residual */
	double t = 0;			  /* initial time */
	int it;					  /* sor iteration counter */
	double mindt=10000;
    
    openLogFile(); // Initialize the log file descriptor.
    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value); 

    double** U = matrix(0, imax+1, 0, jmax+1);
    double** V = matrix(0, imax+1, 0, jmax+1);
    double** F = matrix(0, imax+1, 0, jmax+1);
    double** G = matrix(0, imax+1, 0, jmax+1);
    double** RS = matrix(0, imax+1, 0, jmax+1);
    double** P = matrix(0, imax+1, 0, jmax+1);

    // initialise velocities and pressure
	init_uvp(UI,VI,PI,imax,jmax,U,V,P);
	
	// TODO: Check if this visualization output can be removed!
//	write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
//	n++;
	// simulation interval 0 to t_end
	double currentOutputTime = 0; // For chosing when to output
	while(t < t_end){
		
		// adaptive stepsize control based on stability conditions ensures stability of the method!
		// dt = tau * min(cond1, cond2, cond3) where tau is a safety factor
		// NOTE: if tau<0, stepsize is not adaptively computed!
		if(tau > 0){
			calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
            dt = fmin(dt, dt_value); // test, to avoid a dt bigger than visualization interval
			// Used to check the minimum time-step for convergence
			if (dt < mindt)
				mindt = dt;
		}
		
		// ensure boundary conditions for velocity
		boundaryvalues(imax, jmax, U, V);
//		if(t == 0){
//			write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
//			n++;
//		}
		// momentum equations M1 and M2 - F and G are the terms arising from explicit Euler velocity update scheme
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G);
		
		// momentum equations M1 and M2 are plugged into continuity equation C to produce PPE - depends on F and G - RS is the rhs of the implicit pressure update scheme
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);
		
		// solve the system of eqs arising from implicit pressure uptate scheme using succesive overrelaxation solver
		it = 0;
        res = 1e9;
        while(it < itermax && res > eps){
			sor(omg, dx, dy, imax, jmax, P, RS, &res);
			it++;
		}
        if (it == itermax)
        {
//            printf("[%12.9f] WARNING: max number of iterations reached on SOR. Probably it did not converge!\n", t);
            logEvent(t, "WARNING: max number of iterations reached on SOR. Probably it did not converge!");
        }
		// calculate velocities acc to explicit Euler velocity update scheme - depends on F, G and P
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);
		
		// write visualization file for current iteration (only every dt_value step)
		if (t >= currentOutputTime)
		{
            logEvent(t, "INFO: Writing visualization file n=%d", n);
			write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
			currentOutputTime += dt_value;
			// update output timestep iteration counter
			n++;
		}
        // Recap shell output
//        printf("[%12.9f] INFO: dt=%f, numSorIterations=%d, sorResidual=%f\n", t, dt, it, res);
        logEvent(t, "INFO: dt=%f, numSorIterations=%d, sorResidual=%f", dt, it, res);
		// advance in time
		t += dt;
	}

	// write visualisation file for the last iteration
    logEvent(t, "INFO: Writing visualization file n=%d", n);
    write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

	// Check value of U[imax/2][7*jmax/8] (task6)
    logMsg("Final value for U[imax/2][7*jmax/8] = %16e", U[imax / 2][7 * jmax / 8]);

	free_matrix( U, 0, imax+1, 0, jmax+1);
	free_matrix( V, 0, imax+1, 0, jmax+1);
	free_matrix( F, 0, imax+1, 0, jmax+1);
	free_matrix( G, 0, imax+1, 0, jmax+1);
	free_matrix( RS, 0, imax+1, 0, jmax+1);
	free_matrix( P, 0, imax+1, 0, jmax+1);
    
    logMsg("Min dt value used: %16e", mindt);
    
    closeLogFile(); // Properly close the log file

	return 0;
}

