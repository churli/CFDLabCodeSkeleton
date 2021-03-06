#include "uvp.h"
#include "helper.h"
#include "boundary_val.h"
#include <math.h>

const short XDIR = 0;
const short YDIR = 1;

/**
 * Determines the value of F and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */

void calculate_fg(double Re, double GX, double GY, double alpha, double beta, double dt, double dx, double dy, int imax,
                  int jmax, double **U, double **V, double **F, double **G, double **T, int **Flags)
{
    // Compute F, G on boundaries
    // set boundary conditions for G - see discrete momentum equations - In any case apply Neumann BC - first derivative of pressure must be "zero" - dp/dy = 0
    for (int i = 1; i <= imax; i++)
    {
        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax];
    }
    
    // // set boundary conditions for F - see discrete momentum equations - In any case apply Neumann BC - first derivative of pressure must be "zero" - dp/dx = 0
    for (int j = 1; j <= jmax; j++)
    {
        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];
    }
    
    // calculate F in the domain
    for (int i = 1; i < imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            // We need to compute F only on edges between 2 fluid cells (see p.6 WS2).
            int cell = Flags[i][j];
            if (isObstacle(cell) || isNeighbourObstacle(cell, RIGHT))
            {
                // Boundary condition for F at the obstacle-fluid interface. (or on the obstacle itself)
                F[i][j] = U[i][j];
                continue;
            }
            //
            F[i][j] = computeF(Re, GX, alpha, beta, dt, dx, dy, U, V, T, i, j);
        }
    }
    
    // calculate G in the domain
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j < jmax; j++)
        {
            // We need to compute G only on edges between 2 fluid cells (see p.6 WS2).
            if (isObstacle(Flags[i][j]) || isNeighbourObstacle(Flags[i][j], TOP))
            {
                // Boundary condition for G at the obstacle-fluid interface. (or on the obstacle itself)
                G[i][j] = V[i][j];
                continue;
            }
            //
            G[i][j] = computeG(Re, GY, alpha, beta, dt, dx, dy, U, V, T, i, j);
        }
    }
}

double computeF(double Re, double GX, double alpha, double beta, double dt, double dx, double dy, double **U, double **V, double **T, int i, int j)
{
    return U[i][j] // velocity u
           // diffusive term
           + dt *
             (
                     1 / Re * (secondDerivativeDx(U, i, j, dx) + secondDerivativeDy(U, i, j, dy))
                     // convective term
                     - squareDerivativeDx(U, i, j, dx, alpha)
                     // convective term cont.
                     - productDerivativeDy(U, V, i, j, dy, alpha)
                     // volume force
                     + (1 - beta * T[i][j]) * GX
             );
}

double computeG(double Re, double GY, double alpha, double beta, double dt, double dx, double dy, double **U, double **V, double **T, int i, int j)
{
    return V[i][j] // velocity v
    // diffusive term
    + dt *
      (
              1 / Re * (secondDerivativeDx(V, i, j, dx) + secondDerivativeDy(V, i, j, dy))
              // convective term
              - productDerivativeDx(U, V, i, j, dx, alpha)
              // convective term cont.
              - squareDerivativeDy(V, i, j, dy, alpha)
              // volume force
              + (1 - beta * T[i][j]) * GY
      );
}

double secondDerivativeDx(double **A, int i, int j, double h)
{
    // Approximate the second derivative via central difference.
    // A is the matrix of values.
    // i,j are the coordinates of the central element.
    // h is the discretization step for the chosen direction
    return (A[i - 1][j] - 2 * A[i][j] + A[i + 1][j]) / (h * h);
}

double secondDerivativeDy(double **A, int i, int j, double h)
{
    // Approximate the second derivative via central difference.
    // A is the matrix of values.
    // i,j are the coordinates of the central element.
    // h is the discretization step for the chosen direction
    return (A[i][j - 1] - 2 * A[i][j] + A[i][j + 1]) / (h * h);
}

double productDerivativeDx(double **A, double **B, int i, int j, double h, double alpha)
{
    // Approximate the derivative of the AB product as per formula in the worksheet.
    // A,B are the matrices of values. (Their order is important: A is along x, B along y)
    // i,j are the coordinates of the central element.
    // h is the discretization step for the chosen direction
    return 1 / h *
           (
                   (A[i][j] + A[i][j + 1]) / 2 * (B[i][j] + B[i + 1][j]) / 2
                   - (A[i - 1][j] + A[i - 1][j + 1]) / 2 * (B[i - 1][j] + B[i][j]) / 2
           )
           + alpha / h *
             (
                     fabs(A[i][j] + A[i][j + 1]) / 2 * (B[i][j] - B[i + 1][j]) / 2
                     - fabs(A[i - 1][j] + A[i - 1][j + 1]) / 2 * (B[i - 1][j] - B[i][j]) / 2
             );
}

double productDerivativeDy(double **A, double **B, int i, int j, double h, double alpha)
{
    // Approximate the derivative of the AB product as per formula in the worksheet.
    // A,B are the matrices of values. (Their order is important: A is along x, B along y)
    // i,j are the coordinates of the central element.
    // h is the discretization step for the chosen direction
    return 1 / h *
           (
                   (B[i][j] + B[i + 1][j]) / 2 * (A[i][j] + A[i][j + 1]) / 2
                   - (B[i][j - 1] + B[i + 1][j - 1]) / 2 * (A[i][j - 1] + A[i][j]) / 2
           )
           + alpha / h *
             (
                     fabs(B[i][j] + B[i + 1][j]) / 2 * (A[i][j] - A[i][j + 1]) / 2
                     - fabs(B[i][j - 1] + B[i + 1][j - 1]) / 2 * (A[i][j - 1] - A[i][j]) / 2
             );
}

double squareDerivativeDx(double **A, int i, int j, double h, double alpha)
{
    // Approximate the derivative of the AA product as per formula in the worksheet.
    // A is the matrices of values.
    // i,j are the coordinates of the central element.
    // h is the discretization step for the chosen direction
    return 1 / h *
           (
                   pow((A[i][j] + A[i + 1][j]) / 2, 2)
                   - pow((A[i - 1][j] + A[i][j]) / 2, 2)
           )
           + alpha / h *
             (
                     fabs(A[i][j] + A[i + 1][j]) / 2 * (A[i][j] - A[i + 1][j]) / 2
                     - fabs(A[i - 1][j] + A[i][j]) / 2 * (A[i - 1][j] - A[i][j]) / 2
             );
}

double squareDerivativeDy(double **A, int i, int j, double h, double alpha)
{
    // Approximate the derivative of the AA product as per formula in the worksheet.
    // A is the matrices of values.
    // i,j are the coordinates of the central element.
    // h is the discretization step for the chosen direction
    return 1 / h *
           (
                   pow((A[i][j] + A[i][j + 1]) / 2, 2)
                   - pow((A[i][j - 1] + A[i][j]) / 2, 2)
           )
           + alpha / h *
             (
                     fabs(A[i][j] + A[i][j + 1]) / 2 * (A[i][j] - A[i][j + 1]) / 2
                     - fabs(A[i][j - 1] + A[i][j]) / 2 * (A[i][j - 1] - A[i][j]) / 2
             );
}

/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(double dt, double dx, double dy, int imax, int jmax, double **F, double **G, double **RS, int **Flags)
{
    for (int i = 1; i < imax + 1; i++)
    {
        for (int j = 1; j < jmax + 1; j++)
        {
            if (isFluid(Flags[i][j])) // TODO: double check if this restriction is correct
            {
                RS[i][j] = ((F[i][j] - F[i - 1][j]) / dx + (G[i][j] - G[i][j - 1]) / dy) / dt;
            }
        }
    }
}

/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */

void calculate_dt(
        double Re,
        double Pr,
        double tau,
        double *dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        double **U,
        double **V
)
{
    double u_max = 0, v_max = 0;
    for (int i = 0; i < imax + 1; i++)
    {
        for (int j = 0; j < jmax + 1; j++)
        {
            if (fabs(U[i][j]) > u_max)
            {
                u_max = fabs(U[i][j]);
            }
            if (fabs(V[i][j]) > v_max)
            {
                v_max = fabs(V[i][j]);
            }
        }
    }

    //printf("%f\n", dy / v_max); // todo: can this be removed?
    double minimum = fmin((Re * Pr / 2 / (1 / pow(dx, 2) + 1 / pow(dy, 2))), fmin(dx / u_max, dy / v_max));
    *dt = tau * minimum;
}

/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */

void calculate_uv(double dt, double dx, double dy, int imax, int jmax, double **U, double **V, double **F, double **G,
                  double **P, int **Flags)
{
    for (int i = 1; i < imax; ++i)
    {
        for (int j = 1; j < jmax + 1; ++j)
        {
            int cell = Flags[i][j];
            if (isFluid(cell) && isNeighbourFluid(cell,RIGHT))
            {
                // We need to compute velocity updates only on edges between 2 fluid cells (see p.6 WS2).
                U[i][j] = F[i][j] - (dt / dx * (P[i + 1][j] - P[i][j]));
            }
        }
    }
    for (int i = 1; i < imax + 1; ++i)
    {
        for (int j = 1; j < jmax; ++j)
        {
            int cell = Flags[i][j];
            if (isFluid(cell) && isNeighbourFluid(cell,TOP))
            {
                // We need to compute velocity updates only on edges between 2 fluid cells (see p.6 WS2).
                V[i][j] = G[i][j] - (dt / dy * (P[i][j + 1] - P[i][j]));
            }
        }
    }
}


void calculate_T(double Re, double Pr, double dt, double dx, double dy, double alpha, int imax, int jmax,
                 double **T, double **U, double **V){
    for(int i=0; i < imax+1; ++i){
        for(int j=0; j < jmax+1; ++j){
            T[i][j] = T[i][j] + dt *
                        (
                                - 1/dx * (
                                        U[i][j] * ( T[i][j] + T[i+1][j] ) / 2
                                        - U[i-1][j] * ( T[i-1][j] + T[i][j] ) / 2
                                          )
                                - alpha * 1/dx * (
                                        fabs(U[i][j]) * ( T[i][j] + T[i+1][j] ) / 2
                                        - fabs(U[i-1][j]) * ( T[i-1][j] + T[i][j] ) / 2
                                )

                                - 1/dy * (
                                        V[i][j] * ( T[i][j] + T[i][j+1] ) / 2
                                        - V[i][j-1] * ( T[i][j-1] + T[i][j] ) / 2
                                )
                                - alpha * 1/dy * (
                                        fabs(V[i][j]) * ( T[i][j] + T[i][j+1] ) / 2
                                        - fabs(V[i][j-1]) * ( T[i][j-1] + T[i][j] ) / 2
                                )

                                + 1/(Re * Pr) *
                                (
                                        ( T[i+1][j] - 2 * T[i][j] + T[i-1][j]) / (dx*dx)
                                        + ( T[i][j+1] - 2 * T[i][j] + T[i][j-1]) / (dy*dy)
                                )
                        );
        }
    }
}
