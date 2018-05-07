#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/*
 * Auxiliary data structures to handle the boundary values.
 */
typedef enum BoundarySide
{
    TOPBOUNDARY,
    BOTTOMBOUNDARY,
    LEFTBOUNDARY,
    RIGHTBOUNDARY
} BoundarySide;

typedef enum BoundaryType
{
    DIRICHLET,
    NEUMANN
} BoundaryType;

typedef struct BoundaryInfo
{
    BoundaryType typeU;
    BoundaryType typeV;
    char constU; // 1 means the value to apply is uniform.
    char constV; // 1 means the value to apply is uniform.
    double *valuesU;
    double *valuesV;
} BoundaryInfo;

// Initialize a BoundaryInfo object
void initBoundaryInfo(BoundaryInfo *boundaryInfo, BoundaryType typeU, BoundaryType typeV,
                      int numValuesU, int numValuesV);

/**
 * The boundary values of the problem are set.
 */

void boundaryvalues(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo boundaryInfo[4]);

void setLeftBoundaryVelocities(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo *boundaryInfo);

void setRightBoundaryVelocities(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo *boundaryInfo);

void setTopBoundaryVelocities(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo *boundaryInfo);

void setBottomBoundaryVelocities(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo *boundaryInfo);

#endif
