//
// Created by tommaso on 07/05/18.
//

#ifndef SIM_BOUNDARY_CONFIGURATOR_H
#define SIM_BOUNDARY_CONFIGURATOR_H

void configureBoundary(BoundaryInfo *boundaryInfo, BoundarySide boundarySide,
                       const char *boundaryTypeStr, double u, double v);

#endif //SIM_BOUNDARY_CONFIGURATOR_H
