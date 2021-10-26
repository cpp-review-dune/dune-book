//
// Created by carlosal1015 on 10/26/21.
//
#include <set>

#include <gmsh.h>

int main(int argc, char** argv)
{
    gmsh::initialize();

    gmsh::model::add("t1");

    double lc = 1e-2;
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);

    gmsh::model::geo::addPoint(.1, 0, 0, lc, 2);
    gmsh::model::geo::addPoint(.1, .3, 0, lc, 3);

    int p4 = gmsh::model::geo::addPoint(0, .3, 0, lc);

    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(3, 2, 2);
    gmsh::model::geo::addLine(3, p4, 3);
    gmsh::model::geo::addLine(4, 1, p4);

    gmsh::model::geo::addCurveLoop({4, 1, -2, 3}, 1);

    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(1, {1, 2, 4}, 5);
    int ps = gmsh::model::addPhysicalGroup(2, {1});
    gmsh::model::setPhysicalName(2, ps, "My surface");

    gmsh::model::mesh::generate(2);

    gmsh::write("t1.msh");
    gmsh::finalize();
    return 0;
}