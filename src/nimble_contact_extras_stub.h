/*
 * nimble_contact_extras_stub.h
 */

#ifndef SRC_NIMBLE_CONTACT_EXTRAS_STUB_H_
#define SRC_NIMBLE_CONTACT_EXTRAS_STUB_H_

#ifdef NIMBLE_HAVE_EXTRAS

#include <stdexcept>

#include "Geom_Views.h"

namespace gtk
{

template<typename ExecSpace, bool save_proj>
class ComputeProjections
{
public:
ComputeProjections(gtk::PointsView<ExecSpace> points,
                   gtk::TrianglesView<ExecSpace> tris,
                   gtk::PointsView<ExecSpace> closest_points,
                   Kokkos::View<short*, ExecSpace> proj_types_returned_d)
{
    throw std::runtime_error("Error! gtk::ComputeProjections not currently supported.");
}

~ComputeProjections() = default;

};

}

#endif

#endif /* SRC_NIMBLE_CONTACT_EXTRAS_STUB_H_ */
