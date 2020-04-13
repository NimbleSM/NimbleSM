/*
 * nimble_contact_extras.h
 */

#ifndef SRC_NIMBLE_CONTACT_EXTRAS_H_
#define SRC_NIMBLE_CONTACT_EXTRAS_H_

#include <string>
#include <nimble_kokkos_contact_defs.h>
#include <nimble_kokkos_defs.h>
#include <nimble_extras_contact_includes.h>

namespace nimble {

std::string GTKProjectionTypeToString(short int val);

class ExtrasContactInterface {
 public:
  ExtrasContactInterface()
      : contact_nodes_search_tree_("contact nodes search tree"),
        contact_faces_search_tree_("contact faces search tree") {
  }

  ~ExtrasContactInterface() = default;

  void ComputeContact(nimble_kokkos::DeviceContactEntityArrayView contact_nodes, nimble_kokkos::DeviceContactEntityArrayView contact_faces,
                      nimble_kokkos::DeviceScalarNodeView contact_manager_force, double penalty_parameter);

  stk::search::mas_aabb_tree<double, nimble_kokkos::kokkos_device_execution_space> contact_nodes_search_tree_;
  stk::search::mas_aabb_tree<double, nimble_kokkos::kokkos_device_execution_space> contact_faces_search_tree_;
};

}

#endif /* SRC_NIMBLE_CONTACT_EXTRAS_H_ */
