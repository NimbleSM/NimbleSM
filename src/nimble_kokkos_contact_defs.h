/*
 * nimble_kokkos_contact_defs.h
 */

#ifndef SRC_NIMBLE_KOKKOS_CONTACT_DEFS_H_
#define SRC_NIMBLE_KOKKOS_CONTACT_DEFS_H_

#include <Kokkos_View.hpp>
#include <nimble_contact_entity.h>
#include <nimble_kokkos_defs.h>

namespace nimble_kokkos {

using DeviceContactEntityArrayView = Kokkos::View< nimble::ContactEntity*, nimble_kokkos::kokkos_layout, nimble_kokkos::kokkos_device >;
using HostContactEntityArrayView = Kokkos::View< nimble::ContactEntity*, nimble_kokkos::kokkos_layout, nimble_kokkos::kokkos_host >;

}

#endif /* SRC_NIMBLE_KOKKOS_CONTACT_DEFS_H_ */
