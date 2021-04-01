/*
 * nimble_kokkos_contact_defs.h
 */

#ifndef SRC_NIMBLE_KOKKOS_CONTACT_DEFS_H_
#define SRC_NIMBLE_KOKKOS_CONTACT_DEFS_H_

#include <nimble_contact_entity.h>
#include <nimble_kokkos_defs.h>

#include <Kokkos_View.hpp>

namespace nimble_kokkos {

using DeviceContactEntityArrayView = Kokkos::View<
    nimble::ContactEntity*,
    nimble_kokkos::kokkos_layout,
    nimble_kokkos::kokkos_device>;
using HostContactEntityArrayView = Kokkos::View<
    nimble::ContactEntity*,
    nimble_kokkos::kokkos_layout,
    nimble_kokkos::kokkos_host>;

}  // namespace nimble_kokkos

#endif /* SRC_NIMBLE_KOKKOS_CONTACT_DEFS_H_ */
