/*
 * nimble_contact_extras.cc
 */

#include <stdexcept>
#include <Geom_Common.h>
#include <Geom_NGP_Interaction_List.h>
#include <Geom_NGP_NodeFacePushback.h>
#include <Geom_Views.h>
#include <nimble_contact_entity.h>
#include <nimble_contact_extras.h>
#include <temp_move_to_stk_search/Morton_Boxes_And_Collisions.hpp>
#include <temp_move_to_stk_search/MortonLBVH_Search.hpp>
#include <temp_move_to_stk_search/MortonLBVH_Tree.hpp>
#include <vector/Vec3.h>

namespace nimble {

void load_contact_points_and_tris(
    const int num_collisions, stk::search::CollisionList<nimble_kokkos::kokkos_device_execution_space> collision_list,
    nimble_kokkos::DeviceContactEntityArrayView contact_nodes_d,
    nimble_kokkos::DeviceContactEntityArrayView contact_faces_d,
    gtk::PointsView<nimble_kokkos::kokkos_device_execution_space> points,
    gtk::TrianglesView<nimble_kokkos::kokkos_device_execution_space> triangles) {
  Kokkos::parallel_for(
      "Load contact points and triangles",
      num_collisions,
      KOKKOS_LAMBDA(const int i_collision) {
        int contact_node_index = collision_list.m_data(i_collision, 0);
        int contact_face_index = collision_list.m_data(i_collision, 1);
        const ContactEntity &node = contact_nodes_d(contact_node_index);
        const ContactEntity &face = contact_faces_d(contact_face_index);
        points.setPointValue(i_collision, node.coord_1_x_, node.coord_1_y_, node.coord_1_z_);
        triangles.setVertexValues(i_collision, mtk::Vec3<double>(face.coord_1_x_, face.coord_1_y_, face.coord_1_z_),
                                  mtk::Vec3<double>(face.coord_2_x_, face.coord_2_y_, face.coord_2_z_),
                                  mtk::Vec3<double>(face.coord_3_x_, face.coord_3_y_, face.coord_3_z_));
      });
}

KOKKOS_INLINE_FUNCTION
void scatter_contact_forces(double gap, double magnitude, double direction[3],
                            const nimble_kokkos::DeviceContactEntityArrayView &contact_faces_d, int contact_face_index,
                            const double *closest_pt,
                            const nimble_kokkos::DeviceContactEntityArrayView &contact_nodes_d, int contact_node_index,
                            const nimble_kokkos::DeviceScalarNodeView &contact_manager_force_d) {
  if (gap < 0.0) {
    double contact_force[3];
    for (int i = 0; i < 3; ++i) {
      contact_force[i] = magnitude * direction[i];
    }
    contact_faces_d(contact_face_index).ComputeNodalContactForces(contact_force, closest_pt);

    for (int i = 0; i < 3; ++i) {
      contact_force[i] *= -1.0;
    }
    contact_nodes_d(contact_node_index).ComputeNodalContactForces(contact_force, closest_pt);
    contact_nodes_d(contact_node_index).ScatterForceToContactManagerForceVector(contact_manager_force_d);
    contact_faces_d(contact_face_index).ScatterForceToContactManagerForceVector(contact_manager_force_d);
  }
}

void compute_and_scatter_contact_force(
    nimble_kokkos::DeviceContactEntityArrayView contact_nodes_d,
    nimble_kokkos::DeviceContactEntityArrayView contact_faces_d,
    stk::search::CollisionList<nimble_kokkos::kokkos_device_execution_space> collision_list,
    nimble_kokkos::DeviceScalarNodeView contact_manager_force_d, double penalty_parameter) {
  using namespace gtk::exp_ngp_contact;

  int numPoints = contact_nodes_d.extent(0);
  gtk::PointsView<nimble_kokkos::kokkos_device_execution_space> points("points", numPoints);
  Kokkos::parallel_for("Load contact_nodes_d into PointsView", numPoints, KOKKOS_LAMBDA(const int i_point) {
    const ContactEntity &node = contact_nodes_d(i_point);
    points.setPointValue(i_point, node.coord_1_x_, node.coord_1_y_, node.coord_1_z_);
  });

  int numTris = contact_faces_d.extent(0);
  gtk::TrianglesView<nimble_kokkos::kokkos_device_execution_space> triangles("triangles", numTris);
  Kokkos::parallel_for(
      "Load contact_faces_d into TrianglesView",
      numTris,
      KOKKOS_LAMBDA(const int i_face) {
    const ContactEntity &face = contact_faces_d(i_face);
    triangles.setVertexValues(i_face, mtk::Vec3<double>(face.coord_1_x_, face.coord_1_y_, face.coord_1_z_),
                              mtk::Vec3<double>(face.coord_2_x_, face.coord_2_y_, face.coord_2_z_),
                              mtk::Vec3<double>(face.coord_3_x_, face.coord_3_y_, face.coord_3_z_));
  });

  NodeFaceInteractions interactionLists = compute_node_face_interaction_lists(numPoints, collision_list, points,
                                                                              triangles);
  PushbackDirectionsAndGaps pushbackAndGaps = face_avg_pushback(interactionLists, points, triangles);
  Kokkos::deep_copy(contact_manager_force_d, 0.0);
  auto pushbackDirs = pushbackAndGaps.directions;
  auto gaps = pushbackAndGaps.gaps;
  auto closestPoints = get_projected_points(interactionLists);

  Kokkos::parallel_for(
      "compute_something",
      numPoints,
      KOKKOS_LAMBDA(const int i_node) {
    double gap = gaps(i_node);
    interactionLists.for_each_face_interaction_of_node(i_node,
                                                       [=](int nodeIdx, int contact_face_index, int dataIdx) {
      double closest_pt[3] = { closestPoints(dataIdx, 0),
          closestPoints(dataIdx, 1), closestPoints(dataIdx, 2) };
      double direction[3] = { pushbackDirs(i_node, 0),
          pushbackDirs(i_node, 1), pushbackDirs(i_node, 2) };
      const double scale = penalty_parameter * gap / interactionLists.num_interactions(i_node);

      scatter_contact_forces(gap, scale, direction,
                             contact_faces_d, contact_face_index,
                             closest_pt, contact_nodes_d, i_node,
                             contact_manager_force_d);
    });
  });
}

void ExtrasContactInterface::ComputeContact(nimble_kokkos::DeviceContactEntityArrayView contact_nodes,
                                            nimble_kokkos::DeviceContactEntityArrayView contact_faces,
                                            nimble_kokkos::DeviceScalarNodeView contact_manager_force,
                                            double penalty_parameter) {
  using namespace gtk::exp_ngp_contact;
  stk::search::CollisionList<nimble_kokkos::kokkos_device_execution_space> collision_list("contact_collision_list");
  stk::search::MortonLBVHSearch_Timers timers;
  stk::search::mas_aabb_tree_loader<double, nimble_kokkos::kokkos_device_execution_space> contact_nodes_tree_loader(
      contact_nodes_search_tree_, contact_nodes.size());
  int num_contact_nodes = contact_nodes.size();
  int num_contact_faces = contact_faces.size();
  Kokkos::parallel_for("Load contact nodes search tree", num_contact_nodes, KOKKOS_LAMBDA(const int i_contact_node) {
    double min_x = contact_nodes(i_contact_node).get_x_min();
    double max_x = contact_nodes(i_contact_node).get_x_max();
    double min_y = contact_nodes(i_contact_node).get_y_min();
    double max_y = contact_nodes(i_contact_node).get_y_max();
    double min_z = contact_nodes(i_contact_node).get_z_min();
    double max_z = contact_nodes(i_contact_node).get_z_max();
    contact_nodes_tree_loader.set_box(i_contact_node, min_x, max_x, min_y, max_y, min_z, max_z);
  });

  stk::search::mas_aabb_tree_loader<double, nimble_kokkos::kokkos_device_execution_space> contact_faces_tree_loader(
      contact_faces_search_tree_, contact_faces.size());
  Kokkos::parallel_for("Load contact faces search tree", num_contact_faces, KOKKOS_LAMBDA(const int i_contact_face) {
    double min_x = contact_faces(i_contact_face).get_x_min();
    double max_x = contact_faces(i_contact_face).get_x_max();
    double min_y = contact_faces(i_contact_face).get_y_min();
    double max_y = contact_faces(i_contact_face).get_y_max();
    double min_z = contact_faces(i_contact_face).get_z_min();
    double max_z = contact_faces(i_contact_face).get_z_max();
    contact_faces_tree_loader.set_box(i_contact_face, min_x, max_x, min_y, max_y, min_z, max_z);
  });

  stk::search::TimedMortonLBVHSearch<double, nimble_kokkos::kokkos_device_execution_space>(contact_nodes_search_tree_,
                                                                                           contact_faces_search_tree_,
                                                                                           collision_list, timers);

  compute_and_scatter_contact_force(contact_nodes, contact_faces, collision_list, contact_manager_force,
                                    penalty_parameter);
}

std::string GTKProjectionTypeToString(short int val)  {
  gtk::ProjectionType type = static_cast<gtk::ProjectionType>(val);
  std::string type_str("Unknown");
  switch (type) {
    case gtk::NODE_PROJECTION:
      type_str = "NODE_PROJECTION";
      break;
    case gtk::EDGE_PROJECTION:
      type_str = "EDGE_PROJECTION";
      break;
    case gtk::FACE_PROJECTION:
      type_str = "FACE_PROJECTION";
      break;
    case gtk::NUM_PROJ_TYPE:
      type_str = "NUM_PROJ_TYPE";
      break;
    case gtk::NULL_PROJECTION:
      type_str = "NULL_PROJECTION";
      break;
    default:
      throw std::logic_error("\nError in GTKProjectionTypeToString(), unrecognized ProjectionType.\n");
      break;
  }
  return type_str;
}

}
