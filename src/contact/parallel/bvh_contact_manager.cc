/*
//@HEADER
// ************************************************************************
//
//                                NimbleSM
//                             Copyright 2018
//   National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
// retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
// NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact David Littlewood (djlittl@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include "bvh_contact_manager.h"

#include <vt/transport.h>

#include <fstream>
#include <iomanip>

#include "bvh/collision_object.hpp"
#include "nimble_data_manager.h"
#include "nimble_model_data.h"
#include "nimble_vector_communicator.h"
#ifdef NIMBLE_HAVE_KOKKOS
#include <bvh/narrowphase/kokkos.hpp>
#endif

#ifdef NIMBLE_HAVE_ARBORX
#include <ArborX.hpp>
#include <Kokkos_Core.hpp>

#include "contact/arborx_utils.h"
#endif

namespace nimble {

namespace details {

#ifdef NIMBLE_HAVE_ARBORX

struct ArborX_PairData
{
  int pred_rank_  = 0;
  int pred_index_ = -1;
  int prim_rank_  = 0;
  int prim_index_ = -1;
  //
  ArborX_PairData(int r_, int i_, int f_r_, int f_i_) : pred_rank_(r_), pred_index_(i_), prim_rank_(f_r_), prim_index_(f_i_) {}
  //
  ArborX_PairData() = default;
  //
  ~ArborX_PairData() = default;
  //
  bool
  operator<(const ArborX_PairData& rhs) const
  {
    return (
        std::tie(pred_rank_, pred_index_, prim_rank_, prim_index_) <
        std::tie(rhs.pred_rank_, rhs.pred_index_, rhs.prim_rank_, rhs.prim_index_));
  }
};

struct ArborXCallback
{
  const nimble_kokkos::HostContactEntityUnmanagedConstView &d_contact_faces;
  const nimble_kokkos::HostContactEntityUnmanagedConstView &d_contact_nodes;

  double enforcement_penalty;

  std::vector<NarrowphaseResult> &resa_vec;
  std::vector<NarrowphaseResult> &resb_vec;

  template <typename Query>
  KOKKOS_FUNCTION void
  operator()(Query const& query, int j) const
  {
    const auto& myNode = d_contact_nodes(ArborX::getData(query));
    const auto& myFace = d_contact_faces(j);

    bool   inside               = false;
    double normal[3] = {0.0, 0.0, 0.0};

    //--- Determine whether the node is projected inside the triangular face
    //--- Determine whether the node is projected inside the triangular face
    NarrowphaseResult entry;
    ContactManager::Projection(myNode, myFace, inside, entry.gap, &normal[0], entry.bary);
    if (inside) {
      details::getContactForce(enforcement_penalty, entry.gap, normal, entry.contact_force);
      //
      entry.local_index = myFace.local_id();
      entry.node        = false;
      resa_vec.push_back(entry);
      //
      entry.local_index = myNode.local_id();
      entry.node        = true;
      resb_vec.push_back(entry);
    }
  }
  //
  template <typename Predicate, typename OutputFunctor>
  KOKKOS_FUNCTION void
  operator()(Predicate const& pred, int f_primitive, OutputFunctor const& out) const
  {
    //
    const auto& myFace = d_contact_faces(f_primitive);
    //
    auto const& p_data     = ArborX::getData(pred); // <- type nimble_kokkos::details::OutputData
    auto const& p_geometry = ArborX::getGeometry(pred);  // <- type Box
    ContactEntity myNode(ContactEntity::ContactEntityType::NODE, 0, p_data.coord_, 0.0, 0);
    //
    bool   inside               = false;
    double normal[3] = {0.0, 0.0, 0.0};
    //
    //--- Determine whether the node is projected inside the triangular face
    //
    NarrowphaseResult entry;
    ContactManager::Projection(myNode, myFace, inside, entry.gap, &normal[0], entry.bary);
    if (inside) {
      details::getContactForce(enforcement_penalty, entry.gap, normal, entry.contact_force);
      //
      entry.local_index = myFace.local_id();
      entry.node        = false;
      resa_vec.push_back(entry);
      //
      entry.local_index = p_data.index_;
      entry.node        = true;
      resb_vec.push_back(entry);
    }
  }
};
#endif

}

struct NarrowphaseFunc
{
  explicit NarrowphaseFunc(BvhContactManager* cm) : contact_manager{cm} {}

#ifndef NIMBLE_HAVE_KOKKOS
  bvh::narrowphase_result_pair
  operator()(const bvh::broadphase_collision<ContactEntity>& _a, const bvh::broadphase_collision<ContactEntity>& _b)
  {
    auto res   = bvh::narrowphase_result_pair();
    res.a      = bvh::narrowphase_result(sizeof(NarrowphaseResult));
    res.b      = bvh::narrowphase_result(sizeof(NarrowphaseResult));
    auto& resa = static_cast<bvh::typed_narrowphase_result<NarrowphaseResult>&>(res.a);
    auto& resb = static_cast<bvh::typed_narrowphase_result<NarrowphaseResult>&>(res.b);
    auto  tree = build_snapshot_tree_top_down(_a.elements);

    std::size_t j = 0;
    for (auto&& elb : _b.elements) {
      query_tree_local(tree, elb, [&_a, &_b, &elb, &resa, &resb, this, j](std::size_t _i) {
        const auto&       face = _a.elements[_i];
        const auto&       node = elb;
        NarrowphaseResult entry;

        bool   hit = false;
        double norm[3];
        ContactManager::Projection(node, face, hit, entry.gap, norm, entry.bary);

        if (hit) {
          details::getContactForce(contact_manager->GetPenaltyForceParam(), entry.gap, norm, entry.contact_force);

          entry.local_index = face.local_id();
          entry.node        = false;
          resa.emplace_back(entry);

          entry.local_index = node.local_id();
          entry.node        = true;
          resb.emplace_back(entry);
        }
      });
      ++j;
    }

    return {resa, resb};
  }
#else
  bvh::narrowphase_result_pair
  operator()(const bvh::broadphase_collision<ContactEntity>& _a, const bvh::broadphase_collision<ContactEntity>& _b)
  {
    auto res   = bvh::narrowphase_result_pair();
    res.a      = bvh::narrowphase_result(sizeof(NarrowphaseResult));
    res.b      = bvh::narrowphase_result(sizeof(NarrowphaseResult));
    auto& resa = static_cast<bvh::typed_narrowphase_result<NarrowphaseResult>&>(res.a);
    auto& resb = static_cast<bvh::typed_narrowphase_result<NarrowphaseResult>&>(res.b);

    auto view_a = nimble_kokkos::HostContactEntityUnmanagedConstView( _a.elements.data(), _a.elements.size() );
    auto view_b = nimble_kokkos::HostContactEntityUnmanagedConstView( _b.elements.data(), _b.elements.size() );
    std::vector<NarrowphaseResult> resa_vec, resb_vec;
    using memory_space = nimble_kokkos::kokkos_host_mirror_memory_space;
    //ArborX::BVH<memory_space> a_bvh{nimble_kokkos::kokkos_host_execution_space{}, view_a};
    //
    // ArborX has the option to turning off/on a sorting of the predicates (= contact nodes)
    // It has also the option to provide an estimate of the number of resulting contact to speed up the query.
    // For these two features, a `TraversalPolicy` variable has to be specified as an input
    // parameter for the query.
    //
    // ArborX::Experimental::TraversalPolicy policy;
    //
    //a_bvh.query(nimble_kokkos::kokkos_host_execution_space{}, view_b,
    //            details::ArborXCallback{view_a, view_b, contact_manager->GetPenaltyForceParam(), resa_vec, resb_vec});
    //
    ArborX::DistributedTree<memory_space> a_bvh(MPI_COMM_SELF, nimble_kokkos::kokkos_host_execution_space{}, view_a);
    Kokkos::View<nimble_kokkos::details::OutputData<>*, nimble_kokkos::kokkos_host> results("results", 0);
    Kokkos::View<int*, nimble_kokkos::kokkos_host>                 offset("offset", 0);
    a_bvh.query(
        nimble_kokkos::kokkos_host_execution_space{},
        nimble_kokkos::details::PredicateTypeNodesRank_Host{view_b, 0},
        details::ArborXCallback{view_a, view_b, contact_manager->GetPenaltyForceParam(), resa_vec, resb_vec},
        results,
        offset);
    //
    resa.set_data(resa_vec.data(), resa_vec.size());
    resb.set_data(resb_vec.data(), resb_vec.size());
    
    //-----------------
    /*
    static int iter_count = 0;
    if ( resa.size() > 0 )
      std::cout << iter_count << ". resa = " << resa.size() << "\n";
    if ( resb.size() > 0 )
      std::cout << iter_count << ". resb = " << resb.size() << "\n";
    ++iter_count;
    */
    //-----------------

    return {resa, resb};
  }
#endif

  BvhContactManager* contact_manager;
};

BvhContactManager::BvhContactManager(
    std::shared_ptr<ContactInterface> interface,
    nimble::DataManager&              data_manager,
    std::size_t                       _overdecomposition)
    : ParallelContactManager(interface, data_manager),
      m_world{_overdecomposition},
      m_nodes{&m_world.create_collision_object()},
      m_faces{&m_world.create_collision_object()}
{
  m_world.set_narrowphase_functor<ContactEntity>(NarrowphaseFunc{this});
}

void
BvhContactManager::ComputeBoundingVolumes()
{
  for (auto&& node : contact_nodes_) {
    const double           inflation_length = node.inflation_factor * node.char_len_;
    ContactEntity::vertex* v                = reinterpret_cast<ContactEntity::vertex*>(&node.coord_1_x_);
    node.kdop_                              = bvh::bphase_kdop::from_sphere(*v, inflation_length);
  }
  for (auto&& face : contact_faces_) {
    const double           inflation_length = face.inflation_factor * face.char_len_;
    ContactEntity::vertex* v                = reinterpret_cast<ContactEntity::vertex*>(&face.coord_1_x_);
    face.kdop_                              = bvh::bphase_kdop::from_vertices(v, v + 3, inflation_length);
  }
}

void
BvhContactManager::ComputeParallelContactForce(int step, bool debug_output, nimble::Viewify<2> contact_force)
{
  auto model_ptr = this->data_manager_.GetMacroScaleData();

  auto field_ids    = this->data_manager_.GetFieldIDs();
  auto displacement = model_ptr->GetVectorNodeData(field_ids.displacement);
  this->ApplyDisplacements(displacement.data());

  total_search_time.Start();
  m_world.start_iteration();

  // Update collision objects, this will build the trees
  ComputeBoundingVolumes();

  m_nodes->set_entity_data(contact_nodes_);
  m_faces->set_entity_data(contact_faces_);

  m_faces->broadphase(*m_nodes);

  m_last_results.clear();
  m_faces->for_each_result<NarrowphaseResult>(
      [this](const NarrowphaseResult& _res) { m_last_results.emplace_back(_res); });

  m_world.finish_iteration();
  total_search_time.Stop();

  // event = vt::theTrace()->registerUserEventColl(“name”)
  // tr = vt::trace::TraceScopedNote scope{“my node”, event};
  // theTrace()->addUserEventBracketed(event, start, stop)
  // theTrace()->addUserBracketedNote(start, stop, “my node”, event)
  total_enforcement_time.Start();
  for (auto& f : force_) f = 0.0;

  // Update contact entities
  for (auto&& r : m_last_results) {
    if (r.node) {
      if (r.local_index >= contact_nodes_.size())
        std::cerr << "contact node index " << r.local_index << " is out of bounds (" << contact_nodes_.size() << ")\n";
      auto& node = contact_nodes_.at(r.local_index);
      node.set_contact_status(true);
      node.SetNodalContactForces(r.contact_force);
      node.ScatterForceToContactManagerForceVector(force_);
    } else {
      if (r.local_index >= contact_faces_.size())
        std::cerr << "contact face index " << r.local_index << " is out of bounds (" << contact_faces_.size() << ")\n";
      auto& face = contact_faces_.at(r.local_index);
      face.set_contact_status(true);
      face.SetNodalContactForces(r.contact_force, r.bary);
      face.ScatterForceToContactManagerForceVector(force_);
    }
  }

  total_num_contacts += m_last_results.size();
  total_enforcement_time.Stop();

  this->GetForces(contact_force.data());

  constexpr int vector_dim           = 3;
  auto          myVectorCommunicator = this->data_manager_.GetVectorCommunicator();
  myVectorCommunicator->VectorReduction(vector_dim, contact_force.data());
}

namespace {

void
WriteContactNodesToVTKFile(const std::vector<ContactEntity>& nodes, const std::string& prefix, int step)
{
  std::stringstream file_name_ss;
  file_name_ss << prefix;
  file_name_ss << std::setfill('0') << std::setw(5) << step;
  file_name_ss << ".vtk";
  std::ofstream vis_file;
  vis_file.open(file_name_ss.str().c_str());

  // file version and identifier
  vis_file << "# vtk DataFile Version 3.0" << std::endl;

  // header
  vis_file << "Contact entity visualization" << std::endl;

  // file format ASCII | BINARY
  vis_file << "ASCII" << std::endl;

  // dataset structure STRUCTURED_POINTS | STRUCTURED_GRID | UNSTRUCTURED_GRID | POLYDATA | RECTILINEAR_GRID | FIELD
  vis_file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  std::vector<int> vtk_vertices(nodes.size());

  vis_file << "POINTS " << nodes.size() << " float" << std::endl;

  for (int i = 0; i < nodes.size(); i++) {
    vis_file << nodes[i].coord_1_x_ << " " << nodes[i].coord_1_y_ << " " << nodes[i].coord_1_z_ << std::endl;
    vtk_vertices[i] = i;
  }

  vis_file << "CELLS " << vtk_vertices.size() << " " << 2 * vtk_vertices.size() << std::endl;
  for (const auto& vtk_vertex : vtk_vertices) { vis_file << "1 " << vtk_vertex << std::endl; }

  vis_file << "CELL_TYPES " << nodes.size() << std::endl;
  for (unsigned int i = 0; i < vtk_vertices.size(); ++i) {
    // cell type 1 is VTK_VERTEX
    vis_file << 1 << std::endl;
  }

  vis_file.close();
}

void
WriteContactFacesToVTKFile(const std::vector<ContactEntity>& faces, const std::string& prefix, int step)
{
  std::stringstream file_name_ss;
  file_name_ss << prefix;
  file_name_ss << std::setfill('0') << std::setw(5) << step;
  file_name_ss << ".vtk";
  std::ofstream vis_file;
  vis_file.open(file_name_ss.str().c_str());

  // file version and identifier
  vis_file << "# vtk DataFile Version 3.0" << std::endl;

  // header
  vis_file << "Contact entity visualization" << std::endl;

  // file format ASCII | BINARY
  vis_file << "ASCII" << std::endl;

  // dataset structure STRUCTURED_POINTS | STRUCTURED_GRID | UNSTRUCTURED_GRID | POLYDATA | RECTILINEAR_GRID | FIELD
  vis_file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  std::vector<std::vector<int>> vtk_triangles(faces.size());

  vis_file << "POINTS " << 3 * faces.size() << " float" << std::endl;

  for (int i = 0; i < faces.size(); i++) {
    ContactEntity const& face = faces[i];
    vis_file << face.coord_1_x_ << " " << face.coord_1_y_ << " " << face.coord_1_z_ << std::endl;
    vis_file << face.coord_2_x_ << " " << face.coord_2_y_ << " " << face.coord_2_z_ << std::endl;
    vis_file << face.coord_3_x_ << " " << face.coord_3_y_ << " " << face.coord_3_z_ << std::endl;
    vtk_triangles[i].push_back(3 * i);
    vtk_triangles[i].push_back(3 * i + 1);
    vtk_triangles[i].push_back(3 * i + 2);
  }

  vis_file << "CELLS " << vtk_triangles.size() << " " << 4 * vtk_triangles.size() << std::endl;
  for (const auto& vtk_triangle : vtk_triangles) {
    vis_file << "3 " << vtk_triangle[0] << " " << vtk_triangle[1] << " " << vtk_triangle[2] << std::endl;
  }

  vis_file << "CELL_TYPES " << faces.size() << std::endl;
  for (unsigned int i = 0; i < vtk_triangles.size(); ++i) {
    // cell type 6 is VTK_TRIANGLE
    vis_file << 6 << std::endl;
  }

  vis_file.close();
}

void
WriteContactEntitiesToVTKFile(
    const std::vector<ContactEntity>& faces,
    const std::vector<ContactEntity>& nodes,
    const std::string&                prefix,
    int                               step)
{
  std::stringstream file_name_ss;
  file_name_ss << prefix;
  file_name_ss << std::setfill('0') << std::setw(5) << step;
  file_name_ss << ".vtk";
  std::ofstream vis_file;
  vis_file.open(file_name_ss.str().c_str());

  // file version and identifier
  vis_file << "# vtk DataFile Version 3.0" << std::endl;

  // header
  vis_file << "Contact entity visualization" << std::endl;

  // file format ASCII | BINARY
  vis_file << "ASCII" << std::endl;

  // dataset structure STRUCTURED_POINTS | STRUCTURED_GRID | UNSTRUCTURED_GRID | POLYDATA | RECTILINEAR_GRID | FIELD
  vis_file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  std::vector<int>              vtk_vertices(nodes.size());
  std::vector<std::vector<int>> vtk_triangles(faces.size());

  vis_file << "POINTS " << nodes.size() + 3 * faces.size() << " float" << std::endl;

  for (int i = 0; i < nodes.size(); i++) {
    vis_file << nodes[i].coord_1_x_ << " " << nodes[i].coord_1_y_ << " " << nodes[i].coord_1_z_ << std::endl;
    vtk_vertices[i] = i;
  }
  int offset = static_cast<int>(vtk_vertices.size());
  for (int i = 0; i < faces.size(); i++) {
    ContactEntity const& face = faces[i];
    vis_file << face.coord_1_x_ << " " << face.coord_1_y_ << " " << face.coord_1_z_ << std::endl;
    vis_file << face.coord_2_x_ << " " << face.coord_2_y_ << " " << face.coord_2_z_ << std::endl;
    vis_file << face.coord_3_x_ << " " << face.coord_3_y_ << " " << face.coord_3_z_ << std::endl;
    vtk_triangles[i].push_back(offset + 3 * i);
    vtk_triangles[i].push_back(offset + 3 * i + 1);
    vtk_triangles[i].push_back(offset + 3 * i + 2);
  }

  vis_file << "CELLS " << vtk_vertices.size() + vtk_triangles.size() << " "
           << 2 * vtk_vertices.size() + 4 * vtk_triangles.size() << std::endl;
  for (const auto& vtk_vertex : vtk_vertices) { vis_file << "1 " << vtk_vertex << std::endl; }
  for (const auto& vtk_triangle : vtk_triangles) {
    vis_file << "3 " << vtk_triangle[0] << " " << vtk_triangle[1] << " " << vtk_triangle[2] << std::endl;
  }

  vis_file << "CELL_TYPES " << nodes.size() + faces.size() << std::endl;
  for (unsigned int i = 0; i < vtk_vertices.size(); ++i) {
    // cell type 1 is VTK_VERTEX
    vis_file << 1 << std::endl;
  }
  for (unsigned int i = 0; i < vtk_triangles.size(); ++i) {
    // cell type 6 is VTK_TRIANGLE
    vis_file << 6 << std::endl;
  }

  vis_file.close();
}
}  // namespace
}  // namespace nimble
