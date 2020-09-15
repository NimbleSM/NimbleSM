#include "bvh_contact_manager.h"
#include <vt/transport.h>
#include <fstream>
#include <iomanip>

namespace nimble {

  struct NarrowphaseFunc {
    explicit NarrowphaseFunc( BvhContactManager *cm )
      : contact_manager{cm}
    {}

    bvh::narrowphase_result_pair operator()(const bvh::broadphase_collision< ContactEntity > &_a,
        const bvh::broadphase_collision< ContactEntity > &_b )
    {
      auto res = bvh::narrowphase_result_pair();
      res.a = bvh::narrowphase_result( sizeof( NarrowphaseResult ) );
      res.b = bvh::narrowphase_result( sizeof( NarrowphaseResult ) );
      auto &resa = static_cast< bvh::typed_narrowphase_result< NarrowphaseResult > & >( res.a );
      auto &resb = static_cast< bvh::typed_narrowphase_result< NarrowphaseResult > & >( res.b );
      auto tree = build_snapshot_tree_top_down( _a.elements );

      std::size_t j = 0;
      for ( auto &&elb : _b.elements ) {
        query_tree_local(tree, elb, [&_a, &_b, &elb, &resa, &resb, this, j](std::size_t _i) {
          const auto &face = _a.elements[_i];
          const auto &node = elb;
          NarrowphaseResult entry;

          bool hit = false;
          double norm[3];
          ContactManager::Projection(node, face, hit, entry.gap, norm, entry.bary );

          if ( hit )
          {
            details::getContactForce(contact_manager->GetPenaltyForceParam(), entry.gap, norm, entry.contact_force);

            std::cout << "Hit! force: " << entry.contact_force[0] << ", " << entry.contact_force[1] << ", " << entry.contact_force[2] << '\n';
            std::cout << "\tgap: " << entry.gap << ", norm: " << norm[0] << ", " << norm[1] << ", " << norm[2] << '\n';
            std::cout << "\tcharacteristic length: " << face.char_len_ << '\n';

            entry.local_index = face.local_id();
            entry.node = false;
            resa.emplace_back( entry );

            entry.local_index = node.local_id();
            entry.node = true;
            resb.emplace_back( entry );
          }
        } );
        ++j;
      }

      return { resa, resb };
    }

    BvhContactManager *contact_manager;
  };

  BvhContactManager::BvhContactManager(std::shared_ptr<ContactInterface> interface, std::size_t _overdecomposition)
    : ParallelContactManager(interface),
      m_world{_overdecomposition},
      m_nodes{&m_world.create_collision_object()},
      m_faces{&m_world.create_collision_object()}
  {
    m_world.set_narrowphase_functor< ContactEntity >( NarrowphaseFunc{ this } );
  }

  BvhContactManager::BvhContactManager(BvhContactManager &&) = default;

  BvhContactManager &
  BvhContactManager::operator=( BvhContactManager &&) = default;

  BvhContactManager::~BvhContactManager() = default;

  void BvhContactManager::ComputeParallelContactForce(int step, bool debug_output) {
    m_world.start_iteration();

    // Update collision objects, this will build the trees
    for ( auto &&node : contact_nodes_ )
      node.RecomputeKdop();
    for ( auto &&face : contact_faces_ )
      face.RecomputeKdop();

    m_nodes->set_entity_data(contact_nodes_);
    m_faces->set_entity_data(contact_faces_);

    m_faces->broadphase(*m_nodes);

    m_last_results.clear();
    m_faces->for_each_result< NarrowphaseResult >( [this]( const NarrowphaseResult &_res ){
      m_last_results.emplace_back( _res );
    } );

    m_world.finish_iteration();

    for ( auto &f : force_ )
      f = 0.0;

    // Update contact entities
    for ( auto &&r : m_last_results )
    {
      if ( r.node )
      {
        if ( r.local_index >= contact_nodes_.size() )
          std::cerr << "contact node index " << r.local_index << " is out of bounds (" << contact_nodes_.size() << ")\n";
        auto &node = contact_nodes_.at(r.local_index);
        node.set_contact_status( true );
        std::cout << "(" << m_rank << ") node " << r.local_index << " got force " << r.contact_force[0] << ", " << r.contact_force[1] << ", " << r.contact_force[2] << '\n';
        node.SetNodalContactForces( r.contact_force );
        node.ScatterForceToContactManagerForceVector(force_);
      } else {
        if (r.local_index >= contact_faces_.size())
          std::cerr << "contact face index " << r.local_index << " is out of bounds (" << contact_faces_.size()
                    << ")\n";
        auto &face = contact_faces_.at(r.local_index);
        face.set_contact_status( true );
        std::cout << "(" << m_rank << ") face " << r.local_index << " got force " << r.contact_force[0] << ", " << r.contact_force[1] << ", " << r.contact_force[2] << '\n';
        face.SetNodalContactForces( r.contact_force, r.bary );
        face.ScatterForceToContactManagerForceVector(force_);
      }
    }

    if ( !m_last_results.empty() )
      std::cout << m_last_results.size() << " total collisions on node " << m_rank << '\n';
  }

  namespace
  {

    void
    WriteContactNodesToVTKFile(const std::vector<ContactEntity> &nodes,
                               const std::string &prefix,
                               int step) {
      std::stringstream file_name_ss;
      file_name_ss << prefix;
      file_name_ss << std::setfill( '0' ) << std::setw( 5 ) << step;
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

      for (int i=0 ; i<nodes.size() ; i++) {
        vis_file << nodes[i].coord_1_x_ << " " << nodes[i].coord_1_y_ << " " << nodes[i].coord_1_z_ << std::endl;
        vtk_vertices[i] = i;
      }

      vis_file << "CELLS " << vtk_vertices.size() << " " << 2*vtk_vertices.size() << std::endl;
      for(const auto & vtk_vertex : vtk_vertices) {
        vis_file << "1 " << vtk_vertex << std::endl;
      }

      vis_file << "CELL_TYPES " << nodes.size() << std::endl;
      for(unsigned int i=0 ; i<vtk_vertices.size() ; ++i) {
        // cell type 1 is VTK_VERTEX
        vis_file << 1 << std::endl;
      }

      vis_file.close();
    }

    void
    WriteContactFacesToVTKFile(const std::vector<ContactEntity> &faces,
                                  const std::string &prefix,
                                  int step) {

      std::stringstream file_name_ss;
      file_name_ss << prefix;
      file_name_ss << std::setfill( '0' ) << std::setw( 5 ) << step;
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

      std::vector< std::vector<int> > vtk_triangles(faces.size());

      vis_file << "POINTS " << 3*faces.size() << " float" << std::endl;

      for (int i=0 ; i<faces.size() ; i++) {
        ContactEntity const & face = faces[i];
        vis_file << face.coord_1_x_ << " " << face.coord_1_y_ << " " << face.coord_1_z_ << std::endl;
        vis_file << face.coord_2_x_ << " " << face.coord_2_y_ << " " << face.coord_2_z_ << std::endl;
        vis_file << face.coord_3_x_ << " " << face.coord_3_y_ << " " << face.coord_3_z_ << std::endl;
        vtk_triangles[i].push_back(3*i);
        vtk_triangles[i].push_back(3*i + 1);
        vtk_triangles[i].push_back(3*i + 2);
      }

      vis_file << "CELLS " << vtk_triangles.size() << " " << 4*vtk_triangles.size() << std::endl;
      for(const auto & vtk_triangle : vtk_triangles) {
        vis_file << "3 " << vtk_triangle[0] << " " << vtk_triangle[1] << " " << vtk_triangle[2] << std::endl;
      }

      vis_file << "CELL_TYPES " << faces.size() << std::endl;
      for(unsigned int i=0 ; i<vtk_triangles.size() ; ++i) {
        // cell type 6 is VTK_TRIANGLE
        vis_file << 6 << std::endl;
      }

      vis_file.close();
    }

    void
    WriteContactEntitiesToVTKFile(const std::vector<ContactEntity> &faces,
                                  const std::vector<ContactEntity> &nodes,
                                  const std::string &prefix,
                                  int step) {

      std::stringstream file_name_ss;
      file_name_ss << prefix;
      file_name_ss << std::setfill( '0' ) << std::setw( 5 ) << step;
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
      std::vector< std::vector<int> > vtk_triangles(faces.size());

      vis_file << "POINTS " << nodes.size() + 3*faces.size() << " float" << std::endl;

      for (int i=0 ; i<nodes.size() ; i++) {
        vis_file << nodes[i].coord_1_x_ << " " << nodes[i].coord_1_y_ << " " << nodes[i].coord_1_z_ << std::endl;
        vtk_vertices[i] = i;
      }
      int offset = static_cast<int>(vtk_vertices.size());
      for (int i=0 ; i<faces.size() ; i++) {
        ContactEntity const & face = faces[i];
        vis_file << face.coord_1_x_ << " " << face.coord_1_y_ << " " << face.coord_1_z_ << std::endl;
        vis_file << face.coord_2_x_ << " " << face.coord_2_y_ << " " << face.coord_2_z_ << std::endl;
        vis_file << face.coord_3_x_ << " " << face.coord_3_y_ << " " << face.coord_3_z_ << std::endl;
        vtk_triangles[i].push_back(offset + 3*i);
        vtk_triangles[i].push_back(offset + 3*i + 1);
        vtk_triangles[i].push_back(offset + 3*i + 2);
      }

      vis_file << "CELLS " << vtk_vertices.size() + vtk_triangles.size() << " " << 2*vtk_vertices.size() + 4*vtk_triangles.size() << std::endl;
      for(const auto & vtk_vertex : vtk_vertices) {
        vis_file << "1 " << vtk_vertex << std::endl;
      }
      for(const auto & vtk_triangle : vtk_triangles) {
        vis_file << "3 " << vtk_triangle[0] << " " << vtk_triangle[1] << " " << vtk_triangle[2] << std::endl;
      }

      vis_file << "CELL_TYPES " << nodes.size() + faces.size() << std::endl;
      for(unsigned int i=0 ; i<vtk_vertices.size() ; ++i) {
        // cell type 1 is VTK_VERTEX
        vis_file << 1 << std::endl;
      }
      for(unsigned int i=0 ; i<vtk_triangles.size() ; ++i) {
        // cell type 6 is VTK_TRIANGLE
        vis_file << 6 << std::endl;
      }

      vis_file.close();
    }
  }
}
