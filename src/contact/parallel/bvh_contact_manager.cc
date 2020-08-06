#include "bvh_contact_manager.h"
#include <vt/transport.h>
#include <fstream>
#include <iomanip>

namespace nimble {

  struct NarrowphaseFunc {
    explicit NarrowphaseFunc( BvhContactManager *cm )
      : contact_manager{cm}
    {}

    bvh::narrowphase_result operator()(const bvh::broadphase_collision< ContactEntity > &_a,
        const bvh::broadphase_collision< ContactEntity > &_b )
    {
      auto res = bvh::typed_narrowphase_result< NarrowphaseResult >();
      auto tree = build_snapshot_tree_top_down( _a.elements );

      for ( auto &&elb : _b.elements ) {
        query_tree_local(tree, elb, [&_a, &_b, &elb, &res, this](std::size_t _i) {
          const auto &face = _a.elements[_i];
          const auto &node = elb;
          ContactManager::PROJECTION_TYPE proj_type;
          ContactEntity::vertex closest_point;
          double gap;
          double normal[3];
          contact_manager->SimpleClosestPointProjectionSingle(node, face,
                                                              &proj_type, &closest_point, gap, normal);
          if ( proj_type != ContactManager::PROJECTION_TYPE::UNKNOWN ) {
            if ( gap <= 0.15 * ( node.char_len_ + face.char_len_ ) ) {
              auto &entry = res.emplace_back();
              entry.first_global_id = face.contact_entity_global_id();
              entry.second_global_id = node.contact_entity_global_id();
              entry.gap = gap;
            }
          };
        } );
      }

      return res;
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

    // Update contact entities
    std::size_t fcount = 0;
    for ( auto &&face : contact_faces_ )
    {
      bool hit = std::count_if( m_last_results.begin(), m_last_results.end(),
                                [&face]( NarrowphaseResult &res ){ return res.first_global_id == face.contact_entity_global_id(); } ) > 0;
      face.set_contact_status( hit );
      if ( hit )
        ++fcount;
    }

    for ( auto &&node : contact_nodes_ )
    {
      node.set_contact_status( std::count_if( m_last_results.begin(), m_last_results.end(),
                                              [&node]( NarrowphaseResult &res ){ return res.second_global_id == node.contact_entity_global_id(); } ) > 0 );
    }

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
