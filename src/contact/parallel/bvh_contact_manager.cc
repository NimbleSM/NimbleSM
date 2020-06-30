#include "bvh_contact_manager.h"
#include <vt/transport.h>
#include <fstream>

namespace nimble {

  struct NarrowphaseFunc {
    bvh::narrowphase_result operator()(const bvh::broadphase_collision< ContactEntity > &_a,
        const bvh::broadphase_collision< ContactEntity > &_b )
    {
      auto res = bvh::typed_narrowphase_result< NarrowphaseResult >();
      auto tree = build_snapshot_tree_top_down( _b.elements );

      for ( auto &&el : _a.elements )
        query_tree( tree, el, [&_a, &_b, &res]( std::size_t _ai, std::size_t _bi ) {
          auto &entry = res.emplace_back();
          entry.first_global_id = _ai;
          entry.second_global_id = _bi;
        } );

      return res;
    }
  };

  BvhContactManager::BvhContactManager(std::shared_ptr<ContactInterface> interface, std::size_t _overdecomposition)
    : ParallelContactManager(interface),
      m_world{_overdecomposition},
      m_nodes{&m_world.create_collision_object()},
      m_faces{&m_world.create_collision_object()}
  {
    m_world.set_narrowphase_functor< ContactEntity >( NarrowphaseFunc{} );
  }

  BvhContactManager::BvhContactManager(BvhContactManager &&) noexcept = default;

  BvhContactManager &
  BvhContactManager::operator=( BvhContactManager &&) noexcept = default;

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

    m_nodes->broadphase(*m_faces);

    m_last_results.clear();
    m_nodes->for_each_result< NarrowphaseResult >( [this]( const NarrowphaseResult &_res ){
      m_last_results.emplace_back( _res );
    } );

    m_world.finish_iteration();
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

  void
  BvhContactManager::ContactVisualizationWriteStep( double time )
  {
    std::vector<ContactEntity> colliding_faces;
    std::vector<ContactEntity> noncolliding_faces;
    for ( auto &&face : contact_faces_ )
    {
      if ( std::count_if( m_last_results.begin(), m_last_results.end(),
                          [&face]( NarrowphaseResult &res ){ return res.second_global_id == face.contact_entity_global_id(); } ) )
      {
        colliding_faces.push_back(face);
      } else {
        noncolliding_faces.push_back(face);
      }
    }

    std::vector<ContactEntity> colliding_nodes;
    std::vector<ContactEntity> noncolliding_nodes;
    for ( auto &&node : contact_nodes_ )
    {
      if ( std::count_if( m_last_results.begin(), m_last_results.end(),
                          [&node]( NarrowphaseResult &res ){ return res.first_global_id == node.contact_entity_global_id(); } ) )
      {
        colliding_nodes.push_back(node);
      } else {
        noncolliding_nodes.push_back(node);
      }
    }

    static int step = 0;

    std::stringstream node_colliding_out_name;
    node_colliding_out_name << "nodes_colliding." << Rank() << '.';
    std::stringstream face_colliding_out_name;
    face_colliding_out_name << "faces_colliding." << Rank() << '.';

    std::stringstream node_noncolliding_out_name;
    node_noncolliding_out_name << "nodes_noncolliding." << Rank() << '.';
    std::stringstream face_noncolliding_out_name;
    face_noncolliding_out_name << "faces_noncolliding." << Rank() << '.';

    WriteContactNodesToVTKFile(colliding_nodes, node_colliding_out_name.str(), step);
    WriteContactFacesToVTKFile(colliding_faces, face_colliding_out_name.str(), step);
    WriteContactNodesToVTKFile(noncolliding_nodes, node_noncolliding_out_name.str(), step);
    WriteContactFacesToVTKFile(noncolliding_faces, face_noncolliding_out_name.str(), step);

    ++step;
  }

}