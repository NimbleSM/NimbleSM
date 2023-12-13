#include <Kokkos_Core_fwd.hpp>
#include <stdexcept>
#include <tpl_include/CLI11/CLI11.hpp>
#include <Ionit_Initializer.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Region.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_FaceGenerator.h>
#include <Ioss_NodeBlock.h>
#include <filesystem>
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include "boundary.hpp"

int
main( int argc, char **argv )
{
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
    Ioss::Init::Initializer io;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CLI::App app("obtain load-balanced surface elements from a mesh", "preskin");

    std::string input_mesh_path;

    app.add_option("input", input_mesh_path, "the input mesh");

    CLI11_PARSE(app, argc, argv);

    std::cout << "reading from file \"" << input_mesh_path << "\"\n";

    Ioss::PropertyManager property_manager;

    auto *db = Ioss::IOFactory::create("exodus", input_mesh_path, Ioss::READ_RESTART, Ioss::ParallelUtils::comm_world(), property_manager);

    std::string err_msg = "<unknown error>";
    if ( !db || !db->ok( false, &err_msg ) )
      throw std::runtime_error( err_msg );

    Ioss::Region main_region( db );

    if ( rank == 0 )
    {
      main_region.output_summary(std::cout, false);
    }

    const auto &element_blocks = main_region.get_element_blocks();

    std::cout << "\n\n";
    Ioss::FaceGenerator face_generator(main_region);

    // Get all the faces in the mesh
    face_generator.generate_faces(0, true);

    std::unordered_map<std::string, std::vector<Ioss::Face>> boundaries;

    for ( const auto &eb : element_blocks )
    {
      std::cout << "processing block " << eb->name() << '\n';
      const auto& faces = face_generator.faces(eb->name());
      auto &boundary = boundaries[eb->name()];
      boundary.reserve(faces.size());
      for (auto&& face : faces)
      {
        // Faces on a boundary have only one element associated with them
        if ( face.elementCount_ == 1 )
        {
          boundary.emplace_back(face);
        }
      }
      std::cout << "found " << boundary.size() << " boundary faces\n";
    }

    Kokkos::View<double *> coords("coords", 0);
    Kokkos::View<int *> ids("ids", 0);
    auto num_boundary_nodes = ps::compute_boundary_nodes(main_region, boundaries, coords, ids);

    auto* out_db = Ioss::IOFactory::create(
        "exodus", "tmp.g", Ioss::WRITE_RESTART, Ioss::ParallelUtils::comm_world(), property_manager);
    err_msg = "<unknown error>";
    if ( !out_db || !out_db->ok( false, &err_msg ) )
      throw std::runtime_error( err_msg );

    Ioss::Region out_region(out_db, "skin");

    out_region.begin_mode(Ioss::STATE_DEFINE_MODEL);
    auto *out_node_block = new Ioss::NodeBlock(out_region.get_database(), "nodeblock_1", num_boundary_nodes, 3);
    out_region.add(out_node_block);
    for ( const auto &eb : element_blocks )
    {
      auto *topo = eb->topology()->face_type(0);
      if ( !topo )
        throw std::runtime_error( "invalid topology" );

      const std::string topo_str = topo->name() == "tri3" ? "trishell" : "shell";
      out_region.add(new Ioss::ElementBlock(out_region.get_database(), eb->name(), topo_str, boundaries.at(eb->name()).size()));
    }
    out_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    auto coords_h = Kokkos::create_mirror_view(coords);
    Kokkos::deep_copy(coords_h, coords);

    auto ids_h = Kokkos::create_mirror_view(ids);
    Kokkos::deep_copy(ids_h, ids);

    out_region.begin_mode(Ioss::STATE_MODEL);
    auto num_written = out_node_block->put_field_data("mesh_model_coordinates", coords_h.data(), coords_h.extent(0) * sizeof(double));
    std::cout << "wrote " << num_written << " output node coordinates\n";
    if ( num_written != num_boundary_nodes )
      throw std::runtime_error("could not write all node coord data!");
    num_written = out_node_block->put_field_data("ids", ids_h.data(), ids_h.extent(0) * sizeof(int));
    std::cout << "wrote " << num_written << " output node ids\n";
    if ( num_written != num_boundary_nodes )
      throw std::runtime_error("could not write all node id data!");
    for ( const auto &eb : element_blocks )
    {
      const auto &boundary = boundaries.at(eb->name());
      auto *out_block = out_region.get_element_block(eb->name());
      const auto num_connected_nodes = out_block->topology()->number_corner_nodes();

      std::vector<int> element_ids;
      element_ids.reserve(boundary.size());

      std::vector<int> connectivity;
      connectivity.reserve(num_connected_nodes * boundary.size());
      for (auto&& face : boundary) {
        element_ids.emplace_back(face.element[0]);

        for (std::size_t i = 0; i < num_connected_nodes; ++i)
        {
          connectivity.emplace_back(face.connectivity_[i]);
        }
      }

      out_block->put_field_data("ids", element_ids);
      out_block->put_field_data("connectivity", connectivity);
    }
    out_region.end_mode(Ioss::STATE_MODEL);

    if ( rank == 0 )
    {
      out_region.output_summary(std::cout, false);
    }
  }
  Kokkos::finalize();
  MPI_Finalize();
}
