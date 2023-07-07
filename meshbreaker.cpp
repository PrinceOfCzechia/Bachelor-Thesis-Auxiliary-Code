#include <TNL/Meshes/Mesh.h>
#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Meshes/Geometry/getEntityMeasure.h>
#include <TNL/Meshes/TypeResolver/resolveMeshType.h>
#include <TNL/Meshes/Writers/VTKWriter.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>


using namespace TNL;
using namespace TNL::Meshes;

struct MyConfigTag
{};

namespace TNL {
namespace Meshes {
namespace BuildConfigTags {

/****
 * Turn off all grids.
 */
template<> struct GridRealTag< MyConfigTag, float > { static constexpr bool enabled = false; };
template<> struct GridRealTag< MyConfigTag, double > { static constexpr bool enabled = false; };
template<> struct GridRealTag< MyConfigTag, long double > { static constexpr bool enabled = false; };

template<> struct MeshCellTopologyTag< MyConfigTag, Topologies::Triangle >{ static constexpr bool enabled = true; };

// Meshes are enabled only for the world dimension equal to the cell dimension.
template< typename CellTopology, int WorldDimension >
struct MeshSpaceDimensionTag< MyConfigTag, CellTopology, WorldDimension >
{ static constexpr bool enabled = WorldDimension == CellTopology::dimension; };

// Meshes are enabled only for types explicitly listed below.
template<> struct MeshRealTag< MyConfigTag, float >{ static constexpr bool enabled = true; };
template<> struct MeshRealTag< MyConfigTag, double >{ static constexpr bool enabled = true; };
template<> struct MeshGlobalIndexTag< MyConfigTag, int >{ static constexpr bool enabled = true; };
template<> struct MeshGlobalIndexTag< MyConfigTag, long int >{ static constexpr bool enabled = true; };
template<> struct MeshLocalIndexTag< MyConfigTag, short int >{ static constexpr bool enabled = true; };

}  // namespace BuildConfigTags
}  // namespace Meshes
}  // namespace TNL

const double P = 2.0e-1;

template< typename MeshConfig >
bool breakMesh( Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
{
    using MeshType = Mesh< MeshConfig, Devices::Host >;
    using RealType = typename MeshType::RealType;
    using GlobalIndexType = typename MeshType::GlobalIndexType;
    using LocalIndexType = typename MeshType::LocalIndexType;
    using VectorType = TNL::Containers::Vector< RealType, TNL::Devices::Host, GlobalIndexType >;
    using PointType = typename MeshTraits< MeshConfig >::PointType;

    const int verticesCount = mesh.template getEntitiesCount< 0 >();
    const int facesCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() - 1 >();
    const int cellsCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() >();

    double minL = getEntityMeasure( mesh, mesh.template getEntity< MeshType::getMeshDimension() - 1 >( 0 ) );
    auto getShortestFace = [ &mesh, &minL ] ( GlobalIndexType i ) mutable
    {
        auto measure = getEntityMeasure( mesh, mesh.template getEntity< MeshType::getMeshDimension() - 1 >( 0 ) );
        if( measure < minL ) minL = measure;
    };
    mesh.template forAll< MeshType::getMeshDimension() - 1 >( getShortestFace );

    auto breakMesh = [ &mesh, &minL ] ( GlobalIndexType i )  mutable
    {
        mesh.getPoints()[ i ][ 0 ] += P * minL * ( ((double) rand() / (RAND_MAX)) * 2 - 1 );
        mesh.getPoints()[ i ][ 1 ] += P * minL * ( ((double) rand() / (RAND_MAX)) * 2 - 1 );
    };
    mesh.template forInterior< 0 >( breakMesh );

    // writing the perturbated mesh into a new file
    using VTKWriter = Meshes::Writers::VTKWriter< MeshType >;
    std::ofstream out = std::ofstream( "broken.vtk" );
    VTKWriter writer = VTKWriter( out );
    writer.template writeEntities< MeshType::getMeshDimension() >( mesh );

    std::cout << "Mesh broken succesfully" << "\n";
    return true;
}

int main( int argc, char* argv[] )
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << argv[ 0 ] << " [mesh file adress]" << "\n";
        return EXIT_FAILURE;
    }

    bool result = true;

    for( int i = 1; i < argc; i++ )
    {
        const std::string fileName = argv[ i ];
        auto wrapper = [&]( auto& reader, auto&& mesh ) -> bool
        {
            return breakMesh(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
