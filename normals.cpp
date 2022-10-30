#include <TNL/Meshes/Mesh.h>
#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Meshes/Geometry/getEntityMeasure.h>
#include <TNL/Meshes/TypeResolver/resolveMeshType.h>

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

    // normala: souradnice vrcholu na stene, pravdepodobne getCenter()
    // viz clanek v mailu
    // na priste implementovat a otestovat
    // check python for vtk meshes
    // deformace pro kontrolu spravnosti vypoctu, provest primo v programu (rotace >>> zkoseni)
    // pomoci fce getPoints
    // PointType StaticVector, PointArrayType Array< PointType >, kdyz to jde, pouzit typ auto
    // auto Points = mesh.getPoints()[ pointIdx ];
    // zobrazit v ParaView ?? pomoci VTKwriter.h ?? spis ne
    // kdyz hotovo, pokracovat ve schematu

template< typename V >
V normalize(V v)
{
    V u = (1.0 / l2Norm( v )) * v;
    return u;
}

template< typename MeshConfig >
bool getNormals( const Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
{
    using MeshType = Mesh< MeshConfig, Devices::Host >;
    using RealType = typename MeshType::RealType;
    using GlobalIndexType = typename MeshType::GlobalIndexType;
    using LocalIndexType = typename MeshType::LocalIndexType;
    using VectorType = TNL::Containers::Vector< RealType, TNL::Devices::Host, GlobalIndexType >;
    using PointType = typename MeshTraits< MeshConfig >::PointType;

    const auto verticesCount = mesh.template getEntitiesCount< 0 >();
    const auto facesCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() - 1 >();
    const auto cellsCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() >();

    Containers::StaticVector< 8, PointType > verticesVector;
    for(int i = 0; i < cellsCount; i++)
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType vertex0 = mesh.getPoint( cell.template getSubentityIndex< 0 > ( 0 ));
        PointType vertex1 = mesh.getPoint( cell.template getSubentityIndex< 0 > ( 1 ));
        PointType vertex2 = mesh.getPoint( cell.template getSubentityIndex< 0 > ( 2 ));
        // note:
        // auto vertex = mesh.getPoint ( cell.template getSubentityIndex< 0 > ( 2 ));
        // works just as well
        std::cout << "\nCell " << i << ", vertex 0: " << vertex0;
        std::cout << "\nCell " << i << ", vertex 1: " << vertex1;
        std::cout << "\nCell " << i << ", vertex 2: " << vertex2;
    }

    Containers::Vector< Containers::Vector< PointType, TNL::Devices::Host > > normals ( cellsCount );
    for(int i = 0; i < cellsCount; i++)
    {
        normals[ i ] = Containers::Vector< PointType, TNL::Devices::Host > ( 3 );
    }

    for(int i = 0; i < cellsCount; i++)
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        for(int j = 0; j < 3; j++)
        {
            const auto faceIdx = cell.template getSubentityIndex< MeshType::getMeshDimension() - 1 >( j );
            std::cout << "\nCell " << i << ", local face index: " << j << ", global face index: " << faceIdx /*<< ", face: "
                      << mesh.template getEntity< MeshType::getMeshDimension() - 1 >( faceIdx )*/ << std::endl;
            PointType faceVector = mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2)%3 )) 
                                 - mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1)%3 ));
            std::cout << "Face vector: " << faceVector;
            PointType outwardNormal = { faceVector[1], -faceVector[0] };
            normals[ i ][ j ] = normalize< PointType >(outwardNormal);
            std::cout << ", outward normal vector: " << outwardNormal;
        }
    }

    for(int i = 0; i < verticesCount; i++)
    {
        PointType vertex = mesh.template getPoint(i);
        verticesVector[i] = vertex;
    }

    std::cout << "\n" << verticesVector << std::endl;
    std::cout << "\nOutward normal vectors (row ~ cellIdx, col ~ local faceIdx): " << std::endl;
    for(int i = 0; i < cellsCount; i++)
    {
        std::cout << normals[i] << std::endl;
    }

    return true;
}

template< typename V >
V rotate2D( V v, const double T)
{
    V u;
    auto v0 = v[0];
    auto v1 = v[1];

    u[0] = v0 * cos(T) - v1 * sin(T);
    u[1] = v1 * cos(T) + v0 * sin(T);

    return u;
}

template< typename MeshConfig >
bool rotateNormals( const Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
{
    const double T = M_PI/2;

    using MeshType = Mesh< MeshConfig, Devices::Host >;
    using RealType = typename MeshType::RealType;
    using GlobalIndexType = typename MeshType::GlobalIndexType;
    using LocalIndexType = typename MeshType::LocalIndexType;
    using VectorType = TNL::Containers::Vector< RealType, TNL::Devices::Host, GlobalIndexType >;
    using PointType = typename MeshTraits< MeshConfig >::PointType;

    const auto verticesCount = mesh.template getEntitiesCount< 0 >();
    const auto facesCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() - 1 >();
    const auto cellsCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() >();

    std::cout << "Rotated mesh (row n ~ before, row n+1 ~ after rotation):\n";

    Containers::StaticVector< 8, PointType > verticesVector;
    for(int i = 0; i < cellsCount; i++)
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType vertex0 = mesh.getPoint( cell.template getSubentityIndex< 0 > ( 0 ));
        PointType vertex1 = mesh.getPoint( cell.template getSubentityIndex< 0 > ( 1 ));
        PointType vertex2 = mesh.getPoint( cell.template getSubentityIndex< 0 > ( 2 ));
        
        std::cout << "\nCell " << i << ", vertex 0: " << vertex0;
        std::cout << "\tCell " << i << ", vertex 1: " << vertex1;
        std::cout << "\tCell " << i << ", vertex 2: " << vertex2;

        PointType u0 = rotate2D< PointType >(vertex0, T);
        PointType u1 = rotate2D< PointType >(vertex1, T);
        PointType u2 = rotate2D< PointType >(vertex2, T);

        std::cout << "\nCell " << i << ", vertex 0: " << u0;
        std::cout << "\tCell " << i << ", vertex 1: " << u1;
        std::cout << "\tCell " << i << ", vertex 2: " << u2;
    }
    std::cout << "\n";

    Containers::Vector< Containers::Vector< PointType, TNL::Devices::Host > > normals ( cellsCount );
    for(int i = 0; i < cellsCount; i++)
    {
        normals[ i ] = Containers::Vector< PointType, TNL::Devices::Host > ( 3 );
    }

    for(int i = 0; i < cellsCount; i++)
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        for(int j = 0; j < 3; j++)
        {
            const auto faceIdx = cell.template getSubentityIndex< MeshType::getMeshDimension() - 1 >( j );
            std::cout << "\nCell " << i << ", local face index: " << j << ", global face index: " << faceIdx << std::endl;
            PointType faceVector = rotate2D< PointType >( mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2)%3 )), T ) 
                                 - rotate2D< PointType >( mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1)%3 )), T );
            std::cout << "Face vector: " << faceVector;
            PointType outwardNormal = { faceVector[1], -faceVector[0] };
            normals[ i ][ j ] = normalize< PointType >(outwardNormal);
            std::cout << ", outward normal vector: " << normals[i][j];
        }
    }

    for(int i = 0; i < verticesCount; i++)
    {
        PointType vertex = rotate2D< PointType >( mesh.template getPoint(i), T );
        verticesVector[i] = vertex;
    }

    std::cout << "\nMesh vertices:\n" << verticesVector << std::endl;
    std::cout << "\nOutward normal vectors (row ~ cellIdx, col ~ local faceIdx): " << std::endl;
    for(int i = 0; i < cellsCount; i++)
    {
        std::cout << normals[i] << std::endl;
    }

    return true;
}

int
main( int argc, char* argv[] )
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << argv[ 0 ] << " [mesh file adress]" << std::endl;
        return EXIT_FAILURE;
    }

    bool result = true;

    for( int i = 1; i < argc; i++ )
    {
        const std::string fileName = argv[ i ];
        auto wrapper = [&]( auto& reader, auto&& mesh ) -> bool
        {
            // return iterateMesh( mesh, "" );
            return rotateNormals(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
