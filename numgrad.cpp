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

template< typename V >
double f( V v )
{
    double x = v[0];
    double y = v[1];

    return x*y;
}

// manually computed gradient of the function above
template< typename V >
V angrad( V v )
{
    double x = v[0];
    double y = v[1];

    V grad = { 0, 0 };
    grad[0] = y;
    grad[1] = x;
    return grad;
}

template< typename MeshConfig >
bool grad( const Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
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

    // Containers::StaticVector< 8, PointType > cellCenters;
    Containers::Vector< PointType, Devices::Host > cellCenters ( cellsCount );
    for(int i = 0; i < cellsCount; i++)
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType center = getEntityCenter( mesh, cell );
        // std::cout << "\nCell " << i << ", center: " << center;
    }

    // Containers::StaticVector< 15, PointType > faceCenters;
    Containers::Vector< PointType, Devices::Host > faceCenters ( facesCount );
    for(int i = 0; i < facesCount; i++)
    {
        auto face = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( i );
        PointType center = getEntityCenter( mesh, face );
        // std::cout << "\nFace " << i << ", center: " << center;
    }

    std::cout << std::endl;

    Containers::Vector< PointType, TNL::Devices::Host > grads ( cellsCount );
    Containers::Vector< PointType, TNL::Devices::Host > analytical ( cellsCount );

    for( int i = 0; i < cellsCount; i++ )
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType sum = { 0, 0 };
        for( int j = 0; j < 3; j++ )
        {
            const auto faceIdx = cell.template getSubentityIndex< MeshType::getMeshDimension() - 1 >( j );
            const auto sigma = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( faceIdx );

            //std::cout << "\nCell " << i << ", local face index: " << j << ", global face index: " << faceIdx << std::endl;

            PointType faceVector = mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2)%3 )) 
                                 - mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1)%3 ));

            PointType outwardNormal = normalize< PointType >( { faceVector[1], -faceVector[0] } );
            //std::cout << "Outward normal vector n_sigma: " << outwardNormal;

            PointType x_sigma = getEntityCenter( mesh, sigma );
            //std::cout << ", face center: " << x_sigma << "\n";

            double f_sigma = f< PointType >( x_sigma );
            sum += getEntityMeasure( mesh, sigma ) * f_sigma * outwardNormal;
        }

        PointType cellCenter = getEntityCenter( mesh, cell);
        analytical[ i ] = angrad< PointType >( cellCenter );
        
        PointType grad = ( 1.0 / getEntityMeasure( mesh, cell ) ) * sum;
        grads[ i ] = grad;
    }

    // double compoundError = 0;
    std::cout << "\nNumerical approximation of gradient in cell centers (row ~ cell index): " << std::endl;
    for(int i = 0; i < cellsCount; i++)
    {
        PointType error = grads[ i ] - analytical[ i ];
        std::cout << "Numerical " << ": " << grads[ i ] << ", analytical: " << analytical[ i ]
                  << ", with error: " << l2Norm( error ) << std::endl;
        // compoundError += l2Norm( error );
    }

    // std::cout << "Compound error: " << compoundError / cellsCount << std::endl;

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
            return grad(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
