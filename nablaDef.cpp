#include <TNL/Meshes/Mesh.h>
#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Meshes/Geometry/getEntityMeasure.h>
#include <TNL/Meshes/TypeResolver/resolveMeshType.h>
#include <TNL/Meshes/Writers/VTIWriter.h>
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

    return sin(x);
    //return std::exp( -x*x - y*y );
    // return  x*x*sin(x*x + y*y) + y*y*sin(x*x + y*y) ;
}

// manually computed gradient of the function above
template< typename V >
V angrad( V v )
{
    double x = v[0];
    double y = v[1];

    V grad = { 0, 0 };
    grad[0] = cos(x);
    // grad[ 0 ] = -2 * x * std::exp( -x*x - y*y );
    // grad[ 0 ] = -2 * y * std::exp( -x*x - y*y );
    // grad[0] = 2*x*sin(x*x+y*y) + 2*x*x*x*cos(x*x+y*y) + 2*x*y*cos(x*x+y*y);
    // grad[1] = 2*x*y*cos(x*x+y*y) + 2*y*y*y*cos(x*x+y*y) + 2*y*sin(x*x+y*y);
    return grad;
}

// face center
template< typename V >
V x_sigma( V v1, V v2)
{
    return 0.5 * ( v1 + v2 );
}

template< typename t >
double L( Containers::Vector< t >& nabla_h,
          Containers::Vector< t >& nabla )
{
    double L = 0.0;
    for( int i = 0; i < nabla_h.getSize(); i++ )
    {
        L += l2Norm( nabla_h[ i ] - nabla[ i ] ) * l2Norm( nabla_h[ i ] - nabla[ i ] );
    }
    return L;
}

const double EPSILON = 1e-12;

template< typename MeshConfig >
bool nablaDef( Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
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

    Containers::Vector< PointType, Devices::Host > cellCenters ( cellsCount );
    for(int i = 0; i < cellsCount; i++)
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType center = getEntityCenter( mesh, cell );
        // std::cout << "\nCell " << i << ", center: " << center;
    }

    Containers::Vector< PointType, Devices::Host > faceCenters ( facesCount );
    for(int i = 0; i < facesCount; i++)
    {
        auto face = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( i );
        PointType center = getEntityCenter( mesh, face );
        // std::cout << "\nFace " << i << ", center: " << center;
    }

    Containers::Vector< PointType > nabla_h ( verticesCount );
    Containers::Vector< PointType > nabla ( verticesCount );
    nabla_h = 0;
    nabla = 0;

    auto get_nabla_h = [ &mesh, &nabla_h, &nabla ] ( GlobalIndexType i )  mutable
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType sum = { 0, 0 };
        for( int j = 0; j < 3; j++ )
        {
            const auto faceIdx = cell.template getSubentityIndex< MeshType::getMeshDimension() - 1 >( j );
            const auto sigma = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( faceIdx );

            PointType faceVector = mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2) % 3 ) ) 
                                 - mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1) % 3 ) );

            PointType outwardNormal = normalize< PointType >( { faceVector[ 1 ], -faceVector[ 0 ] } );

            PointType x_sigma = getEntityCenter( mesh, sigma );

            double f_sigma = f< PointType >( x_sigma );
            sum += getEntityMeasure( mesh, sigma ) * f_sigma * outwardNormal;
        }

        PointType grad_h = ( 1.0 / getEntityMeasure( mesh, cell ) ) * sum;
        PointType grad = angrad< PointType >( getEntityCenter( mesh, cell ) );

        for( int j = 0; j < 3; j++ )
        {
        int globalPointIdx = cell.template getSubentityIndex< 0 >( j );
        nabla[ globalPointIdx ] += grad;
        nabla_h[ globalPointIdx ] += grad_h;
        }
    };
    mesh.template forAll< MeshType::getMeshDimension() >( get_nabla_h );

    Containers::Vector< PointType > nabla_h_eps ( verticesCount );
    Containers::Vector< PointType > nabla_eps ( verticesCount );
    nabla_h_eps = 0;
    nabla_eps = 0;

    auto get_nabla_h_eps = [ &mesh, &nabla_h_eps, &nabla_eps ] ( GlobalIndexType i )  mutable
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType sum = { 0, 0 };
        for( int j = 0; j < 3; j++ )
        {
            const auto faceIdx = cell.template getSubentityIndex< MeshType::getMeshDimension() - 1 >( j );
            const auto sigma = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( faceIdx );

            PointType faceVector = mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2) % 3 ) ) 
                                 - mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1) % 3 ) );

            PointType outwardNormal = normalize< PointType >( { faceVector[ 1 ], -faceVector[ 0 ] } );

            PointType x_sigma = getEntityCenter( mesh, sigma );

            double f_sigma = f< PointType >( x_sigma );
            sum += getEntityMeasure( mesh, sigma ) * f_sigma * outwardNormal;
        }

        PointType grad_h = ( 1.0 / getEntityMeasure( mesh, cell ) ) * sum;
        PointType grad = angrad< PointType >( getEntityCenter( mesh, cell ) );

        for( int j = 0; j < 3; j++ )
        {
        int globalPointIdx = cell.template getSubentityIndex< 0 >( j );
        nabla_eps[ globalPointIdx ] += grad;
        nabla_h_eps[ globalPointIdx ] += grad_h;
        }
    };


    // calculating L(mesh)
    double loss = L< PointType >( nabla_h, nabla );

    std::cout << "L = " << loss << "\n";

    // give the points vector value of respective vertices
    Containers::Vector< PointType > nabla_mesh( verticesCount );
    nabla_mesh = 0;

    auto kernel = [ &mesh, &nabla_mesh, &nabla_h, &nabla, &nabla_h_eps, &nabla_eps, get_nabla_h_eps ] ( GlobalIndexType i )  mutable
    {
        // first partial derivative
        PointType eps0 = { EPSILON, 0 };
        mesh.getPoints()[ i ] += eps0;
        mesh.template forAll< MeshType::getMeshDimension() >( get_nabla_h_eps );
        nabla_mesh[ i ][ 0 ] = ( L< PointType >( nabla_h_eps, nabla_eps) - L< PointType >( nabla_h, nabla ) ) / EPSILON;
        mesh.getPoints()[ i ] -= eps0;
        nabla_h_eps = 0;
        nabla_eps = 0;

        // second partial derivative
        PointType eps1 = { 0, EPSILON };
        mesh.getPoints()[ i ] += eps1;
        mesh.template forAll< MeshType::getMeshDimension() >( get_nabla_h_eps );
        nabla_mesh[ i ][ 1 ] = ( L< PointType >( nabla_h_eps, nabla_eps) - L< PointType >( nabla_h, nabla ) ) / EPSILON;
        mesh.getPoints()[ i ] -= eps1;
        nabla_h_eps = 0;
        nabla_eps = 0;
    };
    mesh.template forAll< 0 >( kernel );

    std::cout << nabla_mesh << "\n";

    std::cout << "OK" << "\n";
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
            return nablaDef(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
