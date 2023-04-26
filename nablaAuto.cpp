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

// autodiff
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

using namespace autodiff;

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
real old_L( Containers::Vector< t >& nabla_h,
        Containers::Vector< t >& nabla )
{
    real L = 0.0;
    for( int i = 0; i < nabla_h.getSize(); i++ )
    {
        L += l2Norm( nabla_h[ i ] - nabla[ i ] ) * l2Norm( nabla_h[ i ] - nabla[ i ] );
    }
    return L;
}

real L( const ArrayXreal& x )
{
    real L = 0.0;
    for( int i = 0; i < x.size(); i++ )
    {
        L += ( x(4*i) - x(4*i+2) ) * ( x(4*i) - x(4*i+2) )
           + ( x(4*i+1) - x(4*i+3) ) * ( x(4*i+1) - x(4*i+3) );
    }
    return L;
}

template< typename MeshConfig >
bool nablaAuto( Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
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
    }

    Containers::Vector< PointType, Devices::Host > faceCenters ( facesCount );
    for(int i = 0; i < facesCount; i++)
    {
        auto face = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( i );
        PointType center = getEntityCenter( mesh, face );
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

    using Eigen::VectorXd;
    ArrayXreal var( 4 * verticesCount ); // the input array x with verticesCount variables
    for( int i = 0; i < verticesCount; i++ ) // fill x
    {
        real aux = nabla_h[ i ][ 0 ];
        var( 4*i ) = aux;
        aux = nabla_h[ i ][ 1 ];
        var( 4*i +1 ) = aux;
        aux = nabla[ i ][ 0 ];
        var( 4*i +2 ) = aux;
        aux = nabla[ i ][ 1 ];
        var( 4*i +3 ) = aux;
    }
    std::cout << "input array:\n" << var << "\n";

    // computing L(mesh)
    double loss = L( var );
    

    // give the points vector value of respective vertices
    Containers::Vector< PointType > nabla_mesh( verticesCount );
    nabla_mesh = 0;

    auto kernel = [ &mesh, &nabla_mesh, &nabla_h, &nabla, &var ] ( GlobalIndexType i )  mutable
    {
        
    };
    mesh.template forAll< 0 >( kernel );


    real u;                                     // the output scalar u = L(x) evaluated together with gradient below

    VectorXd g = gradient(L, wrt(var), at(var), u); // evaluate the function value u and its gradient vector g = du/dx

    std::cout << "u = " << u << "\n";      // print the evaluated output u
    std::cout << verticesCount << " vertices" << "\n";
    std::cout << "nabla_mesh of size " << g.size() << " = \n" << g << "\n";

    for( int i = 0; i <= g.size(); i += 2 )
    {
        nabla_mesh[ i/2 ][ 0 ] = g( i );
        nabla_mesh[ i/2 ][ 1 ] = g( i+1 );
    }


    Containers::Array< double > nabla_arr( 3 * verticesCount );
    nabla_arr = 0;
    for( int i = 0; i < 3 * verticesCount; i += 3 )
    {
        nabla_arr[ i ] = nabla_mesh[ i / 3 ][ 0 ];
        nabla_arr[ i + 1 ] = nabla_mesh[ i / 3 ][ 1 ];
        nabla_arr[ i + 2 ] = 0; // redundant
    }

    auto descent = [ &mesh, &nabla_mesh ] ( GlobalIndexType i ) mutable
    {
        mesh.getPoints()[ i ] -= 1e-4 * nabla_mesh[ i ]; // TODO change parameter
    };
    mesh.template forInterior< 0 >( descent );

    // writing the computed gradient into a new mesh
    using VTKWriter = Meshes::Writers::VTKWriter< MeshType >;
    std::ofstream out = std::ofstream( "defGrads.vtk" );
    VTKWriter writer = VTKWriter( out );
    writer.template writeEntities< MeshType::getMeshDimension() >( mesh );
    writer.writePointData( nabla_arr, "meshGrads", 3 );

    nabla_h = 0;
    nabla = 0;
    mesh.template forAll< MeshType::getMeshDimension() >( get_nabla_h );
    double new_loss = old_L< PointType >( nabla_h, nabla );

    double improvement = loss - new_loss;
    
    std::cout << "L(initial mesh) = " << loss << "\n"
              << "L(updated mesh) = " << new_loss << "\n"
              << "\nimprovement of L: " << improvement << "\n\n";

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
            return nablaAuto(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
