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
#include <autodiff/forward/dual.hpp>

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

    // return sin(x);
    return std::exp( - x*x - y*y );
    // return  x*x*sin(x*x + y*y) + y*y*sin(x*x + y*y) ;
}

real f_eigen( real x, real y )
{
    return real( exp( - x*x - y*y ) );
}

// manually computed gradient of the function above
template< typename V >
V angrad( V v )
{
    double x = v[0];
    double y = v[1];

    V grad = { 0, 0 };
    // grad[0] = cos(x);
    grad[ 0 ] = -2 * x * std::exp( -x*x - y*y );
    grad[ 1 ] = -2 * y * std::exp( -x*x - y*y );
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
double old_L( Containers::Vector< t >& nabla_h,
            Containers::Vector< t >& nabla )
{
    double L = 0.0;
    for( int i = 0; i < nabla_h.getSize(); i++ )
    {
        L += l2Norm( nabla_h[ i ] - nabla[ i ] ) * l2Norm( nabla_h[ i ] - nabla[ i ] );
    }
    return L;
}

real L( const ArrayXreal& x )
{
    real L = 0.0;
    for( int i = 0; i < x.size(); i += 2 )
    {
        L += x( i ) * x( i ) + x( i+1 ) * x( i+1 );
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

    using Eigen::VectorXd;
    using Eigen::ArrayXd;

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

    ArrayXreal meshPoints( 2 * verticesCount );
    auto copyPoints = [ &mesh, &meshPoints ] ( GlobalIndexType i ) mutable
    {
        meshPoints[ 2*i ] = mesh.getPoints()[ i ][ 0 ];
        meshPoints[ 2*i + 1] = mesh.getPoints()[ i ][ 1 ];
    };
    mesh.template forAll< 0 >( copyPoints );

    ArrayXreal nabla_h_eigen( 2 * verticesCount );
    ArrayXreal nabla_eigen( 2 * verticesCount );

    auto get_nabla_h = [ &mesh, &nabla_h, &nabla, &meshPoints, &nabla_h_eigen, &nabla_eigen ] ( GlobalIndexType i )  mutable
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType sum = { 0, 0 };
        real sum0 = 0;
        real sum1 = 0;
        real cellMeasure = 0;
        for( int j = 0; j < 3; j++ )
        // TODO: to reals
        {
            const auto faceIdx = cell.template getSubentityIndex< MeshType::getMeshDimension() - 1 >( j );
            const auto sigma = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( faceIdx );
            const int idx0 = cell.template getSubentityIndex< 0 >( (j+1)%3 );
            const int idx1 = cell.template getSubentityIndex< 0 >( (j+2)%3 );
            const int idx2 = cell.template getSubentityIndex< 0 >( j%3 );

            real faceVector0 = meshPoints[ 2*idx1 ] - meshPoints[ 2*idx0 ];
            real faceVector1 = meshPoints[ 2*idx1 + 1 ] - meshPoints[ 2*idx0 + 1 ];
            real faceNorm = sqrt( faceVector0 * faceVector0 + faceVector1 * faceVector1 );

            real otherFace0 = meshPoints[ 2*idx1 ] - meshPoints[ 2*idx2 ];
            real otherFace1 = meshPoints[ 2*idx1 + 1 ] - meshPoints[ 2*idx2 + 1 ];

            cellMeasure += abs( faceVector0 * otherFace1 - faceVector1 * otherFace0 ) / 6;

            real outwardNormal0 = faceVector1/faceNorm;
            real outwardNormal1 = -faceVector0/faceNorm;

            real faceCenter0 = ( meshPoints[ 2*idx1 ] + meshPoints[ 2*idx0 ] ) / 2;
            real faceCenter1 = ( meshPoints[ 2*idx1 + 1 ] + meshPoints[ 2*idx0 + 1 ] ) / 2;

            real f_sigma_eigen = f_eigen( faceCenter0, faceCenter1 );

            sum0 += faceNorm * f_sigma_eigen * outwardNormal0;
            sum1 += faceNorm * f_sigma_eigen * outwardNormal1;

            PointType faceVector = mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2) % 3 ) ) 
                                 - mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1) % 3 ) );
            
            PointType outwardNormal = normalize< PointType >( { faceVector[ 1 ], -faceVector[ 0 ] } );

            PointType x_sigma = getEntityCenter( mesh, sigma ); // vlastni fce

            double f_sigma = f< PointType >( x_sigma ); // vracet dual
            sum += getEntityMeasure( mesh, sigma ) * f_sigma * outwardNormal; // vlastni getEntityMeasure
        }

        sum0 /= cellMeasure;
        sum1 /= cellMeasure;

        PointType grad_h = ( 1.0 / getEntityMeasure( mesh, cell ) ) * sum;
        PointType grad = angrad< PointType >( getEntityCenter( mesh, cell ) );

        for( int j = 0; j < 3; j++ )
        {
            int globalPointIdx = cell.template getSubentityIndex< 0 >( j );
            nabla[ globalPointIdx ] += grad;
            nabla_h[ globalPointIdx ] += grad_h;
            nabla_h_eigen[ 2*globalPointIdx ] += sum0;
            nabla_h_eigen[ 2*globalPointIdx + 1 ] += sum1;
            nabla_eigen[ 2*globalPointIdx ] += grad[ 0 ];
            nabla_eigen[ 2*globalPointIdx ] += grad[ 1 ];
        }
    };
    mesh.template forAll< MeshType::getMeshDimension() >( get_nabla_h );

    // std::cout << "TNL\n" << nabla_h << "\neigen\n" << nabla_h_eigen << std::endl;
    
    ArrayXreal var( 2 * verticesCount ); // the input array var with verticesCount variables
    for( int i = 0; i < 2*verticesCount; i++ ) // fill var
    {
        real aux = nabla_h_eigen[ i ] - nabla_eigen[ i ];
        var( i ) = aux;
    }

    // computing L(initial_mesh)
    double loss = L( var );
        

    Containers::Vector< PointType > nabla_mesh( verticesCount );
    nabla_mesh = 0;

    real u; // the output scalar u = L(x) evaluated together with gradient below

    Eigen::VectorXd g = gradient( L, wrt(var), at(var), u ); // evaluate the function value u and its gradient vector g

    for( int i = 0; i < g.size(); i += 2 )
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
    }

    double lengthConst = TNL::l2Norm( nabla_mesh[ 0 ] );
    auto getLongestGrad = [ &mesh, &nabla_mesh, &lengthConst ] ( GlobalIndexType i ) mutable
    {
        double N = TNL::l2Norm( nabla_mesh[ i ] );
        if( N > lengthConst ) lengthConst = N;
    };
    mesh.template forAll< 0 >( getLongestGrad );

    auto descent = [ &mesh, &nabla_mesh, &lengthConst ] ( GlobalIndexType i ) mutable
    {
        mesh.getPoints()[ i ] -= 1e-4 * ( 1/lengthConst ) * nabla_mesh[ i ]; // TODO change parameter
    };
    mesh.template forAll< 0 >( descent );

    // writing the computed gradient into a new mesh
    using VTKWriter = Meshes::Writers::VTKWriter< MeshType >;
    std::ofstream out = std::ofstream( "autoGrads.vtk" );
    VTKWriter writer = VTKWriter( out );
    writer.template writeEntities< MeshType::getMeshDimension() >( mesh );
        

    nabla_h = 0;
    nabla = 0;
    nabla_h_eigen = 0;
    nabla_eigen = 0;
    mesh.template forAll< MeshType::getMeshDimension() >( get_nabla_h );
    for( int i = 0; i < 2*verticesCount; i++ ) // fill var
    {
        real aux = nabla_h_eigen[ i ] - nabla_eigen[ i ];
        var( i ) = aux;
    }
    double new_loss = L( var );
    // double new_loss = old_L< PointType >( nabla_h, nabla );

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
