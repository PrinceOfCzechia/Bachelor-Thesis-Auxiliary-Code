#include <TNL/Meshes/Mesh.h>
#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Meshes/Geometry/getEntityMeasure.h>
#include <TNL/Meshes/TypeResolver/resolveMeshType.h>
//#include <TNL/Meshes/Writers/WTIWriter.h>
//#include <TNL/Meshes/Writers/WTKWriter.h>
#include <iostream>
#include <vector>

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

    return x*x*y*y;
    // return  x*x*sin(x*x + y*y) + y*y*sin(x*x + y*y) ;
}

// manually computed gradient of the function above
template< typename V >
V angrad( V v )
{
    double x = v[0];
    double y = v[1];

    V grad = { 0, 0 };
    grad[0] = 2*x*y*y;
    grad[1] = 2*x*x*y;
    // grad[0] = 2*x*sin(x*x+y*y) + 2*x*x*x*cos(x*x+y*y) + 2*x*y*cos(x*x+y*y);
    // grad[1] = 2*x*y*cos(x*x+y*y) + 2*y*y*y*cos(x*x+y*y) + 2*y*sin(x*x+y*y);
    return grad;
}

// TODO test
template< typename V >
V x_sigma( V v1, V v2)
{
    return v2 - v1;
}

// face measure
template< typename t >
t m1( t x00, t x01, t x10, t x11 )
{ return pow( ( (x10 - x00) * (x10 - x00) + (x11 - x01) * (x11 - x01) ), 0.5 ); }

// derivative of the face measure by the first component
template< typename t >
t dm10( t x00, t x01, t x10, t x11, bool direction )
{ return 1 / m1( x00, x01, x10, x11 ) * ( x10 - x00 ) * ( direction ? 1 : -1 ); }

// derivative of the face measure by the second component
template< typename t >
t dm11( t x00, t x01, t x10, t x11, bool direction )
{ return 1 / m1( x00, x01, x10, x11 ) * ( x11 - x01 ) * ( direction ? 1 : -1 ); }
// *** direction 1 ~ derivative by x10 or x11, direction 0 ~ derivative by x00 or x01 ***

// first component of a normal vector
template< typename t >
t n0( t x00, t x01, t x10, t x11 )
{ return 1 / m1( x00, x01, x10, x11 ) * ( x11 - x01 ); }

// second component
template< typename t >
t n1( t x00, t x01, t x10, t x11 )
{ return 1 / m1( x00, x01, x10, x11 ) * ( -1 ) * ( x10 - x00 ); }

// normal vector
template< typename t >
std::vector< t > n( t x00, t x01, t x10, t x11 )
{ return { n0( x00, x01, x10, x11 ), n1( x00, x01, x10, x11 ) }; }

// derivatives of the first component by the first component
template< typename t >
t d0n00( t x00, t x01, t x10, t x11 )
{ return 1 / ( pow( m1(x00, x01, x10, x11),3 ) ) * ( x10 - x00 ) * ( x11 - x01 ) + 1 / m1( x00, x01, x10, x11 ) * 0; }
template< typename t >
t d0n10( t x00, t x01, t x10, t x11 )
{ return -1 / ( pow( m1(x00, x01, x10, x11),3 ) ) * ( x10 - x00 ) * ( x11 - x01 ) + 1 / m1( x00, x01, x10, x11 ) * 0; }

// derivatives of the first component by the second component
template< typename t >
t d0n01( t x00, t x01, t x10, t x11 )
{ return 1 / ( pow( m1(x00, x01, x10, x11),3 ) ) * ( x10 - x00 ) * ( x00 - x10 ) + 1 / m1( x00, x01, x10, x11 ) * 1; }
template< typename t >
t d0n11( t x00, t x01, t x10, t x11 )
{ return -1 / ( pow( m1(x00, x01, x10, x11),3 ) ) * ( x10 - x00 ) * ( x00 - x10 ) + 1 / m1( x00, x01, x10, x11 ) * ( -1 ); }

// derivatives of the second component by the first component
template< typename t >
t d1n00( t x00, t x01, t x10, t x11 )
{ return 1 / ( pow( m1(x00, x01, x10, x11), 3 ) ) * ( x11 - x01 ) * ( x11 - x01 ) + 1 / m1( x00, x01, x10, x11 ) * ( -1 ); }
template< typename t >
t d1n10( t x00, t x01, t x10, t x11 )
{ return -1 / ( pow( m1(x00, x01, x10, x11), 3 ) ) * ( x11 - x01 ) * ( x11 - x01 ) + 1 / m1( x00, x01, x10, x11 ) * 1; }

// derivatives of the second component by the second component
template< typename t >
t d1n01( t x00, t x01, t x10, t x11 )
{ return 1 /( pow( m1(x00, x01, x10, x11), 3 ) ) * ( x11 - x01 ) * ( x00 - x10 ) + 1 / m1( x00, x01, x10, x11 ) * 0; }
template< typename t >
t d1n11( t x00, t x01, t x10, t x11 )
{ return -1/( pow( m1(x00, x01, x10, x11), 3 ) ) * ( x11 - x01 ) * ( x00 - x10 ) + 1 / m1( x00, x01, x10, x11 ) * 0; }

// auxiliary functions for the cell measure
template< typename t >
t S0( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return x20*x11 - x00*x11 - x20*x01 - x10*x21 + x10*x01 + x00*x21; }

template< typename t >
t S( t x00, t x01, t x10, t x11, t x20, t x21 )
{
    if(S0(x00, x01, x10, x11, x20, x21) > 0) return 1;
    else if (S0(x00, x01, x10, x11, x20, x21) < 0) return -1;
    else return 0;
}
// cell measure
template< typename t >
t m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return ( S0( x00, x01, x10, x11, x20, x21 ) >= 0 ? 1 : -1 ) * ( S0 ( x00, x01, x10, x11, x20, x21 ) ) / 2; }

template< typename t >
t inv_m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{
    return 1 / m2< t >( x00, x01, x10, x11, x20, x21 );
}

// derivatives of inverted cell measure
template< typename t >
t d00inv_m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return -1 / ( 2 * pow( m2( x00, x01, x10, x11, x20, x21 ), 2) ) * S ( x00, x01, x10, x11, x20, x21 ) * ( -x11 + x21 ); }

template< typename t >
t d10inv_m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return -1 / ( 2 * pow( m2( x00, x01, x10, x11, x20, x21 ), 2) ) * S( x00, x01, x10, x11, x20, x21 ) * ( -x21 + x01 ); }

template< typename t >
t d20inv_m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return -1 / ( 2 * pow( m2( x00, x01, x10, x11, x20, x21 ), 2) ) * S ( x00, x01, x10, x11, x20, x21 ) * ( -x01 + x11 ); }

template< typename t >
t d01inv_m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return -1 / ( 2 * pow( m2( x00, x01, x10, x11, x20, x21 ), 2) ) * S ( x00, x01, x10, x11, x20, x21 ) * ( -x20 + x10 ); }

template< typename t >
t d11inv_m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return -1 / ( 2 * pow( m2( x00, x01, x10, x11, x20, x21 ), 2) ) * S ( x00, x01, x10, x11, x20, x21 ) * ( -x00 + x20 ); }

template< typename t >
t d21inv_m2( t x00, t x01, t x10, t x11, t x20, t x21 )
{ return -1 / ( 2 * pow( m2( x00, x01, x10, x11, x20, x21 ), 2) ) * S ( x00, x01, x10, x11, x20, x21 ) * ( -x10 + x00 ); }

template< typename MeshConfig >
bool numScheme( const Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
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

    std::cout << std::endl;

    Containers::Vector< PointType, Devices::Host > nabla_h ( cellsCount );
    Containers::Vector< PointType, Devices::Host > nabla ( cellsCount );

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
            // std::cout << "Outward normal vector n_sigma: " << outwardNormal;

            PointType x_sigma = getEntityCenter( mesh, sigma );
            // std::cout << ", face center: " << x_sigma << "\n";

            double f_sigma = f< PointType >( x_sigma );
            sum += getEntityMeasure( mesh, sigma ) * f_sigma * outwardNormal;
        }

        PointType cellCenter = getEntityCenter( mesh, cell);
        nabla[ i ] = angrad< PointType >( cellCenter );
        
        PointType grad = ( 1.0 / getEntityMeasure( mesh, cell ) ) * sum;
        nabla_h[ i ] = grad;
        // std::cout << i << ": " << nabla_h[ i ] << std::endl;
    }

    // derivatives of nabla_h f(x^i), see LaTeX
    // derivatives by the first component
    Containers::Vector< PointType, Devices::Host > dk1_nabla_h ( cellsCount );
    dk1_nabla_h = 0;
    //derivatives by the second component
    Containers::Vector< PointType, Devices::Host > dk2_nabla_h ( cellsCount );
    dk2_nabla_h = 0;

    // point 0
    for( int j = 0; j < cellsCount; j++ )
    {
        
        Containers::Vector< PointType, Devices::Host > cp ( 3 ); // stands for cellPoints
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( j );
        for(int i = 0; i < 3; i++)
        {
            PointType point = mesh.getPoint( cell.template getSubentityIndex< 0 > ( i ));
            cp[ i ] = point;
        }

        // get global index of vertex 0 in cell j
        int globalPointIndex0 = cell.template getSubentityIndex< 0 >( 0 );

        // points of the chosen cell to iterate over
        double x00 = cp[ 0 ][ 0 ];
        double x01 = cp[ 0 ][ 1 ];
        double x10 = cp[ 1 ][ 0 ];
        double x11 = cp[ 1 ][ 1 ];
        double x20 = cp[ 2 ][ 0 ];
        double x21 = cp[ 2 ][ 1 ];

        // by the first component
        double d_nabla_h = 0;

        double m2 = inv_m2< double >( x00, x01, x10, x11, x20, x21 );

        double auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            // m1 in points (i+1)%3 and (i+2)%3 becuase of mesh indexing, v^0 is opposite to sigma^0
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1) ], cp[ (i+1)%3 ]  ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d00inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += dm10< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 0 ] * 0.5 /*1/2 from d x_sigma*/ *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d0n00< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIndex0 ][ 0 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << std::endl;

        // by the second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d01inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += dm11< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 0 ] * 0 /*0 from d x_sigma*/ *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d0n01< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIndex0 ][ 1 ] += d_nabla_h;
    }

    // point 1
    for( int j = 0; j < cellsCount; j++ )
    {
        
        Containers::Vector< PointType, Devices::Host > cp ( 3 ); // stands for cellPoints
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( j );
        for(int i = 0; i < 3; i++)
        {
            PointType point = mesh.getPoint( cell.template getSubentityIndex< 0 > ( i ));
            cp[ i ] = point;
        }

        int globalPointIndex1 = cell.template getSubentityIndex< 0 >( 1 );

        double x00 = cp[ 0 ][ 0 ];
        double x01 = cp[ 0 ][ 1 ];
        double x10 = cp[ 1 ][ 0 ];
        double x11 = cp[ 1 ][ 1 ];
        double x20 = cp[ 2 ][ 0 ];
        double x21 = cp[ 2 ][ 1 ];

        // by the first component

        double d_nabla_h = 0;

        double m2 = inv_m2< double >( x00, x01, x10, x11, x20, x21 );

        double auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d10inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += dm10< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 0 ] * 0.5 /*d x_sigma*/ *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d0n00< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIndex1 ][ 0 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << std::endl;

        // by the second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d01inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += dm11< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 0 ] * 0 /*0 from d x_sigma*/ *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d0n01< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIndex1 ][ 1 ] += d_nabla_h;
    }

    // point 2
    for( int j = 0; j < cellsCount; j++ )
    {
        
        Containers::Vector< PointType, Devices::Host > cp ( 3 ); // stands for cellPoints
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( j );
        for(int i = 0; i < 3; i++)
        {
            PointType point = mesh.getPoint( cell.template getSubentityIndex< 0 > ( i ));
            cp[ i ] = point;
        }

        int globalPointIndex2 = cell.template getSubentityIndex< 0 >( 2 );

        double x00 = cp[ 0 ][ 0 ];
        double x01 = cp[ 0 ][ 1 ];
        double x10 = cp[ 1 ][ 0 ];
        double x11 = cp[ 1 ][ 1 ];
        double x20 = cp[ 2 ][ 0 ];
        double x21 = cp[ 2 ][ 1 ];

        // by the first component

        double d_nabla_h = 0;

        double m2 = inv_m2< double >( x00, x01, x10, x11, x20, x21 );

        double auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d20inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += dm10< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 0 ] * 0.5 /*d x_sigma*/ *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d0n00< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIndex2 ][ 0 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << std::endl;

        // by the second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d01inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += dm11< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 0 ] * 0 /*0 from d x_sigma*/ *
                      n0< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d0n01< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIndex2 ][ 1 ] += d_nabla_h;
    }

    // second component
    // point 0
    for( int j = 0; j < cellsCount; j++ )
    {
        
        Containers::Vector< PointType, Devices::Host > cp ( 3 ); // stands for cellPoints
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( j );
        for(int i = 0; i < 3; i++)
        {
            PointType point = mesh.getPoint( cell.template getSubentityIndex< 0 > ( i ));
            cp[ i ] = point;
        }

        int globalPointIndex0 = cell.template getSubentityIndex< 0 >( 0 );

        double x00 = cp[ 0 ][ 0 ];
        double x01 = cp[ 0 ][ 1 ];
        double x10 = cp[ 1 ][ 0 ];
        double x11 = cp[ 1 ][ 1 ];
        double x20 = cp[ 2 ][ 0 ];
        double x21 = cp[ 2 ][ 1 ];

        // first component
        
        double d_nabla_h = 0;

        double m2 = inv_m2< double >( x00, x01, x10, x11, x20, x21 );

        double auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d20inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += dm10< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 1 ] * 0 /*0 from d x_sigma*/ *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d1n00< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIndex0 ][ 0 ] += d_nabla_h;

        // second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d01inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += dm11< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 1 ] * 0.5 /*1/2 from d x_sigma*/ *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d1n01< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIndex0 ][ 1 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << std::endl;
    }

    // point 1
    for( int j = 0; j < cellsCount; j++ )
    {
        
        Containers::Vector< PointType, Devices::Host > cp ( 3 ); // stands for cellPoints
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( j );
        for(int i = 0; i < 3; i++)
        {
            PointType point = mesh.getPoint( cell.template getSubentityIndex< 0 > ( i ));
            cp[ i ] = point;
        }

        int globalPointIndex1 = cell.template getSubentityIndex< 0 >( 1 );

        double x00 = cp[ 0 ][ 0 ];
        double x01 = cp[ 0 ][ 1 ];
        double x10 = cp[ 1 ][ 0 ];
        double x11 = cp[ 1 ][ 1 ];
        double x20 = cp[ 2 ][ 0 ];
        double x21 = cp[ 2 ][ 1 ];

        // first component
        
        double d_nabla_h = 0;

        double m2 = inv_m2< double >( x00, x01, x10, x11, x20, x21 );

        double auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d20inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += dm10< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 1 ] * 0 /*0 from d x_sigma*/ *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d1n00< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIndex1 ][ 0 ] += d_nabla_h;

        // second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d11inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += dm11< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 1 ] * 0.5 /*1/2 from d x_sigma*/ *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d1n01< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIndex1 ][ 1 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << std::endl;
    }

    // point 2
    for( int j = 0; j < cellsCount; j++ )
    {
        
        Containers::Vector< PointType, Devices::Host > cp ( 3 ); // stands for cellPoints
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( j );
        for(int i = 0; i < 3; i++)
        {
            PointType point = mesh.getPoint( cell.template getSubentityIndex< 0 > ( i ));
            cp[ i ] = point;
        }

        int globalPointIndex2 = cell.template getSubentityIndex< 0 >( 2 );

        double x00 = cp[ 0 ][ 0 ];
        double x01 = cp[ 0 ][ 1 ];
        double x10 = cp[ 1 ][ 0 ];
        double x11 = cp[ 1 ][ 1 ];
        double x20 = cp[ 2 ][ 0 ];
        double x21 = cp[ 2 ][ 1 ];

        // first component
        
        double d_nabla_h = 0;

        double m2 = inv_m2< double >( x00, x01, x10, x11, x20, x21 );

        double auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d20inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += dm10< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 1 ] * 0 /*0 from d x_sigma*/ *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d1n00< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIndex2 ][ 0 ] += d_nabla_h;

        // second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ i%3 ], cp[ (i+1)%3 ]  ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += d21inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += dm11< double >( cp[ i%3 ][ 0 ], cp[ (i+1)%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ], 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            // TODO autodiff instead of angrad
            auxSum += m1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ (i+1)%3 ], cp[ i%3 ]) )[ 0 ] * 0.5 /*1/2 from d x_sigma*/ *
                      n1< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            auxSum += m1< double >( cp[ (i+1)%3 ][ 0 ], cp[ (i+2)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ], cp[ (i+2)%3 ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ (i+1)%3 ], cp[ i%3 ] ) ) *
                      d1n01< double >( cp[ i%3 ][ 0 ], cp[ i%3 ][ 1 ], cp[ (i+1)%3 ][ 0 ], cp[ (i+1)%3 ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIndex2 ][ 1 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << std::endl;
    }

    std::cout << "dk1_nabla_h" << std::endl; 
    for( int i = 0; i < cellsCount; i++ )
    {
        std::cout << dk1_nabla_h[ i ];
    }
    std::cout << std::endl;

    std::cout << "dk2_nabla_h" << std::endl; 
    for( int i = 0; i < cellsCount; i++ )
    {
        std::cout << dk2_nabla_h[ i ];
    }
    std::cout << std::endl;

    Containers::Vector< PointType, Devices::Host > nabla_mesh( cellsCount );
    nabla_mesh = 0;
    for( int i = 0; i < cellsCount; i++ )
    {
        nabla_mesh[ i ][ 0 ] += dk1_nabla_h[ i ][ 0 ] * ( nabla_h[ i ][ 0 ] - nabla[ i ][ 0 ] )
                             + dk1_nabla_h[ i ][ 1 ] * ( nabla_h[ i ][ 1 ] - nabla[ i ][ 1 ] );
        nabla_mesh[ i ][ 1 ] += dk2_nabla_h[ i ][ 0 ] * ( nabla_h[ i ][ 0 ] - nabla[ i ][ 0 ] )
                             + dk2_nabla_h[ i ][ 1 ] * ( nabla_h[ i ][ 1 ] - nabla[ i ][ 1 ] );
    }

    std::cout << "mesh gradient" << std::endl; 
    for( int i = 0; i < cellsCount; i++ )
    {
        std::cout << nabla_mesh[ i ];
    }
    std::cout << std::endl;

    return true;
}

int main( int argc, char* argv[] )
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
            return numScheme(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
