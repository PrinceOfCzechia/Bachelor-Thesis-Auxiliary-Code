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

template< typename t >
t df( t x, t y, bool component )
{
    return (component ? 2*x*x*y : 2*x*y*y);
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

// face center
template< typename V >
V x_sigma( V v1, V v2)
{
    V sum = v1 + v2;
    return 0.5 * sum;
}

template< typename t >
Containers::Vector< t > x_sigma( t x00, t x01, t x10, t x11 )
{
    Containers::Vector< t > result( 2 );
    result[ 0 ] = 0.5 * ( x00 + x10 );
    result[ 1 ] = 0.5 * ( x01 + x11 );
    return result;
}

template< typename t >
t x_sigma_0( t x00, t x10 ) { return 0.5 * ( x00 + x10 ); }

template< typename t >
t x_sigma_1(t x01, t x11 ) { return 0.5 * ( x01 + x11 ); }

// TODO maybe remove arguments, but this is more instructive
template< typename t >
t d0_x_sigma_0( t x00, t x10 ) { return 0.5; }

template< typename t >
t d0_x_sigma_1( t x00, t x10 ) { return 0; }

template< typename t >
t d1_x_sigma_0( t x00, t x10 ) { return 0; }

template< typename t >
t d1_x_sigma_1( t x00, t x10 ) { return 0.5; }

template< typename t >
t d_x_sigma( t x00, t x01, t x10, t x11, int point, bool component )
{
    if( point == 0 && component == 0 ) return d0_x_sigma_0( x00, x10 );
    else if( point == 1 && component == 0  ) return d1_x_sigma_0( x00, x10 );
    else if( point == 0 && component == 1 ) return d0_x_sigma_1( x10, x11 );
    else if( point == 1 && component == 1 ) return d1_x_sigma_1( x10, x11 );
    else return 0;
}

// face measure
template< typename t >
t m1( t x00, t x01, t x10, t x11 )
{ return pow( ( (x10 - x00) * (x10 - x00) + (x11 - x01) * (x11 - x01) ), 0.5 ); }

// derivative of the face measure by the first component
template< typename t >
t d_m1_0( t x00, t x01, t x10, t x11, bool direction )
{ return 1 / m1( x00, x01, x10, x11 ) * ( x10 - x00 ) * ( direction ? 1 : -1 ); }

// derivative of the face measure by the second component
template< typename t >
t d_m1_1( t x00, t x01, t x10, t x11, bool direction )
{ return 1 / m1( x00, x01, x10, x11 ) * ( x11 - x01 ) * ( direction ? 1 : -1 ); }
// *** direction 1 ~ derivative by x10 or x11, direction 0 ~ derivative by x00 or x01 ***

template< typename t >
t d_m1( t x00, t x01, t x10, t x11, int point, bool component  )
{
    if( point == 0 && component == 0 ) return d_m1_0( x00, x01, x10, x11, component );
    else if( point == 1 && component == 0 ) return d_m1_0( x00, x01, x10, x11, !component );
    else if( point == 0 && component == 1 ) return d_m1_1( x00, x01, x10, x11, component );
    else if( point == 1 && component == 1 ) return d_m1_1( x00, x01, x10, x11, !component );
    else return 0;
}

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

template< typename t >
t d_n0( t x00, t x01, t x10, t x11, int point, bool component)
{
    if( point == 0 && component == 0 ) return d0n00( x00, x01, x10, x11 );
    else if( point == 1 && component == 0 ) return d1n00( x00, x01, x10, x11 );
    else if( point == 0 && component == 1 ) return d0n01( x00, x01, x10, x11 );
    else if( point == 1 && component == 1 ) return d1n01( x00, x01, x10, x11 );
    else return 0;
}

template< typename t >
t d_n1( t x00, t x01, t x10, t x11, int point, bool component)
{
    if( point == 0 && component == 0 ) return d0n10( x00, x01, x10, x11 );
    else if( point == 1 && component == 0 ) return d1n10( x00, x01, x10, x11 );
    else if( point == 0 && component == 1 ) return d0n11( x00, x01, x10, x11 );
    else if( point == 1 && component == 1 ) return d1n11( x00, x01, x10, x11 );
    else return 0;
}

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

template< typename t >
t d_inv_m2( t x00, t x01, t x10, t x11, t x20, t x21, int point, bool component )
{
    if( point == 0 & component == 0) return d00inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 0 & component == 1) return d01inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 1 & component == 0) return d10inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 1 & component == 1) return d11inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 2 & component == 0) return d20inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 2 & component == 1) return d21inv_m2( x00, x01, x10, x11, x20, x21 );
    else return 0;
}

template< typename PointType >
auto diff_num( const Containers::Vector< PointType >& y, int point, bool component ) -> Containers::Vector< typename PointType::RealType >
{
    using t = typename PointType::RealType;
    Containers::Vector< Containers::Vector< t > > x ( 3 );
    for( int i = 0; i < 3; i++ )
    {
        x[ i ][ 0 ] = y[ i ][ 0 ];
        x[ i ][ 1 ] = y[ i ][ 1 ];
    }
    Containers::Vector< t > sum ( 2 );
    sum = 0;
    Containers::Vector< t > auxSum ( 2 );
    auxSum = 0;
    for( int i = 0; i < 3; i++ )
    {
        int k = i%3;
        int l = (i+1)%3;
        int m = (i+2)%3;
        Containers::Vector< t > n ( 2 );
        n[ 0 ] = n0< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] );
        n[ 1 ] = n1< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] );
        auxSum += m1< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] ) *
                  f< Containers::Vector< t > >(
                    x_sigma< t >( x[ l ][ 0 ],x[ l ][ 1 ], x[ m ][ 0 ],x[ m ][ 1 ] ) ) *
                  n;
    }
    t d = d_inv_m2< t >( x[ 0 ][ 0 ], x[ 0 ][ 1 ], x[ 1 ][ 0 ], x[ 1 ][ 1 ], x[ 2 ][ 0 ], x[ 2 ][ 1 ], point, component );
    auxSum *= d;
    sum += auxSum;
    auxSum = 0;

    for( int i = 0; i < 3; i++ )
    {
        int k = i%3;
        int l = (i+1)%3;
        int m = (i+2)%3;
        Containers::Vector< t > d_n ( 2 );
        d_n[ 0 ] = d_n0 < t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ], point, component ),
        d_n[ 1 ] = d_n1< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ], point, component );
        auxSum += m1< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] ) *
                  f< Containers::Vector< t > >(
                   x_sigma< Containers::Vector< t > >( {x[ l ][ 0 ],x[ l ][ 1 ]},
                                                            {x[ m ][ 0 ],x[ m ][ 1 ]} )
                                                   ) *
                  d_n;
    }
    d = inv_m2< t >( x[ 0 ][ 0 ], x[ 0 ][ 1 ], x[ 1 ][ 0 ], x[ 1 ][ 1 ], x[ 2 ][ 0 ], x[ 2 ][ 1 ] );
    auxSum *= d;
    sum += auxSum;
    auxSum = 0;

    for( int i = 0; i < 3; i++ )
    {
        int k = i%3;
        int l = (i+1)%3;
        int m = (i+2)%3;
        Containers::Vector< t > n ( 2 );
        n[ 0 ] = n0< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] ),
        n[ 1 ] = n1< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] );
        auxSum += d_m1< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ], point, component ) *
                  f< Containers::Vector< t > >(
                   x_sigma< Containers::Vector< t > >( {x[ l ][ 0 ],x[ l ][ 1 ]},
                                                            {x[ m ][ 0 ],x[ m ][ 1 ]} )
                                                   ) *
                  n;
    }
    d = inv_m2< t >( x[ 0 ][ 0 ], x[ 0 ][ 1 ], x[ 1 ][ 0 ], x[ 1 ][ 1 ], x[ 2 ][ 0 ], x[ 2 ][ 1 ] );
    auxSum *= d;
    sum += auxSum;
    auxSum = 0;

    for( int i = 0; i < 3; i++ )
    {
        int k = i%3;
        int l = (i+1)%3;
        int m = (i+2)%3;
        Containers::Vector< t > dx_sigma ( 2 );
        dx_sigma[ 0 ] = d_x_sigma< t >( x[ l ][ 0 ],
                                 x[ l ][ 1 ],
                                 x[ m ][ 0 ],
                                 x[ m ][ 1 ], point, 0 );
        dx_sigma[ 1 ] = d_x_sigma< t >( x[ l ][ 0 ],
                                 x[ l ][ 1 ],
                                 x[ m ][ 0 ],
                                 x[ m ][ 1 ], point, 1 );
        Containers::Vector< t > n ( 2 );
        n[ 0 ] = n0( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] );
        n[ 1 ] = n1( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] );
        auxSum += m1< t >( x[ l ][ 0 ], x[ l ][ 1 ], x[ m ][ 0 ], x[ m ][ 1 ] ) *
                  df< t >( x_sigma< t >( x[l][0], x[l][1], x[m][0], x[m][1] )[0],
                                x_sigma< t >( x[l][0], x[l][1], x[m][0], x[m][1] )[1],
                                component ) *
                  dx_sigma *
                  n;
    }
    d = inv_m2< t >( x[ 0 ][ 0 ], x[ 0 ][ 1 ], x[ 1 ][ 0 ], x[ 1 ][ 1 ], x[ 2 ][ 0 ], x[ 2 ][ 1 ] );
    auxSum *= d;
    sum += auxSum;
    auxSum = 0;

    return sum;
}

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

    std::cout << "Cells: " << cellsCount << ", faces: " << facesCount << ", vertices: " << verticesCount << std::endl;

    Containers::Vector< PointType, Devices::Host > nabla_h ( verticesCount );
    Containers::Vector< PointType, Devices::Host > nabla ( verticesCount );
    nabla_h = 0;
    nabla = 0;

    for( int i = 0; i < cellsCount; i++ )
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType sum = { 0, 0 };
        for( int j = 0; j < 3; j++ )
        {
            const auto faceIdx = cell.template getSubentityIndex< MeshType::getMeshDimension() - 1 >( j );
            const auto sigma = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( faceIdx );

            //std::cout << "\nCell " << i << ", local face index: " << j << ", global face index: " << faceIdx << std::endl;

            PointType faceVector = mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2) % 3 ) ) 
                                 - mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1) % 3 ) );

            PointType outwardNormal = normalize< PointType >( { faceVector[ 1 ], -faceVector[ 0 ] } );
            // std::cout << "Outward normal vector n_sigma: " << outwardNormal;

            PointType x_sigma = getEntityCenter( mesh, sigma );
            // std::cout << ", face center: " << x_sigma << "\n";

            double f_sigma = f< PointType >( x_sigma );
            sum += getEntityMeasure( mesh, sigma ) * f_sigma * outwardNormal;
        }

        PointType grad = ( 1.0 / getEntityMeasure( mesh, cell ) ) * sum;

        for( int j = 0; j < 3; j++ )
        {
        // PointType cellCenter = getEntityCenter( mesh, cell);
        // in center or in each individual point?
        // so far i chose to compute the angrad in each point
        int globalPointIdx = cell.template getSubentityIndex< 0 >( j );
        nabla[ globalPointIdx ] += angrad< PointType >( mesh.template getPoint( globalPointIdx ) );
        nabla_h[ globalPointIdx ] += grad;
        }
    }

    // derivatives of nabla_h f(x^i), see LaTeX
    // derivatives by the first component
    Containers::Vector< PointType, Devices::Host > dk1_nabla_h ( verticesCount );
    dk1_nabla_h = 0;
    //derivatives by the second component
    Containers::Vector< PointType, Devices::Host > dk2_nabla_h ( verticesCount );
    dk2_nabla_h = 0;

    // point 0
    for( int j = 0; j < 1/*TODO cellsCount*/; j++ )
    {
        Containers::Vector< PointType > cp ( 3 ); // cellPoints
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( j );
        for(int i = 0; i < 1/*TODO 3*/; i++)
        {
            int globalPointIdx = cell.template getSubentityIndex< 0 >( i );
            PointType point = mesh.getPoint( globalPointIdx );
            cp[ i ] = point;
        

            Containers::Vector< double > diff0( 2 );
            diff0 = diff_num( cp, i%3, 0 );
            Containers::Vector< double > diff1( 2 );
            diff1 = diff_num( cp, i%3, 1 );

            std::cout << diff0 << std::endl;
            std::cout << diff1 << std::endl;

            dk1_nabla_h[ globalPointIdx ] += diff0;
            dk2_nabla_h[ globalPointIdx ] += diff1;
        }
    }

    Containers::Vector< PointType, Devices::Host > nabla_mesh( verticesCount );
    nabla_mesh = 0;
    for( int i = 0; i < verticesCount; i++ )
    {
        nabla_mesh[ i ][ 0 ] += dk1_nabla_h[ i ][ 0 ] * ( nabla_h[ i ][ 0 ] - nabla[ i ][ 0 ] )
                             + dk1_nabla_h[ i ][ 1 ] * ( nabla_h[ i ][ 1 ] - nabla[ i ][ 1 ] );
        nabla_mesh[ i ][ 1 ] += dk2_nabla_h[ i ][ 0 ] * ( nabla_h[ i ][ 0 ] - nabla[ i ][ 0 ] )
                             + dk2_nabla_h[ i ][ 1 ] * ( nabla_h[ i ][ 1 ] - nabla[ i ][ 1 ] );
    }

    // save vector of N PointTypes as an array of 2N doubles for writing
    // add zeros as a third component because of .vtk format
    Containers::Array< double, Devices::Host > nabla_arr( 3 * verticesCount );
    nabla_arr = 0;
    for( int i = 0; i < 3 * verticesCount; i += 3 )
    {
        nabla_arr[ i ] = nabla_mesh[ i / 3 ][ 0 ];
        nabla_arr[ i + 1 ] = nabla_mesh[ i / 3 ][ 1 ];
        nabla_arr[ i + 2 ] = 0; // redundant
    }
 
    std::cout << "mesh gradient" << std::endl; 
    for( int i = 0; i < verticesCount; i++ )
    {
        std::cout << nabla_mesh[ i ];
    }
    std::cout << std::endl;

    // writing the result into a new mesh
    using VTKWriter = Meshes::Writers::VTKWriter< MeshType >;
    std::ofstream out = std::ofstream( "new.vtk" );
    VTKWriter writer = VTKWriter( out );
    writer.template writeEntities< MeshType::getMeshDimension() >( mesh );
    writer.writePointData( nabla_arr, "meshGrads", 3 );

    std::cout << "OK" << std::endl;

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
