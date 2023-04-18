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

template< typename t >
Containers::StaticVector< 2, t > dxsv( t x00, t x01, t x10, t x11, int point )
{
    Containers::StaticVector< 2, t > result;
    result[ 0 ] = d_x_sigma( x00, x01, x10, x11, point, 0 );
    result[ 1 ] = d_x_sigma( x00, x01, x10, x11, point, 1 );
    return result;
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
    if( point == 0 && component == 0) return d00inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 0 && component == 1) return d01inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 1 && component == 0) return d10inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 1 && component == 1) return d11inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 2 && component == 0) return d20inv_m2( x00, x01, x10, x11, x20, x21 );
    else if( point == 2 && component == 1) return d21inv_m2( x00, x01, x10, x11, x20, x21 );
    else return 0;
}

template< typename MeshConfig >
bool numScheme( Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
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

            //std::cout << "\nCell " << i << ", local face index: " << j << ", global face index: " << faceIdx << "\n";

            PointType faceVector = mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+2) % 3 ) ) 
                                 - mesh.getPoint( cell.template getSubentityIndex< 0 > ( (j+1) % 3 ) );

            PointType outwardNormal = normalize< PointType >( { faceVector[ 1 ], -faceVector[ 0 ] } );
            // std::cout << "Outward normal vector n_sigma: " << outwardNormal;

            PointType x_sigma = getEntityCenter( mesh, sigma );
            // std::cout << ", face center: " << x_sigma << "\n";

            double f_sigma = f< PointType >( x_sigma );
            sum += getEntityMeasure( mesh, sigma ) * f_sigma * outwardNormal;
        }

        PointType grad_h = ( 1.0 / getEntityMeasure( mesh, cell ) ) * sum;
        PointType grad = angrad< PointType >( getEntityCenter( mesh, cell ) );

        for( int j = 0; j < 3; j++ )
        {
        // PointType cellCenter = getEntityCenter( mesh, cell);
        // in center or in each individual point?
        // so far i chose to compute the angrad in each point
        int globalPointIdx = cell.template getSubentityIndex< 0 >( j );
        nabla[ globalPointIdx ] += grad;
        nabla_h[ globalPointIdx ] += grad_h;
        }
    }

    // std::cout << "nabla: " << nabla << "\n";
    // std::cout << "nabla_h: " << nabla_h << "\n";

    // calculating L(mesh)
    double L = 0.0;
    for( int i = 0; i < verticesCount; i++ )
    {
        L += l2Norm( nabla_h[i] - nabla[i] ) * l2Norm( nabla_h[i] - nabla[i] );
    }


    // derivatives of nabla_h f(x^i), see LaTeX
    // derivatives by the first component
    Containers::Vector< PointType, Devices::Host > dk1_nabla_h ( verticesCount );
    dk1_nabla_h = 0;
    //derivatives by the second component
    Containers::Vector< PointType, Devices::Host > dk2_nabla_h ( verticesCount );
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
        int globalPointIdx0 = cell.template getSubentityIndex< 0 >( 0 );

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
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 0, 0 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], 
                                      k, 0 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      ( 
                        angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                        dxsv< double >( cp[l][0],cp[l][1],cp[m][0],cp[m][1],k )
                      ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n0< double >( cp[ l ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 0 );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIdx0 ][ 0 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << "\n";

        // by the second component
        
        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 0, 1 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) )[ 0 ] *
                      d_x_sigma< double >(cp[l][0],cp[l][1],cp[m][0],cp[m][1],k,1) * //0 from d x_sigma
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n0< double >( cp[ l ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 1 );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIdx0 ][ 1 ] += d_nabla_h;
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

        int globalPointIdx1 = cell.template getSubentityIndex< 0 >( 1 );

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
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 1, 0 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 0 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                        angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                        dxsv< double >( cp[l][0],cp[l][1],cp[m][0],cp[m][1],k )
                      ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n0< double >( cp[ l ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 0 );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIdx1 ][ 0 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << "\n";

        // by the second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 1, 1 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                        cp[ m ][ 0 ],
                                        cp[ l ][ 1 ],
                                        cp[ m ][ 1 ],  k, 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[l][0], cp[l][1], cp[m][0], cp[m][1], k )
                      )* //0 from d x_sigma
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n0< double >( cp[ l ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 1 );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIdx1 ][ 1 ] += d_nabla_h;
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

        int globalPointIdx2 = cell.template getSubentityIndex< 0 >( 2 );

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
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 2, 0 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 0 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[l][0], cp[l][1], cp[m][0], cp[m][1], k )
                      ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n0< double >( cp[ l ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 0 );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIdx2 ][ 0 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << "\n";

        // by the second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 0 ],
                                    cp[ k ][ 1 ],
                                    cp[ l ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 2, 1 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                        cp[ m ][ 0 ],
                                        cp[ l ][ 1 ],
                                        cp[ m ][ 1 ],  k, 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[l][0], cp[l][1], cp[m][0], cp[m][1], k )
                      ) *
                      n0< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n0< double >( cp[ l ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 1 );
        }
        d_nabla_h += m2 * auxSum;

        dk1_nabla_h[ globalPointIdx2 ][ 1 ] += d_nabla_h;
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

        int globalPointIdx0 = cell.template getSubentityIndex< 0 >( 0 );

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
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 0, 0 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 0 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 0 ],
                                    cp[ m ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[l][0], cp[l][1], cp[m][0], cp[m][1], k )
                      )* // 0 from d x_sigma
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n1< double >( cp[ k ][ 0 ],
                                      cp[ k ][ 1 ],
                                      cp[ l ][ 0 ],
                                      cp[ l ][ 1 ], k, 0 );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIdx0 ][ 0 ] += d_nabla_h;

        // second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 0, 1 ) * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                        cp[ m ][ 0 ],
                                        cp[ l ][ 1 ],
                                        cp[ m ][ 1 ],  k, 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[ 0 ][ 0 ], cp[ 0 ][ 1 ],cp[ 1 ][ 0 ], cp[ 1 ][ 1 ], k )
                      ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 0; i < 3; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n1< double >( cp[ k ][ 0 ], // TODO
                                      cp[ k ][ 1 ],
                                      cp[ l ][ 0 ],
                                      cp[ l ][ 1 ], k, 1 );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIdx0 ][ 1 ] += d_nabla_h;
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

        int globalPointIdx1 = cell.template getSubentityIndex< 0 >( 1 );

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
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 1, 0 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 0 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ k ][ 0 ],
                                    cp[ l ][ 0 ],
                                    cp[ k ][ 1 ],
                                    cp[ l ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ k ], cp[ l ]) ),
                       dxsv< double >( cp[ 0 ][ 0 ], cp[ 0 ][ 1 ], cp[ 1 ][ 0 ], cp[ 1 ][ 1 ], k )
                      ) * // 0 from d x_sigma
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n1< double >( cp[ k ][ 0 ], // TODO
                                      cp[ k ][ 1 ],
                                      cp[ l ][ 0 ],
                                      cp[ l ][ 1 ], k, 0 );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIdx1 ][ 0 ] += d_nabla_h;

        // second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 1, 1 ) * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[ l ][ 0 ], cp[ l ][ 1 ], cp[ m ][ 0 ], cp[ m ][ 1 ], k )
                      ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 1; i < 4; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n1< double >( cp[ l ][ 0 ], // TODO
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 1 );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIdx1 ][ 1 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << "\n";
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

        int globalPointIdx2 = cell.template getSubentityIndex< 0 >( 2 );

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
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 2, 0 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 0 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[ l ][ 0 ], cp[ l ][ 1 ], cp[ m ][ 0 ], cp[ m ][ 1 ], k )
                      ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n1< double >( cp[ l ][ 0 ], // TODO
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 0 );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIdx2 ][ 0 ] += d_nabla_h;

        // second component

        d_nabla_h = 0;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ]  ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += d_inv_m2< double >( x00, x01, x10, x11, x20, x21, 2, 1 ) * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += d_m1< double >( cp[ l ][ 0 ],
                                      cp[ m ][ 0 ],
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 1 ], k, 1 ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      (
                       angrad< PointType >( x_sigma< PointType >(cp[ l ], cp[ m ]) ),
                       dxsv< double >( cp[ l ][ 0 ], cp[ l ][ 1 ], cp[ m ][ 0 ], cp[ m ][ 1 ], k )
                      ) *
                      n1< double >( cp[ l ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 0 ],
                                    cp[ m ][ 1 ] );
        }
        d_nabla_h += m2 * auxSum;

        auxSum = 0;
        for( int i = 2; i < 5; i++ )
        {
            int k = i % 3;
            int l = ( i + 1 ) % 3;
            int m = ( i + 2 ) % 3;
            auxSum += m1< double >( cp[ l ][ 0 ],
                                    cp[ m ][ 0 ],
                                    cp[ l ][ 1 ],
                                    cp[ m ][ 1 ] ) *
                      f< PointType >( x_sigma< PointType >( cp[ l ], cp[ m ] ) ) *
                      d_n1< double >( cp[ l ][ 0 ], // TODO
                                      cp[ l ][ 1 ],
                                      cp[ m ][ 0 ],
                                      cp[ m ][ 1 ], k, 1 );
        }
        d_nabla_h += m2 * auxSum;

        dk2_nabla_h[ globalPointIdx2 ][ 1 ] += d_nabla_h;

        // std::cout << "Cell " << j << ": " << d_nabla_h << "\n";
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

    // save vector of N PointTypes as an array of 3N doubles for writing
    // add zeros as a third component because of .vtk format
    Containers::Array< double > nabla_arr( 3 * verticesCount );
    nabla_arr = 0;
    for( int i = 0; i < 3 * verticesCount; i += 3 )
    {
        nabla_arr[ i ] = nabla_mesh[ i / 3 ][ 0 ];
        nabla_arr[ i + 1 ] = nabla_mesh[ i / 3 ][ 1 ];
        nabla_arr[ i + 2 ] = 0; // redundant
    }

    // writing the computed gradient into a new mesh
    using VTKWriter = Meshes::Writers::VTKWriter< MeshType >;
    std::ofstream out = std::ofstream( "new.vtk" );
    VTKWriter writer = VTKWriter( out );
    writer.template writeEntities< MeshType::getMeshDimension() >( mesh );
    writer.writePointData( nabla_arr, "meshGrads", 3 );



    // norm of the longest interior gradient (0 is a lower bound)
    double maxNorm = 0;
    auto getLongestGrad = [ &nabla_mesh, &maxNorm ] ( GlobalIndexType i )  // mutable, uncomment for possible future use
    {
        double aux = l2Norm( nabla_mesh[ i ] );
        if( aux > maxNorm ) maxNorm = aux;
    };
    mesh.template forAll< 0 >( getLongestGrad );

    // getting minimal face measure
    double minFace = std::numeric_limits< double >::max();
    auto getShortestFace = [ &mesh, &minFace ] ( GlobalIndexType i )  // mutable, uncomment for possible future use
    {
        auto face = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( i );
        double aux = getEntityMeasure( mesh, face );
        if( aux < minFace ) minFace = aux;
    };
    mesh.template forAll< MeshType::getMeshDimension() - 1 >( getShortestFace );

    // update interior vertices by -nabla_mesh
    auto updatePoints = [ &mesh, &nabla_mesh, &maxNorm, &minFace ] ( GlobalIndexType i ) mutable
    {
        mesh.getPoints()[ i ] -= 1e-9 * minFace * maxNorm * nabla_mesh[ i ];
    };
    mesh.template forAll< 0 >( updatePoints );

    // writing the updated mesh into a new file
    std::ofstream output = std::ofstream( "iteration.vtk" );
    VTKWriter writer2 = VTKWriter( output, TNL::Meshes::VTK::FileFormat::ascii );
    writer2.template writeEntities< MeshType::getMeshDimension() >( mesh );

    Containers::Vector< PointType > new_nabla_h( verticesCount );
    Containers::Vector< PointType > new_nabla( verticesCount );

    // give the points vector value of respective vertices
    Containers::Vector< PointType > points( verticesCount );
    auto fillPoints = [ &mesh, &points ] ( GlobalIndexType i )  // mutable, uncomment for possible future use
    {
        points[ i ] = mesh.getPoints( )[ i ];
    };
    mesh.template forAll< 0 >( fillPoints );

    // calculate nabla_h and nabla after the first iteration
    auto next_iteration = [ &mesh, &points, &new_nabla_h, &new_nabla ] ( GlobalIndexType i )  mutable
    {
        auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
        PointType sum = { 0, 0 };
        for( int j = 0; j < 3; j++ )
        {
            PointType faceVector = points[ cell.template getSubentityIndex< 0 > ( (j+2) % 3 ) ] 
                                 - points[ cell.template getSubentityIndex< 0 > ( (j+1) % 3 ) ];

            PointType outwardNormal = normalize< PointType >( { faceVector[ 1 ], -faceVector[ 0 ] } );

            PointType x_sigma = 0.5 *
                                ( points[ cell.template getSubentityIndex< 0 > ( (j+1) % 3 ) ] 
                                + points[ cell.template getSubentityIndex< 0 > ( (j+2) % 3 ) ] );

            double f_sigma = f< PointType >( x_sigma );
            sum += l2Norm( faceVector ) * f_sigma * outwardNormal;
        }
        double x00 = points[ cell.template getSubentityIndex< 0 >( 0 ) ][ 0 ];
        double x01 = points[ cell.template getSubentityIndex< 0 >( 0 ) ][ 1 ];
        double x10 = points[ cell.template getSubentityIndex< 0 >( 1 ) ][ 0 ];
        double x11 = points[ cell.template getSubentityIndex< 0 >( 1 ) ][ 1 ];
        double x20 = points[ cell.template getSubentityIndex< 0 >( 2 ) ][ 0 ];
        double x21 = points[ cell.template getSubentityIndex< 0 >( 2 ) ][ 1 ]; 
        PointType grad_h = inv_m2< double >( x00, x01, x10, x11, x20, x21 ) * sum;
        PointType center = 0;
        center[ 0 ] = (x00 + x10 + x20) / 3;
        center[ 1 ] = (x01 + x11 + x21) / 3;
        PointType grad = angrad< PointType >( center );
        
        for( int j = 0; j < 3; j++ )
        {
        int globalPointIdx = cell.template getSubentityIndex< 0 >( j );
        new_nabla[ globalPointIdx ] += grad;
        new_nabla_h[ globalPointIdx ] += grad_h;
        }
    };
    mesh.template forAll< MeshType::getMeshDimension() >( next_iteration );

    // std::cout << "new nabla: " << new_nabla << "\n";
    // std::cout << "new nabla_h: " << new_nabla_h << "\n";

    // calculating L(after optimization step)
    double new_L = 0.0;
    for( int i = 0; i < verticesCount; i++ )
    {
        new_L += l2Norm( new_nabla_h[ i ] - new_nabla[ i ] ) * l2Norm( new_nabla_h[ i ] - new_nabla[ i ] );
    }

    std::cout << nabla_mesh << "\n";

    std::cout << "L(initial mesh) = " << L << "\n";
    std::cout << "L(after update) = " << new_L << "\n";
    std::cout << "longest interior gradient: " << maxNorm << "\n";
    std::cout << "shortest face in the mesh: " << minFace << "\n";

    std::cout << "\nimprovement of L: " << L - new_L << "\n";



    std::cout << "OK" << "\n";
    return true; 
}

// TODO:
// more iterations of GD
// Tikhonov regularization
// bet border vertices, forAll cycle
// use getPoints()

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
            return numScheme(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
