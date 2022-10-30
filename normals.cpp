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


template< typename MeshConfig >
bool iterateMesh( Mesh< MeshConfig, Devices::Host >& mesh, const std::string& fileName )
{
    using MeshType = Mesh< MeshConfig, Devices::/*Host*/Host >; // Host bezi jednovlaknove
    using RealType = typename MeshType::RealType;
    using GlobalIndexType = typename MeshType::GlobalIndexType;
    using LocalIndexType = typename MeshType::LocalIndexType;
    using VectorType = TNL::Containers::Vector< RealType, TNL::Devices::Host, GlobalIndexType >;
    using PointType = typename MeshTraits< MeshConfig >::PointType;

    const auto verticesCount = mesh.template getEntitiesCount< 0 >();
    const auto facesCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() - 1 >();
    const auto cellsCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() >();

    using StaticArrayType = TNL::Containers::StaticArray< 10, VectorType >;

    VectorType cellMeasures( cellsCount );
    VectorType interiorMeasures ( cellsCount );
    VectorType boundaryMeasures ( cellsCount );
    RealType Smeasures; RealType Sinterior; RealType Sboundary;
    VectorType faceMeasures( facesCount );

    mesh.template forAll< MeshType::getMeshDimension() >( [&] (GlobalIndexType i)
        {
            const auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
            cellMeasures[ i ] = ( getEntityMeasure( mesh, cell ) );
            Smeasures += cellMeasures[ i ];
        }
    );

    // pro paralelizaci pouzit normalni for cyklus !!!!

    mesh.template forAll< MeshType::getMeshDimension() - 1 >( [&] (LocalIndexType j)
        {
            const auto face = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( j );
            faceMeasures[ j ] = ( getEntityMeasure( mesh, face ) );
        }
    );

    mesh.template forInterior< MeshType::getMeshDimension() >( [&] (GlobalIndexType i)
        {
            const auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
            interiorMeasures[ i ] = ( getEntityMeasure( mesh, cell ) );
            Sinterior += interiorMeasures[ i ];
        }
    );

    mesh.template forBoundary< MeshType::getMeshDimension() >( [&] (GlobalIndexType i)
        {
            const auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
            boundaryMeasures[ i ] = ( getEntityMeasure( mesh, cell ) );
            Sboundary += boundaryMeasures[ i ];
        }
    );

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

    std::cout << "\nVertices count:" << verticesCount << std::endl;
    std::cout << "Faces count:" << facesCount << std::endl;
    std::cout << "Cells count:" << cellsCount << std::endl;

    std::cout << "\nCells and their respective measures:";
    for(int i = 0; i < cellsCount; i++)
    {   
        if (i%4==0) std::cout << std::endl;
        std::cout << "cell " << i << ": " << cellMeasures[i] << "\t";
    }
    std::cout << "\nFaces and their respective measures:";
    for(int i = 0; i < facesCount; i++)
    {
        if (i%4==0) std::cout << std::endl;
        std::cout << "face " << i << ": " << faceMeasures[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "Measure of all the elements is " << Smeasures << ", where measure of boundary elements is "
              << Sboundary << " and measure of interior elements is " << Sinterior
              << ", difference between measure of all cells and sum of boundary and interior measures being "
              << std::abs(Smeasures - Sboundary - Sinterior) << std::endl;
    return true;
}



template< typename V >
V normalize(V vector)
{
    V u = (1.0 / l2Norm( vector )) * vector;
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

    //using Vector = StaticVector< 3, double >;
    /*using StaticVectorType = TNL::Containers::StaticVector< 8, PointType >;
    StaticVectorType points = 0;
    double measures = 0;
    mesh.template forAll< 0 >( [&] (GlobalIndexType i)
        {
            const auto vertex = mesh.template getEntity< 0 >( i );
            points[i] = ( mesh.getPoint(i) );
        }
    );
    mesh.template forAll< MeshType::getMeshDimension() >( [&] (GlobalIndexType i)
        {
            const auto cell = mesh.template getEntity< MeshType::getMeshDimension() >( i );
            measures += ( getEntityMeasure( mesh, cell ) );
        }
    );
    std::cout << "* " << measures << " *\n";
    std::cout << "Mesh points: " << points << std::endl;*/

    /*PointType vertex = mesh.template getPoint(0);
    std::cout << vertex << std::endl;*/
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

    /*auto kernel_measures = [] __host__ __device__
    ( GlobalIndex cell_idx, const Mesh& mesh, Real* output_array )
    {
    const auto& cell = mesh.template getEntity<2>( cell_idx );
    const auto& v0 = mesh.getPoint( cell.template getSubentityIndex<0>( 0 ) );
    const auto& v1 = mesh.getPoint( cell.template getSubentityIndex<0>( 1 ) );
    const auto& v2 = mesh.getPoint( cell.template getSubentityIndex<0>( 2 ) );
    output_array[ cell_idx ] = getTriangleArea( v2 - v0, v1 - v0 );
    };*/

    std::cout << "\n" << verticesVector << std::endl;
    std::cout << "\nOutward normal vectors (row ~ cellIdx, col ~ local faceIdx): " << std::endl;
    for(int i = 0; i < cellsCount; i++)
    {
        std::cout << normals[i] << std::endl;
    }

    return true;


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
}

/*Containers::Vector< double, Devices::Host > normalize(TNL::Containers::Vector< double > v)
{
    return (1/TNL::l2norm(v)) * v;
}*/

/*template< typename Device >
void normalizeVector(Containers::Vector< double, Device > v)
{
	v *= ( 1.0 / l2Norm( v ) );
}*/

int
main( int argc, char* argv[] )
{
    /*Containers::Vector< double, Devices::Host > v = { 4, 3 };
	std::cout << v << "of 2-norm " << l2Norm( v ) << std::endl;
    std::cout << v.getSize() << std::endl;
    Containers::Vector< double, Devices::Host > u;
    ( 1.0 / l2Norm( v ) ) * v normalizeVector( v );
    std::cout << v << "of 2-norm " << l2Norm( v ) <<  std::endl;*/
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
            return getNormals(mesh, "");
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
