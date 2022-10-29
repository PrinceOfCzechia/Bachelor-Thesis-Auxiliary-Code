#include <TNL/Meshes/Mesh.h>
#include <TNL/Meshes/Geometry/getEntityMeasure.h>
#include <TNL/Meshes/Geometry/getEntityCircumradius.h>
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
    using MeshType = Mesh< MeshConfig, Devices::/*Host*/Sequential >; // sequential bezi jednovlaknove
    using RealType = typename MeshType::RealType;
    using GlobalIndexType = typename MeshType::GlobalIndexType;
    using LocalIndexType = typename MeshType::LocalIndexType;
    using VectorType = TNL::Containers::Vector< RealType, TNL::Devices::Host, GlobalIndexType >;
    using PointType = typename MeshTraits< MeshConfig >::PointType;

    const auto verticesCount = mesh.template getEntitiesCount< 0 >();
    const auto facesCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() - 1 >();
    const auto cellsCount = mesh.template getEntitiesCount< MeshType::getMeshDimension() >();

    VectorType cellMeasures( cellsCount );
    VectorType interiorMeasures ( cellsCount );
    VectorType boundaryMeasures ( cellsCount );
    RealType Smeasures; RealType Sinterior; RealType Sboundary;
    VectorType faceMeasures( facesCount );
    //TNL::Containers::Vector< PointType, TNL::Devices::Host > normals;

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

    /*for(LocalIndexType j = 0; j < facesCount; j++)
    {
        const auto face = mesh.template getEntity< MeshType::getMeshDimension() - 1 >( j );
        faceMeasures[ j ] = ( getEntityMeasure( mesh, face ) );
    }*/

    std::cout << "Vertices count:" << verticesCount << std::endl;
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
            return iterateMesh( mesh, "/home/kral/Documents/TNL_work/unicorn.vtk" );
        };
        result &= resolveAndLoadMesh< MyConfigTag, Devices::Host >( wrapper, fileName );
    }

   return 0;
}
