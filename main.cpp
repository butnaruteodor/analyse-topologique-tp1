#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/Color.h"
#include <DGtal/geometry/curves/GreedySegmentation.h>
#include <iterator>
#include <numeric>

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef PointVector<2,int> Point;
typedef std::vector<Point> Range;
typedef Range::const_iterator ConstIterator;
typedef ArithmeticalDSSComputer<ConstIterator, int, 4> SegmentComputer;
typedef GreedySegmentation<SegmentComputer> Segmentation;

typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet; // Digital set type
typedef Object<DT4_8, DigitalSet> ObjectType; // Digital object with (4,8)-adjacency pair

template<class T>
Curve getBoundary(T & object)
{
    Curve boundaryCurve;
    // Safety check for empty objects
    if (object.pointSet().size() == 0) return Curve();
    //Khalimsky space (cubical complex)
    try {
        KSpace kSpace;
        // we need to add a margin to prevent situations such that an object touch the border of the domain
        kSpace.init( object.domain().lowerBound() - Point(1,1),
                    object.domain().upperBound() + Point(1,1), true );

        // 1) Call Surfaces::findABel() to find a cell which belongs to the border
        const DigitalSet& set = object.pointSet();
        SCell s = Surfaces<KSpace>::findABel(kSpace, set);

        // // 2) Call Surfece::track2DBoundaryPoints to extract the boundary of the object (boundaryPoints)
        std::vector<SCell> boundarySurfels;
        SurfelAdjacency<2> sAdj(true);
        
        Surfaces<KSpace>::track2DBoundary(boundarySurfels, kSpace, sAdj, set, s);

        // 3) Create "boundaryCurve" from "boundarySurfels"
        boundaryCurve.initFromSCellsVector(boundarySurfels);
    }catch(...)
    {
        return Curve();
    }

    return boundaryCurve;
}
double calculatePolygonArea(const std::vector<Point>& vertices) {
    double area = 0.0;
    int n = vertices.size();
    // Apply formula: sum (x_i * y_i+1 - x_i+1 * y_i)
    for (int i = 0; i < n; ++i) {
        Point p1 = vertices[i];
        Point p2 = vertices[(i + 1) % n];
        area += (p1[0] * p2[1] - p2[0] * p1[1]);
    }
    return std::abs(area) / 2.0;
}

void getAreaDistribution(std::vector<ObjectType>& objects) 
{
    std::cout << "Analyzing " << objects.size() << " objects..." << std::endl;
    std::vector<double> areas;
    // Z2i::Curve c = getBoundary(objects[0]);
    for(auto& obj : objects)
    {
        Z2i::Curve c = getBoundary(obj);
        Range r; 
        auto curvePoints = c.getPointsRange(); 
        std::copy(curvePoints.begin(), curvePoints.end(), std::back_inserter(r));
        SegmentComputer recognitionAlgorithm;
        Segmentation theSegmentation(r.begin(), r.end(), recognitionAlgorithm);
            
        Segmentation::SegmentComputerIterator beg = theSegmentation.begin();
        Segmentation::SegmentComputerIterator end = theSegmentation.end();

        std::vector<Point> polygonVertices;
        for ( ; beg != end; ++beg) {
            Point vertex = *beg->begin(); 
            polygonVertices.push_back(vertex);
        }

        // Calculate Area
        double polyArea = calculatePolygonArea(polygonVertices);
        if (polyArea > 0.0)
        {
            areas.push_back(polyArea);
        }
    }
    auto mean = std::accumulate(areas.begin(), areas.end(), 0.0) / areas.size(); 
    double sq_sum = std::inner_product(areas.begin(), areas.end(), areas.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / areas.size() - mean * mean);
    cout << mean << " "<< " "<<stdev;
}

int main(int argc, char** argv)
{
    Board2D aBoard;
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings

    //typedef Object<DT8_4, DigitalSet> ObjectType; // Digital object with (8,4)-adjacency pair
    typedef ImageSelector<Domain, unsigned char >::Type Image; // type of image
    // 1) read an image
    Image image = PGMReader<Image>::importPGM ("../RiceGrains/Rice_japonais_seg_bin.pgm");

    // 2) Use SetFromImage::append() to convert the image to a "DigitalSet" of the proper size
    Z2i::DigitalSet set2d (image.domain());
    SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 1, 255);

    // 3) Create a digital object from the "DigitalSet"
    std::vector< ObjectType > objects; // All connected components are going to be stored in this vector
    std::back_insert_iterator< std::vector< ObjectType > > inserter( objects ); // Iterator used to populated "objects".
    ObjectType obj(Z2i::dt4_8, set2d);
    
    // 4) Use the method "writeComponents to obtain connected components with the proper adjacency pair
    obj.writeComponents(inserter);
    std::vector<ObjectType> completeObjects;
    Domain domain = image.domain();
    Point pMin = domain.lowerBound();
    Point pMax = domain.upperBound();

    for (const auto& obj : objects) {
        bool touchesBoundary = false;
        
        // Iterate through all points in the digital set of this component
        for (auto it = obj.pointSet().begin(); it != obj.pointSet().end(); ++it) {
            Point p = *it;
            
            // Check if the point is on the edge of the image frame
            if (p[0] == pMin[0] || p[0] == pMax[0] || 
                p[1] == pMin[1] || p[1] == pMax[1]) {
                touchesBoundary = true;
                break; // No need to check other points for this grain
            }
        }

        if (!touchesBoundary) {
            completeObjects.push_back(obj);
        }
    }

    // STEP 2

    std::cout << " number of components : " << completeObjects.size() << endl;

    // DISPLAY STEP 2

    // for (const auto& obj : completeObjects) {
    // aBoard << obj; 
    // }

    // STEP 3

    // Z2i::Curve c = getBoundary(completeObjects[0]);

    // DISPLAY STEP 3

    // aBoard << completeObjects[0];
    // aBoard << c;

    // complete the function "getBoundary" written above and add/modify your codes
    // STEP 4
    // Range r; 
    // auto curvePoints = c.getPointsRange(); 
    // std::copy(curvePoints.begin(), curvePoints.end(), std::back_inserter(r));
    // SegmentComputer recognitionAlgorithm;
    // Segmentation theSegmentation(r.begin(), r.end(), recognitionAlgorithm);
        
    // Segmentation::SegmentComputerIterator i = theSegmentation.begin();
    // Segmentation::SegmentComputerIterator end = theSegmentation.end();

    // DISPLAY STEP 4
    // aBoard << completeObjects[0];
    // aBoard << c;
    // for ( ; i != end; ++i) {
    //     // This holds the actual shape (ArithmeticalDSS) that supports "BoundingBox"
    //     SegmentComputer::Primitive currentShape = i->primitive();
        
    //     // Draw the Bounding Box using the Shape object
    //     aBoard << SetMode(currentShape.className(), "BoundingBox")
    //            << CustomStyle(currentShape.className() + "/BoundingBox", new CustomPen(Color::Red, Color::None, 2.0))
    //            << currentShape;
    // }

    // STEP 5
    // unsigned int digitalArea = completeObjects[0].pointSet().size();
    // std::cout << "Digital Area (Pixel Count): " << digitalArea << std::endl;

    // std::vector<Point> polygonVertices;
    // for ( ; i != end; ++i) {
    //     Point vertex = *i->begin(); 
    //     polygonVertices.push_back(vertex);
    // }

    // // Calculate Area
    // double polyArea = calculatePolygonArea(polygonVertices);
    // std::cout << "Polygon Area (Shoelace):    " << polyArea << std::endl;

    getAreaDistribution(completeObjects);
    // DISPLAY STEP 5

    // aBoard << completeObjects[0];
    // aBoard << c;
    // for ( ; i != end; ++i) {
    //     // This holds the actual shape (ArithmeticalDSS) that supports "BoundingBox"
    //     SegmentComputer::Primitive currentShape = i->primitive();
        
    //     // Draw the Bounding Box using the Shape object
    //     aBoard << currentShape;
        
    // }

    aBoard.saveCairo("./output.pdf", Board2D::CairoPDF);

    return 0;
}
