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

struct RiceFeatures {
    std::vector<double> polyAreas;      // Polygon Area (convergent)
    std::vector<double> polyPerimeters; // Polygon Perimeter (convergent)
    std::vector<double> digitalAreas;   // Pixel count
    std::vector<double> digitalPerimeters; // 1-cells count
    std::vector<double> circularities;  // 4*pi*A / P^2
};

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

double calculatePolygonPerimeter(const std::vector<Point>& vertices) {
    double perimeter = 0.0;
    int n = vertices.size();
    for (int i = 0; i < n; ++i) {
        Point p1 = vertices[i];
        Point p2 = vertices[(i + 1) % n];
        double dx = p1[0] - p2[0];
        double dy = p1[1] - p2[1];
        perimeter += std::sqrt(dx * dx + dy * dy);
    }
    return perimeter;
}

RiceFeatures analyzeRiceGrains(const std::vector<ObjectType>& completeObjects)
{
    // These are guaranteed to stay in sync
    RiceFeatures features;

    std::cout << "Starting analysis on " << completeObjects.size() << " objects..." << std::endl;

    for (const auto& obj : completeObjects)
    {
        Z2i::Curve c = getBoundary(obj);

        if (c.isValid() && c.size() > 2)
        {
            std::vector<Point> polygon;
            double currentArea = 0.0; 
            double digitalArea = (double)obj.pointSet().size();
            double currentPerim = 0.0;
            double digitalPerimeter = (double)c.size();
            
            Range r; 
            auto curvePoints = c.getPointsRange(); 
            std::copy(curvePoints.begin(), curvePoints.end(), std::back_inserter(r));

            SegmentComputer recognitionAlgorithm;
            Segmentation theSegmentation(r.begin(), r.end(), recognitionAlgorithm);

            std::vector<Point> polygonVertices;
            for (auto it = theSegmentation.begin(); it != theSegmentation.end(); ++it) {
                polygonVertices.push_back(*it->begin());
            }
            
            currentArea = calculatePolygonArea(polygonVertices);
            currentPerim = calculatePolygonPerimeter(polygonVertices);

            if (currentPerim > 0) {
                 features.polyAreas.push_back(currentArea);
                 features.polyPerimeters.push_back(currentPerim);
                 
                 features.digitalAreas.push_back(digitalArea);
                 features.digitalPerimeters.push_back(digitalPerimeter);
                 
                 // Calculate Circularity
                 double circ = (4.0 * M_PI * currentArea) / (currentPerim * currentPerim);
                 features.circularities.push_back(circ);
            }
        }
    }

    std::cout << "Successfully analyzed " << features.polyAreas.size() << " valid grains." << std::endl;
    return features;
}

void exportToCSV(const RiceFeatures& features, const std::string& filename)
{
    std::ofstream file(filename);
    
    // 1. Write the Header
    file << "ID,PolyArea,PolyPerimeter,DigitalArea,DigitalPerimeter,Circularity\n";

    // 2. Write the Data
    // (We assume all vectors in 'features' are the same size because of your synchronization logic)
    size_t count = features.polyAreas.size();
    
    for (size_t i = 0; i < count; ++i) {
        file << i << "," 
             << features.polyAreas[i] << "," 
             << features.polyPerimeters[i] << "," 
             << features.digitalAreas[i] << "," 
             << features.digitalPerimeters[i] << "," 
             << features.circularities[i] << "\n";
    }

    file.close();
    std::cout << "Data successfully exported to " << filename << std::endl;
}

void getIQR(std::vector<double> data, double& lowThreshold, double& highThreshold) {
    if (data.empty()) return;
    
    // Sort a copy to find percentiles
    std::sort(data.begin(), data.end());
    
    size_t n = data.size();
    double q1 = data[n / 4];
    double q3 = data[n * 3 / 4];
    double iqr = q3 - q1;
    
    lowThreshold = q1 - 1.5 * iqr;
    highThreshold = q3 + 1.5 * iqr;
}

RiceFeatures cleanOutliers(const RiceFeatures& input) {
    RiceFeatures clean;
    
    double minArea, maxArea;
    getIQR(input.circularities, minArea, maxArea);
    
    size_t count = input.circularities.size();
    int removed = 0;

    for (size_t i = 0; i < count; ++i) {
        double area = input.circularities[i];
        
        // 2. Check if this specific grain is an outlier
        if (area >= minArea && area <= maxArea) {
            // KEEP IT: Copy data to the clean struct
            clean.polyAreas.push_back(input.polyAreas[i]);
            clean.polyPerimeters.push_back(input.polyPerimeters[i]);
            clean.digitalAreas.push_back(input.digitalAreas[i]);
            clean.digitalPerimeters.push_back(input.digitalPerimeters[i]);
            clean.circularities.push_back(input.circularities[i]);
        } else {
            removed++;
        }
    }
    
    std::cout << "Outlier Removal: Kept " << clean.polyAreas.size() 
              << " grains. Removed " << removed << " outliers." << std::endl;
              
    return clean;
}

int main(int argc, char** argv)
{
    Board2D aBoard;
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings

    //typedef Object<DT8_4, DigitalSet> ObjectType; // Digital object with (8,4)-adjacency pair
    typedef ImageSelector<Domain, unsigned char >::Type Image; // type of image
    // 1) read an image
    Image image = PGMReader<Image>::importPGM ("../RiceGrains/Rice_camargue_seg_bin.pgm"); // Change based on desired image

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

        if (!touchesBoundary && obj.pointSet().size() > 20) { 
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

    Z2i::Curve c = getBoundary(completeObjects[99]);

    // DISPLAY STEP 3

    // aBoard << completeObjects[0];
    // aBoard << c;

    // complete the function "getBoundary" written above and add/modify your codes
    // STEP 4
    Range r; 
    auto curvePoints = c.getPointsRange(); 
    std::copy(curvePoints.begin(), curvePoints.end(), std::back_inserter(r));
    SegmentComputer recognitionAlgorithm;
    Segmentation seg(r.begin(), r.end(), recognitionAlgorithm);
        
    Segmentation::SegmentComputerIterator seg_i = seg.begin();

    // DISPLAY STEP 4
    aBoard << completeObjects[99];
    aBoard << c;
    for (seg_i = seg.begin(); seg_i != seg.end(); ++seg_i) {
        // This holds the actual shape (ArithmeticalDSS) that supports "BoundingBox"
        SegmentComputer::Primitive currentShape = seg_i->primitive();
        
        // Draw the Bounding Box using the Shape object
        aBoard << SetMode(currentShape.className(), "BoundingBox")
               << CustomStyle(currentShape.className() + "/BoundingBox", new CustomPen(Color::Red, Color::None, 2.0))
               << currentShape;
    }

    // STEP 5,6,7
    RiceFeatures features = analyzeRiceGrains(completeObjects);

    auto computeStats = [](const std::vector<double>& data, const std::string& name) {
        if (data.empty()) {
            std::cout << name << ": No valid data." << std::endl;
            return;
        }
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
        double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / data.size() - mean * mean);
        std::cout << name << " (Mean/Std): " << mean << " / " << stdev << std::endl;
    };

    computeStats(features.digitalAreas, "Digital Area");
    computeStats(features.polyAreas, "Area");
    computeStats(features.digitalPerimeters, "Digital Perimeter");
    computeStats(features.polyPerimeters, "Perimeter");
    computeStats(features.circularities, "Circularity");

    // STEP 8
    // exportToCSV(features, "basmati.csv"); // Mutliple images can be analyzed by changing the input image path
    RiceFeatures cleanFeatures = cleanOutliers(features);
    exportToCSV(cleanFeatures, "camargue_cleaned.csv");


    aBoard.saveSVG("./output.svg", 600, 600);

    return 0;
}
