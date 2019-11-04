#include <string>
#include <array>
#include <vector>
#include <deque>
#include <map>

#include <QGLViewer/qglviewer.h>
#include "geodesic_algorithm_dijkstra.h"
#include "geodesic_algorithm_exact.h"
#include <ANN/ANN.h>

#include <QKeyEvent>

class Node;
typedef std::shared_ptr<Node> NodePtr;

class Edge;
typedef std::shared_ptr<Edge> EdgePtr;

class Face;
typedef std::shared_ptr<Face> FacePtr;

struct Node
{
    static NodePtr newObject();
    int  id;
    bool active   = 1;
    bool boundary = 0;
    bool visit    = 0;
    std::array<float,3> xyz;
    std::vector<EdgePtr> edges;
    std::vector<FacePtr> faces;
};

inline double length2( const NodePtr n0, const NodePtr &n1)
{
    double dx = n0->xyz[0] - n1->xyz[0];
    double dy = n0->xyz[1] - n1->xyz[1];
    double dz = n0->xyz[2] - n1->xyz[2];

    double l2 = dx*dx + dy*dy + dz*dz;
    return l2;
}

inline NodePtr Node:: newObject()
{
    NodePtr v(new Node);
    return v;
}

struct Edge {
    static EdgePtr newObject( NodePtr &n0, NodePtr &n1);

    Edge() {};
    Edge( const NodePtr &v0, const NodePtr &v1) {
        nodes[0] = v0;
        nodes[1] = v1;
    }

    bool hasNodes( const NodePtr &n0, const NodePtr &n1) const {
        if( (nodes[0]->id == n0->id) && (nodes[1]->id == n1->id) ) return 1;
        if( (nodes[0]->id == n1->id) && (nodes[1]->id == n0->id) ) return 1;
        return 0;
    }

    bool isBoundary() const {
        if( faces[1] == nullptr) return 1;
    }

    std::array<float,3> getCenter() const
    {
	    std::array<float,3> p;
	    p[0] = 0.5*(nodes[0]->xyz[0] + nodes[1]->xyz[0]);
	    p[1] = 0.5*(nodes[0]->xyz[1] + nodes[1]->xyz[1]);
	    p[2] = 0.5*(nodes[0]->xyz[2] + nodes[1]->xyz[2]);
	    return p;
    }

    bool active    = 1;
    bool interface = 0;
    NodePtr steinerNode = nullptr;
    std::array<NodePtr,2>  nodes = {nullptr,nullptr};
    std::array<FacePtr,2>  faces = {nullptr,nullptr};
};

inline EdgePtr Edge:: newObject(NodePtr &n0, NodePtr &n1)
{
    EdgePtr e(new Edge(n0,n1));
    return e;
}

struct Face {
    static FacePtr newObject( NodePtr &n0, NodePtr &n1, NodePtr &n2);

    Face() {};
    Face( NodePtr &v0, NodePtr &v1, NodePtr &v2) {
        nodes[0] = v0;
        nodes[1] = v1;
        nodes[2] = v2;
    }

    NodePtr getOpposite( const NodePtr &n0, const NodePtr &n1) const
    {
        for( int i = 0; i < 3; i++) {
            if( nodes[i] == n0 && nodes[(i+1)%3] == n1) return nodes[(i+2)%3];
            if( nodes[i] == n1 && nodes[(i+1)%3] == n0) return nodes[(i+2)%3];
        }
        return nullptr;
    }

    int getPositionOf( const NodePtr &v) const
    {
        for( int i = 0; i < 3; i++) {
            if( nodes[i] == v ) return i;
        }
        return -1;
    }

    float getAngleAt( const NodePtr &v0) const;

    std::array<float,3> getCentroid() const;
    float getArea() const;

    bool active  = 1;
    bool visited = 0;
    int  id;
    std::array<NodePtr,3> nodes;
    std::array<EdgePtr,3> edges;
};

inline std::array<float,3> Face :: getCentroid() const
{
    std::array<float,3> center = {0.0, 0.0, 0.0};

    for( int i = 0; i < 3; i++) {
        auto p = nodes[i]->xyz;
        center[0] += p[0];
        center[1] += p[1];
        center[2] += p[2];
    }
    center[0] /= 3.0;
    center[1] /= 3.0;
    center[2] /= 3.0;

    return center;
}

inline float Face :: getArea() const
{
    float len[3];

    for( int i = 0; i < 3; i++) {
        auto p0 = nodes[(i+1)%3]->xyz;
        auto p1 = nodes[(i+2)%3]->xyz;
        float dx = p1[0] - p0[0];
        float dy = p1[1] - p0[1];
        float dz = p1[2] - p0[2];
        len[i]   = sqrt(dx*dx + dy*dy + dz*dz);
    }
    float a = len[0];
    float b = len[1];
    float c = len[2];
    float s = (a+b+c)/2.0;

    float ar = sqrt(s*(s-a)*(s-b)*(s-c));
    return ar;
}

inline float Face :: getAngleAt( const NodePtr &n0) const
{
    int pos = getPositionOf(n0);
    assert(pos >= 0);

    auto n1 = nodes[(pos+1)%3];
    auto n2 = nodes[(pos+2)%3];

    double a2 = length2(n1,n2);
    double b2 = length2(n2,n0);
    double c2 = length2(n0,n1);

    double cosA = (b2 + c2 - a2) /(2*sqrt(b2*c2));
    if( cosA > 1.0) cosA =  1.0;
    if( cosA <-1.0) cosA = -1.0;

    double A  = 180.0*acos( cosA )/M_PI;
    return A;
}

inline FacePtr Face:: newObject(NodePtr &n0, NodePtr &n1, NodePtr &n2)
{
    FacePtr f(new Face(n0,n1,n2));
    return f;
}

struct Mesh
{
    void refine( FacePtr &f, int type = 14);
    void flip(   EdgePtr &e);
    void makeConsistent( FacePtr &e);

    EdgePtr addEdge( NodePtr &n0, NodePtr &n1, FacePtr &f);
    void    addFace( FacePtr &f);
    void    remove(FacePtr &f);
    void    pruneAll();

    std::vector<NodePtr> nodes;
    std::vector<EdgePtr> edges;
    std::vector<FacePtr> faces;
    std::array<double,3> center = {0.0, 0.0, 0.0};

    double radius;
    void saveAs( const std::string &s);
};

struct Contour
{
    std::vector<std::array<float,3>>  offSurfPoints;
    std::vector<NodePtr>  landmarks;
    std::vector<NodePtr>  boundnodes;
    std::vector<EdgePtr>  boundedges;
    std::vector<FacePtr>  internalfaces;

    void clear() {
        for( auto e: boundedges) e->active = 1;
        offSurfPoints.clear();
        landmarks.clear();
        boundnodes.clear();
        boundedges.clear();
        internalfaces.clear();
    }

    void smooth()
    {
	    std::map<NodePtr,std::array<float,3>> vmap;

       std::array<float,3> pos = {0.0, 0.0, 0.0};

       for( int i = 0; i < boundedges.size(); i++) {
            auto v0 =  boundedges[i]->nodes[0];
            auto v1 =  boundedges[i]->nodes[1];
	    vmap[v0] = pos;
	    vmap[v1] = pos;
       }

       for( int i = 0; i < boundedges.size(); i++) {
            auto v0 =  boundedges[i]->nodes[0];
            auto v1 =  boundedges[i]->nodes[1];
	    auto p0 =  v0->xyz;
	    auto p1 =  v1->xyz;

	    vmap[v0][0] += p1[0];
	    vmap[v0][1] += p1[1];
	    vmap[v0][2] += p1[2];

	    vmap[v1][0] += p0[0];
	    vmap[v1][1] += p0[1];
	    vmap[v1][2] += p0[2];
       }
       for( auto keyval : vmap){
            auto v = keyval.first;
	    pos    = keyval.second;
	    pos[0] *= 0.5;
	    pos[1] *= 0.5;
	    pos[2] *= 0.5;
	    v->xyz = pos;
       }
    }
};

class MeshCutter : public QGLViewer
{
public:
    void readMesh( const std::string &s);

protected:
    virtual void draw();
    virtual void init();
    virtual void drawWithNames();
    virtual void keyPressEvent( QKeyEvent *e);
    virtual void mousePressEvent( QMouseEvent *e);
    virtual void mouseReleaseEvent( QMouseEvent *e);

private:
    std::vector<Contour> contours;

    Contour newContour;

    ANNkd_tree *kdtree = nullptr;

    geodesic::Mesh geomesh;
    geodesic::GeodesicAlgorithmDijkstra *algorithm = nullptr;

    Mesh mesh;

    int pickEntity = 0;
    bool displayWires = 1;
    bool displaySurface = 1;
    bool useLights      = 0;
    bool displayIDs     = 0;

    size_t nodeSelected;
    void removeHole();

    void drawIDs();
    void setContourEdges();
    void selectNode( size_t id);
    void drawNodes();
    void drawFaces(int filled);
    void drawContours();
    void getLandMarksPath();
    void refinePath();

    std::vector<int> getComponent(int fid);
};
