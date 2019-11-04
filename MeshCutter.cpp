#include "MeshCutter.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////
void MeshCutter::init()
{
    setSelectRegionWidth(10);
    setSelectRegionHeight(10);
}
////////////////////////////////////////////////////////////////////////////////

EdgePtr Mesh:: addEdge( NodePtr &n0, NodePtr &n1, FacePtr &face)
{
    NodePtr vmin = std::min(n0,n1);

    for( auto oldedge : vmin->edges) {
        if( oldedge->hasNodes(n0,n1)) {
            oldedge->faces[1] = face;
            return oldedge;
        }
    }

    EdgePtr newedge(new Edge(n0,n1));
    newedge->faces[0] = face;
    vmin->edges.push_back(newedge);
    edges.push_back(newedge);
    return newedge;
}

////////////////////////////////////////////////////////////////////////////////

void Mesh:: makeConsistent( FacePtr &f)
{
    if( !f->active) return; 

    int nCount = 0;

    NodePtr steiner[3];
    for( int i = 0; i < 3; i++) {
        if(f->edges[i]->steinerNode ) nCount++;
    }
    if( nCount == 0) return;
    cout << "Consistent " << nCount << endl;

    NodePtr nodes[5];
    if( nCount == 1) {
        for( int i = 0; i < 3; i++) {
            if( f->edges[i]->steinerNode ) {
                nodes[0] = f->nodes[i];
                nodes[1] = f->nodes[(i+1)%3];
                nodes[2] = f->nodes[(i+2)%3];
                nodes[3] = f->edges[i]->steinerNode;
            }
        }
        remove(f);
        FacePtr f0 = Face::newObject(nodes[0], nodes[3], nodes[2]);
        FacePtr f1 = Face::newObject(nodes[1], nodes[2], nodes[3]);
	addFace(f0);
	addFace(f1);
        return;
    }

    if( nCount == 2) {
        for( int i = 0; i < 3; i++) {
            if( f->edges[i]->steinerNode == nullptr) {
                nodes[0] = f->nodes[i];
                nodes[1] = f->nodes[(i+1)%3];
                nodes[2] = f->nodes[(i+2)%3];
                nodes[3] = f->edges[(i+1)%3]->steinerNode;
                nodes[4] = f->edges[(i+2)%3]->steinerNode;
                break;
            }
        }
        remove(f);
        FacePtr f0 = Face::newObject( nodes[2], nodes[4], nodes[3]);
        FacePtr f1 = Face::newObject( nodes[0], nodes[1], nodes[3]);
        FacePtr f2 = Face::newObject( nodes[0], nodes[3], nodes[4]);
	addFace(f0);
	addFace(f1);
	addFace(f2);
        return;
    }
}

void MeshCutter:: readMesh( const string &filename)
{
    std::vector<float>  nodes;
    std::vector<size_t> faces;

    ifstream ifile( filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Warning: Input file not read " << endl;
        return;
    }
    string str;
    ifile >> str;
    if( str != "OFF") {
        cout << "Warning: Input file not in Off format" << endl;
        return;
    }

    size_t numNodes, numFaces, numEdges;
    ifile >> numNodes >> numFaces >> numEdges;

    double minval[3], maxval[3];

    minval[0] = std::numeric_limits<double>::max();
    minval[1] = std::numeric_limits<double>::max();
    minval[2] = std::numeric_limits<double>::max();

    maxval[0] =-std::numeric_limits<double>::max();
    maxval[1] =-std::numeric_limits<double>::max();
    maxval[2] =-std::numeric_limits<double>::max();

    double x, y, z;
    nodes.resize(3*numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        ifile >> x >> y >> z;
        minval[0] = min( minval[0], x);
        minval[1] = min( minval[1], y);
        minval[2] = min( minval[2], z);

        maxval[0] = max( maxval[0], x);
        maxval[1] = max( maxval[1], y);
        maxval[2] = max( maxval[2], z);

        nodes[3*i+0] = x;
        nodes[3*i+1] = y;
        nodes[3*i+2] = z;
    }

    double xc = 0.5*(maxval[0] + minval[0]);
    double yc = 0.5*(maxval[1] + minval[1]);
    double zc = 0.5*(maxval[2] + minval[2]);

    double xlen = maxval[0] - minval[0];
    double ylen = maxval[1] - minval[1];
    double zlen = maxval[2] - minval[2];

    mesh.center[0] = xc;
    mesh.center[1] = yc;
    mesh.center[2] = zc;
    mesh.radius    = sqrt(xlen*xlen + ylen*ylen + zlen*zlen);

    mesh.nodes.resize(numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        NodePtr v = Node::newObject();
        v->xyz[0] = nodes[3*i];
        v->xyz[1] = nodes[3*i+1];
        v->xyz[2] = nodes[3*i+2];
        v->id     = i;
        mesh.nodes[i] = v;
    }

    faces.resize(3*numFaces);
    size_t index = 0;
    int dummy, v0, v1, v2;

    for( size_t i = 0; i < numFaces; i++) {
        ifile >> dummy >> v0 >> v1 >> v2;
        assert( dummy == 3);
        faces[3*i+0] = v0;
        faces[3*i+1] = v1;
        faces[3*i+2] = v2;
        NodePtr n0   = mesh.nodes[v0];
        NodePtr n1   = mesh.nodes[v1];
        NodePtr n2   = mesh.nodes[v2];
        FacePtr newface = Face::newObject(n0,n1,n2);
        newface->id     = i;
        mesh.addFace(newface);
    }

    geomesh.initialize_mesh_data(nodes, faces);

    algorithm = new geodesic::GeodesicAlgorithmDijkstra(&geomesh);
    assert(algorithm);

    ANNpointArray  dataPts;
    dataPts = annAllocPts(numNodes, 3);
    for( int i = 0; i < numNodes; i++) {
        for( int j = 0; j < 3; j++) {
            dataPts[i][j] = nodes[3*i+j];
        }
    }
    kdtree = new ANNkd_tree(dataPts, numNodes, 3);
    assert( kdtree);

    ifile.close();

    ifile.open("landmarks.dat", ios::in);
    if( ifile.fail() ) return;

    for( auto v: mesh.nodes) v->visit = 0;

    int nContours;

    ANNcoord qpoint[3];
    double   eps = 0.0;
    ANNidx   nnId;
    ANNdist  dist;

    ifile >> nContours;

    int npoints;
    for( int j = 0; j < nContours; j++) {
        ifile >> npoints; assert( npoints >= 3);
        newContour.clear();
        newContour.offSurfPoints.resize(npoints);
        for( int i = 0; i < npoints; i++)  {
            ifile >> qpoint[0] >> qpoint[1] >> qpoint[2];
            newContour.offSurfPoints[i][0] = qpoint[0];
            newContour.offSurfPoints[i][1] = qpoint[1];
            newContour.offSurfPoints[i][2] = qpoint[2];
            kdtree->annkSearch(qpoint, 1, &nnId, &dist);
            if( mesh.nodes[nnId]->visit == 0) {
                newContour.landmarks.push_back(mesh.nodes[nnId]);
                mesh.nodes[nnId]->visit = 1;
            }
        }
        getLandMarksPath();
        contours.push_back(newContour);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Mesh::addFace( FacePtr &newface)
{
    assert( newface );
    auto n0 = newface->nodes[0]; assert( n0 );
    auto n1 = newface->nodes[1]; assert( n1 );
    auto n2 = newface->nodes[2]; assert( n2 );
    assert((n0 != n1) && (n1 != n2) && (n2 != n0));

    newface->edges[0] = addEdge(n0,n1,newface);
    newface->edges[1] = addEdge(n1,n2,newface);
    newface->edges[2] = addEdge(n2,n0,newface);

    n0->faces.push_back(newface);
    n1->faces.push_back(newface);
    n2->faces.push_back(newface);

    faces.push_back(newface);
}

////////////////////////////////////////////////////////////////////////////////
void Mesh::remove( FacePtr &oldface)
{

    for( int i = 0; i < 3; i++) {
        auto e = oldface->edges[i];
        if( e->faces[0] == oldface) {
            e->faces[0] = e->faces[1];
            e->faces[1] = nullptr;
        }
        if( e->faces[1] == oldface) {
            e->faces[1] = nullptr;
        }
    }

    for( int i = 0; i < 3; i++) {
        auto v  = oldface->nodes[i];
        if( !v->faces.empty() ) {
            auto it = std::remove(v->faces.begin(), v->faces.end(), oldface);
            v->faces.erase(it, v->faces.end() );
        }
    }

    oldface->active = 0;

}
////////////////////////////////////////////////////////////////////////////////
void Mesh::flip( EdgePtr &e)
{

    if( !e->active ) return;
    if( e->faces[0] == nullptr || e->faces[0] == nullptr) return;

    auto n0  = e->nodes[0];
    auto n1  = e->nodes[1];
    auto f0  = e->faces[0];
    auto f1  = e->faces[1];
    auto on0 = f0->getOpposite( n0, n1);
    auto on1 = f1->getOpposite( n0, n1);
    double d0 = length2(n0,n1);
    double d1 = length2(on0,on1);

    if( d0 < d1) return;

    double theta0 = f0->getAngleAt(on0);
    double theta1 = f1->getAngleAt(on1);

    if( theta0 + theta1 < 180) return;

    remove(f0);
    remove(f1);

    f0 = Face::newObject(n0, on0, on1);
    f1 = Face::newObject(n1, on1, on0);

    addFace(f0);
    addFace(f1);

    e->active = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Mesh::refine(FacePtr &f, int type)
{
    assert(f);

    if(f->active) {
        if(f->getArea() > 1.0E-10) {
            remove(f);

            if( type == 14) {
            NodePtr steiner[3];
            for( int i = 0; i < 3; i++) {
                auto edge = f->edges[i];
                if( edge->steinerNode == nullptr) {
                    auto newnode = Node::newObject();
                    newnode->xyz = edge->getCenter();
                    newnode->id  = nodes.size();
                    nodes.push_back(newnode);
                    edge->steinerNode = newnode;
                }
                steiner[i] = edge->steinerNode;
                assert( steiner[i] );
            }
            FacePtr f0 = Face::newObject( f->nodes[0], steiner[0], steiner[2]);
            FacePtr f1 = Face::newObject( f->nodes[1], steiner[1], steiner[0]);
            FacePtr f2 = Face::newObject( f->nodes[2], steiner[2], steiner[1]);
            FacePtr f3 = Face::newObject( steiner[0],  steiner[1], steiner[2]);
            addFace(f0);
            addFace(f1);
            addFace(f2);
            addFace(f3);
	    }

	    if( type == 13) {
                NodePtr newnode = Node::newObject();
                newnode->xyz    = f->getCentroid();
                newnode->id     = nodes.size();
                nodes.push_back(newnode);

                FacePtr f0 = Face::newObject( f->nodes[0], f->nodes[1], newnode);
                FacePtr f1 = Face::newObject( f->nodes[1], f->nodes[2], newnode);
                FacePtr f2 = Face::newObject( f->nodes[2], f->nodes[0], newnode);
                addFace(f0);
                addFace(f1);
                addFace(f2);
	    }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Mesh::saveAs( const std::string &filename)
{
    ofstream ofile(filename.c_str(), ios::out);

    size_t  numnodes = 0;
    for( auto v: nodes) {
        if( v->active ) v->id = numnodes++;
    }

    set<NodePtr> vSet;
    size_t  numfaces = 0;
    for( auto f: faces) {
        if( f->active ) {
            vSet.insert(f->nodes[0]);
            vSet.insert(f->nodes[1]);
            vSet.insert(f->nodes[2]);
            numfaces++;
        }
    }

    size_t index = 0;
    for( auto v: vSet) {
        v->id = index++;
        ofile << "v " << v->xyz[0] << " " << v->xyz[1] << " " << v->xyz[2] << endl;
    }

    for( auto f: faces) {
        if( f->active )
            ofile << "f " << f->nodes[0]->id +1 << " "
                  << f->nodes[1]->id +1 << " "
                  << f->nodes[2]->id +1 << endl;
    }

    cout << "Info: Cutout written to file: " << filename << endl;
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter::keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_0) {
        pickEntity = 0;
        this->setSelectedName(-1);
    }

    if( e->key() == Qt::Key_2) {
        newContour.clear();
        pickEntity = 2;
        this->setSelectedName(-1);
    }

    if( e->key() == Qt::Key_D) {
        refinePath();
        update();
        return;
    }

    if( e->key() == Qt::Key_F) {
        mesh.saveAs("cutout.obj");
        update();
        return;
    }

    if( e->key() == Qt::Key_H) {
        removeHole();
        return;
    }

    if( e->key() == Qt::Key_N) {
        displayIDs = !displayIDs;
    }

    if( e->key() == Qt::Key_R) {
        for( auto v : mesh.nodes) v->active = 1;
        for( auto e : mesh.edges) {
            e->active = 1;
            e->interface = 0;
        }
        for( auto f : mesh.faces) f->active = 1;
        contours.clear();
        newContour.clear();
    }

    if( e->key() == Qt::Key_S) {
        displaySurface = !displaySurface;
        update();
        return;
    }

    if( e->key() == Qt::Key_X) {
        if( !contours.empty() ) {
            contours.back().clear();
            contours.pop_back();
        }
    }

    if( e->key() == Qt::Key_W) {
        displayWires = !displayWires;
    }

    if( e->key() == Qt::Key_L) {
        useLights = !useLights;
        update();
        return;
    }

    if( e->key() == Qt::Key_X) {
        if( !newContour.boundnodes.empty() ) newContour.boundnodes.pop_back();
        update();
        return;
    }

    if( e->key() == Qt::Key_O) {
       for( int i = 0; i < contours.size(); i++)
	       contours[i].smooth();
        update();
        return;

    }

    if( e->key() == Qt::Key_Home) {
        qglviewer::Vec pos;
        pos[0]  = mesh.center[0];
        pos[1]  = mesh.center[1];
        pos[2]  = mesh.center[2];
        camera()->setSceneCenter(pos);
        camera()->setSceneRadius(mesh.radius);
        camera()->centerScene();
        camera()->showEntireScene();
        update();
        return;
    }

    QGLViewer::keyPressEvent(e);

    update();
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter:: mousePressEvent( QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter:: setContourEdges()
{
    int nSize = newContour.boundnodes.size();
    for( int i = 0; i < nSize; i++) {
        auto n0 = newContour.boundnodes[i];
        auto n1 = newContour.boundnodes[(i+1)%nSize];
        auto vm = min(n0,n1);
        for( auto e: vm->edges ) {
            if( e->hasNodes(n0,n1) ) {
                e->interface = 1;
                newContour.boundedges.push_back(e);
                break;
            } else {
                cout << "Warning: No edge found " << endl;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter:: getLandMarksPath()
{
    cout << "Connecting landmarks ..." << endl;

    int nsize = newContour.landmarks.size();
    assert( nsize );
    newContour.boundnodes.clear();
    newContour.boundedges.clear();

    ANNcoord qpoint[3];
    double   eps = 0.0;
    ANNidx   nnId;
    ANNdist  dist;

    vector<int> pathnodes;

    std::vector<geodesic::SurfacePoint> path;  //geodesic path is a sequence of SurfacePoints
    for( int i = 0; i < nsize; i++) {
        int srcid = newContour.landmarks[i]->id;
        int dstid = newContour.landmarks[(i+1)%nsize]->id;
        geodesic::SurfacePoint source(&geomesh.vertices()[srcid]);
        geodesic::SurfacePoint target(&geomesh.vertices()[dstid]);
        path.clear();
        algorithm->geodesic(source, target, path); //find a single source-target path
        std::reverse(path.begin(), path.end());

        pathnodes.clear();
        for( int j = 0; j < path.size(); j++) {
            geodesic::SurfacePoint& s = path[j];
            qpoint[0] = s.x();
            qpoint[1] = s.y();
            qpoint[2] = s.z();
            kdtree->annkSearch(qpoint, 1, &nnId, &dist);
            assert( fabs(dist) < 1.0E-15);
            pathnodes.push_back(nnId);
        }

        assert( pathnodes.front() == srcid);
        assert( pathnodes.back()  == dstid);

        for( int j = 0; j < path.size()-1; j++) {
            auto n0 = mesh.nodes[pathnodes[j]];
            auto n1 = mesh.nodes[pathnodes[j+1]];
            auto vm = min(n0,n1);
            for( auto e: vm->edges ) {
                if( e->hasNodes(n0,n1) ) {
                    e->interface = 1;
                    newContour.boundedges.push_back(e);
                    break;
                }
            }
        }
        for( auto e : newContour.boundedges)
            newContour.boundnodes.push_back( e->nodes[0] );
    }
    cout << "Bounded Edges  " << newContour.boundedges.size() << endl;
}

///////////////////////////////////////////////////////////////////////////////

void MeshCutter:: refinePath()
{
    set<NodePtr> vset;
    for( int i = 0; i < contours.size(); i++) {
        for( auto e: contours[i].boundedges) {
            vset.insert(e->nodes[0] );
            vset.insert(e->nodes[1] );
        }
    }

    set<FacePtr> fset;
    for( auto v: vset) {
        for( auto f: v->faces) fset.insert(f);
    }
    vset.clear();

    set<EdgePtr> eset;
    for( auto f: fset) {
        eset.insert( f->edges[0] );
        eset.insert( f->edges[1] );
        eset.insert( f->edges[2] );
    }

    for( auto f: fset) mesh.refine(f);
    fset.clear();

//  for( auto e: eset) mesh.flip(e);
    eset.clear();

    size_t numFaces = mesh.faces.size();
    for( int i = 0; i < numFaces; i++) mesh.makeConsistent( mesh.faces[i] );

    for( auto v: mesh.nodes) v->visit = 0;
    for( auto e: mesh.edges) e->interface = 0;

    int numnodes =  mesh.nodes.size();
    std::vector<float>  nodes;
    nodes.reserve(3*numnodes);

    for( auto v: mesh.nodes) {
        if( v->active ) {
            nodes.push_back( v->xyz[0] );
            nodes.push_back( v->xyz[1] );
            nodes.push_back( v->xyz[2] );
        }
    }

    int numfaces =  mesh.faces.size();
    std::vector<size_t> faces;
    faces.reserve(3*numfaces);

    for( auto f: mesh.faces) {
        if( f->active) {
            faces.push_back(f->nodes[0]->id);
            faces.push_back(f->nodes[1]->id);
            faces.push_back(f->nodes[2]->id);
        }
    }

    geomesh.initialize_mesh_data(nodes, faces);

    if( algorithm ) delete algorithm;

    algorithm = new geodesic::GeodesicAlgorithmDijkstra(&geomesh);
    assert(algorithm);

    ANNpointArray  dataPts;
    dataPts = annAllocPts(numnodes, 3);
    for( int i = 0; i < numnodes; i++) {
        for( int j = 0; j < 3; j++) {
            dataPts[i][j] = nodes[3*i+j];
        }
    }

    if( kdtree ) delete kdtree;
    kdtree = new ANNkd_tree(dataPts, numnodes, 3);
    assert( kdtree);


    ANNcoord qpoint[3];
    double   eps = 0.0;
    ANNidx   nnId;
    ANNdist  dist;

    for( int j = 0; j < contours.size(); j++) {
        newContour = contours[j];
        newContour.landmarks.clear();
        int npoints = newContour.offSurfPoints.size();
        for( int i = 0; i < npoints; i++)  {
            qpoint[0] = newContour.offSurfPoints[i][0];
            qpoint[1] = newContour.offSurfPoints[i][1];
            qpoint[2] = newContour.offSurfPoints[i][2];
            kdtree->annkSearch(qpoint, 1, &nnId, &dist);
            if( mesh.nodes[nnId]->visit == 0) {
                newContour.landmarks.push_back(mesh.nodes[nnId]);
                mesh.nodes[nnId]->visit = 1;
            }
        }
        getLandMarksPath();
    }

}
///////////////////////////////////////////////////////////////////////////////

void MeshCutter:: selectNode( size_t id)
{
    if( pickEntity != 0) return;

    if( newContour.boundnodes.empty() ) {
        newContour.boundnodes.push_back(mesh.nodes[id]);
        return;
    }

    auto lastnode = newContour.boundnodes.back();
    int srcid = lastnode->id;
    geodesic::SurfacePoint source(&geomesh.vertices()[srcid]);
    geodesic::SurfacePoint target(&geomesh.vertices()[id]);

    std::vector<geodesic::SurfacePoint> path;  //geodesic path is a sequence of SurfacePoints
    algorithm->geodesic(source, target, path); //find a single source-target path

    std::reverse(path.begin(), path.end());

    ANNcoord qpoint[3];
    double   eps = 0.0;
    ANNidx   nnId;
    ANNdist  dist;

    if( path.size() > 2) {
        for( int i = 1; i <path.size()-1; i++) {
            geodesic::SurfacePoint& s = path[i];
            qpoint[0] = s.x();
            qpoint[1] = s.y();
            qpoint[2] = s.z();
            kdtree->annkSearch(qpoint, 1, &nnId, &dist, eps);
            newContour.boundnodes.push_back(mesh.nodes[nnId]);
        }
    }

    if( newContour.boundnodes[0]->id == id) {
        setContourEdges();
        contours.push_back(newContour);
        newContour.clear();
        setSelectedName(-1);
        pickEntity = 2;
        return;
    }
    newContour.boundnodes.push_back(mesh.nodes[id]);
}
////////////////////////////////////////////////////////////////////////////////

void MeshCutter::mouseReleaseEvent( QMouseEvent *e)
{
    int id = this->selectedName();

    if( id >= 0) selectNode(id);

    QGLViewer::mouseReleaseEvent(e);

    update();
}

////////////////////////////////////////////////////////////////////////////////
void MeshCutter::drawNodes()
{
    glDisable( GL_LIGHTING);
    glPointSize(2);
    glColor3f( 0.0, 0.0, 1.0);

    /*
    glBegin(GL_POINTS);
    for( auto v: mesh.nodes) {
        if(v->active) glVertex3fv( &v->xyz[0]);
    }
    glEnd();

    if( !contours.empty() ) {
        glPointSize(5);
        glColor3f( 1.0, 0.0, 0.0);
        glBegin(GL_POINTS);
        for ( size_t i = 0; i < contours.size(); i++) {
            for ( auto v: contours[i].boundnodes) {
                glVertex3fv( &v->xyz[0] );
            }
        }
        glEnd();
    }

    if( !newContour.boundnodes.empty() ) {
        glPointSize(5);
        glColor3f( 1.0, 0.0, 0.0);
        glBegin(GL_POINTS);
        for ( auto v: newContour.boundnodes) {
            glVertex3fv( &v->xyz[0] );
        }
        glEnd();
    }

    glPointSize(5);
    glColor3f( 1.0, 1.0, 0.0);
    glBegin(GL_POINTS);
    for ( size_t j = 0; j < contours.size(); j++) {
        for( auto v: contours[j].landmarks)
            glVertex3fv( &v->xyz[0] );
    }
    */

    glPointSize(10);
    glColor3f( 0.5, 1.0, 1.0);
    glBegin(GL_POINTS);
    for ( size_t j = 0; j < contours.size(); j++) {
        for( auto v: contours[j].offSurfPoints)
            glVertex3fv( &v[0] );
    }

    glEnd();

}

////////////////////////////////////////////////////////////////////////////////
//
void MeshCutter:: removeHole()
{
    if( pickEntity == 0) return;

    int id = this->selectedName();
    if( id < 0) {
        cout << "Warning: No Seed face selected: Press 2 and select a face" << endl;
        return;
    }

    auto seedface = mesh.faces[id];

    for( auto f : mesh.faces ) f->visited = 0;

    deque<FacePtr> faceQ;
    faceQ.push_back( seedface );

    while( !faceQ.empty() ) {
        auto currface = faceQ.front(); faceQ.pop_front();
        if( !currface->visited ) {
            currface->visited = 1;
            newContour.internalfaces.push_back(currface);
            for( int i = 0; i < 3; i++) {
                auto edge = currface->edges[i];
                if( !edge->interface ) {
                    for( int j = 0; j < 2;  j++) {
                        auto nextface = edge->faces[j];
                        if( nextface )  {
                            if( nextface->visited == 0) faceQ.push_back(nextface);
                        }
                    }
                }
            }
        }
    }

    set<EdgePtr> alledges;
    set<NodePtr> allnodes;
    for( auto f : newContour.internalfaces)  {
        f->active = 0;
        for(int i = 0; i < 3; i++) {
            alledges.insert(f->edges[i]);
            allnodes.insert(f->nodes[i]);
        }
    }
    set<NodePtr> boundnodes;
    for( auto v : newContour.boundnodes) boundnodes.insert(v);

    set<NodePtr> internalnodes;
    std::set_difference( allnodes.begin(), allnodes.end(),
                         boundnodes.begin(), boundnodes.end(),
                         std::inserter(internalnodes, internalnodes.end()));


    for( auto v: internalnodes) v->active = 0;

    set<EdgePtr> boundedges;
    for( auto e : newContour.boundedges) boundedges.insert(e);

    set<EdgePtr> internaledges;
    std::set_difference( alledges.begin(), alledges.end(),
                         boundedges.begin(), boundedges.end(),
                         std::inserter(internaledges, internaledges.end()));

    for( auto e: internaledges) e->active = 0;

    pickEntity = 0;
    setSelectedName(-1);

    update();
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter::drawContours()
{
    glDisable( GL_LIGHTING);
    glLineWidth(5);
    glColor3f( 1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    for ( size_t i = 0; i < contours.size(); i++) {
        for ( auto e : contours[i].boundedges) {
            auto v0 = e->nodes[0];
            auto v1 = e->nodes[1];
            glVertex3fv( &v0->xyz[0] );
            glVertex3fv( &v1->xyz[0] );
        }
    }
    glEnd();
    glLineWidth(1);
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter::drawFaces(int filled)
{
    if( useLights ) glEnable(GL_LIGHTING);

    if( filled ) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f( 0.8, 0.9, 0.9);
        for ( auto f : mesh.faces) {
            if(f->active) {
                glBegin(GL_TRIANGLES);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
            }
        }
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f( 0.2, 0.2, 0.2);
        for ( auto f : mesh.faces) {
            if(f->active) {
                glBegin(GL_LINE_LOOP);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
            }
        }
    }



    if( pickEntity == 2) {
        int id = this->selectedName();
        if( id >= 0) {
            auto f = mesh.faces[id];
            if(f->active) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glColor3f( 0.0, 1.0, 0.0);
                glBegin(GL_TRIANGLES);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void MeshCutter:: drawIDs()
{
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    startScreenCoordinatesSystem();

    qglviewer::Vec pos, proj;

    for( int i = 0; i < contours.size(); i++) {
        int index = 0;
        for( auto v : contours[i].landmarks) {
            pos[0] = v->xyz[0];
            pos[1] = v->xyz[1];
            pos[2] = v->xyz[2];
            proj = camera()->projectedCoordinatesOf(pos);
            QFont font("Times", 12);
            drawText(int(proj.x), int(proj.y), QString::number(index), font);
            index++;
        }
    }
    stopScreenCoordinatesSystem();

    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
}
////////////////////////////////////////////////////////////////////////////////

void MeshCutter::drawWithNames()
{
    glDisable( GL_LIGHTING);
    glPointSize(2);
    glColor3f( 0.0, 0.0, 1.0);

    if( pickEntity == 0) {
        for ( auto v: mesh.nodes) {
            if( v->active) {
                glPushName(v->id);
                glBegin(GL_POINTS);
                glVertex3fv( &v->xyz[0] );
                glEnd();
                glPopName();
            }
        }
    }

    if( pickEntity == 2) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f( 0.8, 0.9, 0.9);
        for ( auto f: mesh.faces) {
            if( f->active) {
                glPushName(f->id);
                glBegin(GL_TRIANGLES);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
                glPopName();
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter::draw()
{
    glPolygonOffset(1.0,1.0);
    glEnable(GL_POLYGON_OFFSET_FILL);

    drawNodes();
    drawContours();
    if( displayWires)   drawFaces(0);
    if( displaySurface) drawFaces(1);
    if( displayIDs)     drawIDs();
}

////////////////////////////////////////////////////////////////////////////////
