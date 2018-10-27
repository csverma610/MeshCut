#include "MeshCutter.h"
#include <qapplication.h>

int main(int argc, char **argv) {


  assert( argc == 2);
  QApplication application(argc, argv);

  // Instantiate the viewer.
  MeshCutter viewer;

  viewer.setWindowTitle("MeshCutter");
  viewer.readMesh( argv[1] );

  // Make the viewer window visible on screen.
  viewer.show();

  // Run main loop.
  return application.exec();
}
