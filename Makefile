OBJS = main.o MeshCutter.o

CPPFLAGS = -O3 -fPIC -g
CPPFLAGS += -I$(QGLVIEWER_DIR)/include
CPPFLAGS += -I$(ANN_DIR)/include
CPPFLAGS += -I$(QTDIR)/include -I$(QTDIR)/include/QtCore -I$(QTDIR)/include/QtWidgets -I$(QTDIR)/include/QtXml -I$(QTDIR)/include/QtOpenGL -I$(QTDIR)/include/QtGui

CPPFLAGS += -I$(GEODESIC_DIR)/include
CPPFLAGS += -I$(BOOST_DIR)/include

LIBS += -L$(QTDIR)/lib -lQt5Core -lQt5Xml -lQt5OpenGL -lQt5Widgets -lQt5Gui -lGL -lGLU
LIBS += -L$(QGLVIEWER_DIR)/lib -lQGLViewer
LIBS += -L$(ANN_DIR)/lib -lANN

meshcutter:$(OBJS)
	g++ -o meshcutter $(OBJS) $(LIBS)

.o:.cpp
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o meshcutter
