This is a private project for dental segmentation project. 
We use openmesh to represent 3d mesh and use itk to deal with dicom data.

The main open source libraries we depend on are  OpenMesh, CGAL,QGLViewer, Boost, Eigen ,QT , OpenCV ,ann and itk.  All the codes should be compiled with MSVC 2015 on Windows. Other platforms and compilers are not guaranteed.

To compile our project with no trouble, please add some environment variables to your system (you may need replace the paths below to your actual ones):

TANN_DIR: D:\CodeLibrary\Ann1_1_2

TCGAL_BUILD: D:\CodeLibrary\CGAL4.9\winbuild64

TCGAL_DIR: D:\CodeLibrary\CGAL4.9\ When you cmake CGAL and generated Qt5 related solutions, set Qt5_DIR to ../YOURPATH/Qt/Qt5.7.0/5.7/msvc2015_64/lib/cmake/Qt5

TEIGEN_ROOT: D:\CodeLibrary\eigen3_3

TITK_BUILD: D:\CodeLibrary\InsightToolkit-4.10.1\Build

TITK_DIR: D:\CodeLibrary\InsightToolkit-4.10.1

TLIBIGL_DIR: D:\CodeLibrary\libigl

TLIBIGL_INCLDUEDIR: D:\CodeLibrary\libigl\include

TOPEN_MESH_DIR: D:\CodeLibrary\OpenMesh_6_3\

TOPENCV_DIR: D:\CodeLibrary\opencv310\build\install

TQGLVIEWER_DIR: D:\CodeLibrary\libQGLViewer-2.6.4_vs2015

TQTDIR: D:\CodeLibrary\QT_5_7\5.7\msvc2015_64

BOOST_INCLUDEDIR: D:\CodeLibrary\boost_1_62_0_vs2015

BOOST_LIBRARYDIR: D:\CodeLibrary\boost_1_62_0_vs2015\lib64-msvc-14.0



And further you may need to add your third-party dlls paths to your path in system variables. Alternatively, putting them together with the exe file is also a good choice.