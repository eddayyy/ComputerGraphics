#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <chrono>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>
#include <GL/gl.h>



//------------------------------------------------------------------------------------------------------------------------
// GLFW and GLUT windows
//
const GLint WIDTH = 800, HEIGHT = 600;

int glfw_open_window() {    // new window with GLFW
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
  GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Learn OpenGL using GLFW window", nullptr, nullptr);
  int screenWidth, screenHeight;
  glfwGetFramebufferSize(window, &screenWidth, &screenHeight);
  if (window == nullptr) {
    std::cout << "Failed to create GLFW window\n";
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {  std::cout << "Failed to initialize GLEW\n";  return -1;  }

  glViewport(0, 0, screenWidth, screenHeight);
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();  glClearColor(0.8f, 0.3f, 0.3f, 1.0f);  glClear(GL_COLOR_BUFFER_BIT); glfwSwapBuffers(window);
  }
  glfwTerminate();
  return 0;
}


void RenderSceneCB() {  glClear(GL_COLOR_BUFFER_BIT);  glutSwapBuffers();  }
void InitializeGlutCallbacks() { glutDisplayFunc(RenderSceneCB); }

int glut_open_window(int argc, char* argv[]) {       // new window with GLUT
  glutInit(&argc, argv);  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);  glutInitWindowSize(1024, 768);  glutInitWindowPosition(100, 100);
  glutCreateWindow("Learn OpenGL using GLUT window");
  InitializeGlutCallbacks();  glClearColor(0.0f, 0.0f, 1.0f, 1.0f);  glutMainLoop();  return 0;
}
//
// end of GLFW and GLUT windows
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// display solid sphere with reflected light
//
void display_sphere_cool() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  glPushMatrix();
  GLfloat color[] = {1.0, 0.0, 0.0, 1.0};  glMaterialfv(GL_FRONT, GL_DIFFUSE, color);
  glutSolidSphere(2, 20, 20);  glPopMatrix();  glutSwapBuffers();
}

void glut_sphere_cool(int argc, char* argv[]) {
  glutInit(&argc, argv);  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);  glutInitWindowSize(400, 400);  glutCreateWindow("Sphere");
  glClearColor(0.0, 0.0, 0.0, 1.0);  glShadeModel(GL_SMOOTH);  glEnable(GL_CULL_FACE);  glEnable(GL_DEPTH_TEST);  glEnable(GL_LIGHTING);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};   GLfloat lightZeroColor[] = {0.8, 1.0, 0.8, 1.0};
                    // green tinged
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
  glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1);     glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05);  glEnable(GL_LIGHT0);
  glutDisplayFunc(display_sphere_cool);
  glMatrixMode(GL_PROJECTION);  gluPerspective(40.0, 1.0, 1.0, 40.0);  glMatrixMode(GL_MODELVIEW);  gluLookAt(0, 0, 10, 0, 0, 0, 0, 1, 0);
  glPushMatrix();  glutMainLoop();
}


//------------------------------------------------------------------------------------------------------------------------
// display segment through sphere
//
void display_sphere() {
  glClear(GL_COLOR_BUFFER_BIT);  glColor3f(1.0, 0.0, 0.0);  glLoadIdentity();
//  glutSolidSphere(5.0, 20.0, 20.0);
//  glutSolidSphere(5.0, 20, 20);
  glutSolidSphere(1.1, 20, 20);
//  glutSolidSphere(1.0, 20, 20);
//    glutSolidSphere(0.8, 80, 80);
  glFlush();
}

void myInit() {
  glClearColor(1.0, 1.0, 1.0, 1.0);  glColor3f(1.0, 0.0, 0.0);  glMatrixMode(GL_PROJECTION);  glLoadIdentity();
//  gluOrtho2D(0.0, 499.0, 0.0, 499.0);
//  gluOrtho2D(-5.0, 5.0, -5.0, 5.0);
  gluOrtho2D(-1.1, 1.1, -1.1, 1.1);
//  gluOrtho2D(-1.0, 1.0, -1.0, 1.0);

  glMatrixMode(GL_MODELVIEW);
}

void glut_sphere(int argc,char* argv[]) {
  GLUquadric* qobj = gluNewQuadric();   glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);  glutInitWindowSize(500, 500);  glutCreateWindow("pendulum");    glutDisplayFunc(display_sphere);
  myInit();
  glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
  glutMainLoop();
}


void resize(int width, int height) {
  const float ar = (float) width / (float) height;
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);  glLoadIdentity();  glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);
  glMatrixMode(GL_MODELVIEW);   glLoadIdentity() ;
}

void display_shapes() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColor3d(1, 0, 0);
  glPushMatrix();  glTranslated(0.0, 1.8, -6);                                      glutSolidSphere(1, 50, 50);       glPopMatrix();
  glPushMatrix();  glTranslated(-1.5, -1.2, -6);                                    glutWireSphere(1, 16, 16);        glPopMatrix();
  glPushMatrix();  glTranslated(-2.7, 1.2, -6);   glRotated(-20.0, 0.0, 0.0, 0.0);  glutWireIcosahedron();            glPopMatrix();
  glPushMatrix();  glTranslated(2.7, 1.2, -6);                                      glutWireTorus(0.4, 0.6, 12, 12);  glPopMatrix();
  glPushMatrix();  glTranslated(1.5, -1.2, -6);   glRotated(-30.0, 0.0, 0.0, 0.0);  glutWireTeapot(1.1);              glPopMatrix();
  glutSwapBuffers();
}

const GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };  const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };  const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f };
const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };  const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };  const GLfloat high_shininess[] = { 100.0f };


void glut_sphere_shaded(int argc, char* argv[]) {
  glutInit(&argc, argv);          glutInitWindowSize(640,480);
  glutInitWindowPosition(10,10);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

  glutCreateWindow("Programming Techniques - 3D Spheres");

  glutReshapeFunc(resize);
  glutDisplayFunc(display_shapes);

  glClearColor(1,1,1,1);
  glEnable(GL_CULL_FACE);  glCullFace(GL_BACK);  glEnable(GL_DEPTH_TEST);      glDepthFunc(GL_LESS);
  glEnable(GL_LIGHT0);  glEnable(GL_NORMALIZE);  glEnable(GL_COLOR_MATERIAL);  glEnable(GL_LIGHTING);

  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);   glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);  glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);

  glutMainLoop();
}


//===============================================================================
// rotating colored cube
//
//  This is the include statement I need for Mac OS X.
//
GLfloat vertices[][3] = {
  { -1.0, -1.0, -1.0 },
  {  1.0, -1.0, -1.0 },
  {  1.0,  1.0, -1.0 },
  { -1.0,  1.0, -1.0 },
  { -1.0, -1.0,  1.0 },
  {  1.0, -1.0,  1.0 },
  {  1.0,  1.0,  1.0 },
  { -1.0,  1.0,  1.0 } };

GLfloat normals[][3] = {
  { -1.0, -1.0, -1.0 },
  {  1.0, -1.0, -1.0 },
  {  1.0,  1.0, -1.0 },
  { -1.0,  1.0, -1.0 },
  { -1.0, -1.0,  1.0 },
  {  1.0, -1.0,  1.0 },
  {  1.0,  1.0,  1.0 },
  { -1.0,  1.0,  1.0 } };

GLfloat colors[][3] = {
  { 0.0, 0.0, 0.0 },
  { 1.0, 0.0, 0.0 },
  { 1.0, 1.0, 0.0 },
  { 0.0, 1.0, 0.0 },
  { 0.0, 0.0, 1.0 },
  { 1.0, 0.0, 1.0 },
  { 1.0, 1.0, 1.0 },
  { 0.0, 1.0, 1.0 } };

static GLint axis = 2;
static GLfloat theta[3] = { 0.0, 0.0, 0.0 };

//****************************************************************************80
//    Edward Angel
//    COLORCUBE defines the 6 faces of the color cube object.
//
//****************************************************************************80
//    POLYGON defines the colors, vertices and normals for a quadrilateral.
//
void polygon ( int a, int b, int c, int d ) {
  glBegin ( GL_POLYGON );

  glColor3fv ( colors[a] );   glNormal3fv ( normals[a] );   glVertex3fv ( vertices[a] );
  glColor3fv ( colors[b] );   glNormal3fv ( normals[b] );   glVertex3fv ( vertices[b] );
  glColor3fv ( colors[c] );   glNormal3fv ( normals[c] );   glVertex3fv ( vertices[c] );
  glColor3fv ( colors[d] );   glNormal3fv ( normals[d] );   glVertex3fv ( vertices[d] );

  glEnd ( );
}

void colorcube ( ) {
  polygon ( 0, 3, 2, 1 );      polygon ( 2, 3, 7, 6 );      polygon ( 0, 4, 7, 3 );
  polygon ( 1, 2, 6, 5 );      polygon ( 4, 5, 6, 7 );      polygon ( 0, 1, 5, 4 );
}

//****************************************************************************80
//    DISPLAY generates the graphics output.
//
void display ( ) {
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );   glLoadIdentity ( );
  glRotatef ( theta[0], 1.0, 0.0, 0.0 );   glRotatef ( theta[1], 0.0, 1.0, 0.0 );   glRotatef ( theta[2], 0.0, 0.0, 1.0 );
  colorcube();   glFlush ( );   glutSwapBuffers ( );
}

//****************************************************************************80
//    MOUSE determines the response to mouse input.
//
void mouse ( int btn, int state, int x, int y ) {
  std::cout << "Mouse function entered -- axis value is: " << axis << "\n";
  if ( btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) { ++axis; }
  if ( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN ) { ++axis; }
  if ( btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN ) { ++axis; }
  axis %= 3;
}

//****************************************************************************80
//    MYRESHAPE determines the window mapping.
//
void myReshape ( int w, int h ) {
  glViewport ( 0, 0, w, h );
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );

  if ( w <= h ) {
    glOrtho (
      -2.0, 2.0,
      -2.0 * ( GLfloat ) h / ( GLfloat ) w, 2.0 * ( GLfloat ) h / ( GLfloat ) w,
      -10.0, 10.0 );
  } else {
    glOrtho (
      -2.0 * ( GLfloat ) h / ( GLfloat ) w, 2.0 * ( GLfloat ) h / ( GLfloat ) w,
      -2.0, 2.0,
      -10.0, 10.0 );
  }
  glMatrixMode ( GL_MODELVIEW );
}

//****************************************************************************80
//    SPINCUBE adjusts the angle of rotation and redisplays the picture.
//
void spinCube ( ) {
  std::cout << "entered spinCube()...\n";
//  theta[axis] = theta[axis] + 0.020;
  theta[axis] += 1.0;
  if ( 360.0 < theta[axis] ) {
    theta[axis] = theta[axis] - 360.0;
  }
  glutPostRedisplay ( );
  std::cout << "\t\texiting spinCube()...\n\n";
}

//****************************************************************************80
//    This program constructs a cube, each of whose vertices is given a
//    different color, and displays the cube.  The cube rotates slowly
//    about the X, Y or Z axis.  Each time the user clicks the mouse, the
//    "next" axis is used for rotation.
//
void rotating_cube(int argc, char* argv[]) {
  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Rotating cube" );
  glutReshapeFunc ( myReshape );
  glutDisplayFunc ( display );
  glutIdleFunc ( spinCube );
  glutMouseFunc ( mouse );

//  Enable hidden surface removal.
  glEnable ( GL_DEPTH_TEST );
  glutMainLoop ( );

  std::cout << "axis: " << axis << "\n";
}


#define ColoredVertex(c,v) do{glColor3fv(c);glVertex3fv(v);}while(0)
#define WIDTH 400
#define HEIGHT 400

GLfloat angle = 0.0f;

float square_root(float x) { return (float) sqrt(x); }
//****************************************************************************80
//    This program constructs a cube, each of whose vertices is given a
//    different color, and displays the cube.  The cube rotates slowly
//    about the X, Y or Z axis.  Each time the user clicks the mouse, the
//    "next" axis is used for rotation.
//
void display_tetra(void) {
  static int list = 0;
  if (list == 0) {
    // If the display list does not exist, create
    GLfloat
      PointA[] = { 0.5f, -square_root(6.0f) / 12, -square_root(3.0f) / 6 },
      PointB[] = { -0.5f, -square_root(6.0f) / 12, -square_root(3.0f) / 6 },
      PointC[] = { 0.0f, -square_root(6.0f) / 12, square_root(3.0f) / 3 },
      PointD[] = { 0.0f, square_root(6.0f) / 4,0 };
    GLfloat
      ColorR[] = { 1, 0, 0 },
      ColorG[] = { 0, 1, 0 },
      ColorB[] = { 0, 0, 1 },
      ColorY[] = { 1, 1, 0 };
    list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glBegin(GL_TRIANGLES);
    //plane ABC
    ColoredVertex(ColorR,PointA);  ColoredVertex(ColorG, PointB);  ColoredVertex(ColorB, PointC);
    // flat ACD
    ColoredVertex(ColorR, PointA);  ColoredVertex(ColorB, PointC);  ColoredVertex(ColorY, PointD);
    // Plane CBD
    ColoredVertex(ColorB, PointC);  ColoredVertex(ColorG, PointB);  ColoredVertex(ColorY, PointD);
    // Plane BAD
    ColoredVertex(ColorG, PointB);  ColoredVertex(ColorR, PointA);  ColoredVertex(ColorY, PointD);
    glEnd();
    glEndList();
    glEnable(GL_DEPTH_TEST);
  }
  // A display list has been created, which will be called each time a regular tetrahedron is drawn
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();  glRotatef(angle, 1, 0.5, 0);  glCallList(list);  glPopMatrix();  glutSwapBuffers();
}

void idle(void) {
  ++angle;
  if (angle >= 360.0f) { angle = 0.0f; }
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  display_tetra();
}

void rotate_tetra(int argc, char* argv[]) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowPosition(200, 200);
  glutInitWindowSize(WIDTH, HEIGHT);
  glutCreateWindow("OpenGL Window");
  glutDisplayFunc(&display_tetra);

  glutIdleFunc(&idle);
  glutMainLoop();
}

int main ( int argc, char *argv[] ) {
  std::cout << "Starting program\n";
//  glut_sphere(argc, argv);
//  glut_sphere_cool(argc, argv);
//  glut_sphere_shaded(argc, argv);
//  rotating_cube(argc, argv);
//  glfw_open_window();
//  glut_open_window(argc, argv);
  rotate_tetra(argc, argv);

  return 0;
}






