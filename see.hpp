#ifndef SEE_HPP
#define SEE_HPP

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#include <GL/freeglut.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "PBC.hpp"

// toofus
#include "AABB.hpp"
#include "ColorTable.hpp"
#include "fileTool.hpp"
#include "glTools.hpp"

PBC RefConf;
PBC Conf;
int refConfNum = 1;
int confNum = 1;

AABB worldBox;

int main_window;

// flags
int ref_fixed = 1;
int show_background = 1;
int show_particles = 1;
int show_ghosts = 1;
int show_displacements = 0;
int show_fluctuations = 0;
int show_cell = 1;
int show_forces = 0;
int showOrientations = 0;

int color_option = 0;
ColorTable colorTable;

GLfloat alpha_particles = 1.0f;
GLfloat alpha_ghosts = 0.3f;

double ghost_width = 0.05;
double arrowSize = 0.0005;
double arrowAngle = 0.7;
double vScale = 0.01;

double forceTubeFactor = 1.0;

double radiusMin;
double radiusMax;
double radiusMean;

int width = 800;
int height = 800;
float wh_ratio = (float)width / (float)height;
glTextZone textZone(1, &width, &height);

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int display_mode = 0;  // sample or slice rotation
int mouse_start[2];

// Drawing functions
void setColor(int i, GLfloat alpha);
void drawForces();
void drawFluctuations();
void drawDisplacements();
void drawBox();
void drawParticles();
void drawGhosts();

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void buildMenu();
void printHelp();
void fit_view();
bool try_to_readConf(int num, PBC& CF, int& OKNum);
void updateTextLine();
void add_ghost_pos(int i, double mn, double mx, std::vector<vec2r> & lst);

#endif /* end of include guard: SEE_HPP */
