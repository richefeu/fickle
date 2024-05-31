#include "see.hpp"

void printHelp() {
  using namespace std;
  cout << endl;
  cout << "+         load next configuration file" << endl;
  cout << "-         load previous configuration file" << endl;
  cout << "=         fit the view" << endl;
  cout << "q         quit" << endl;
  // cout << "" << endl;
  cout << endl;
}

void printInfo() {
  using namespace std;

  cout << "Reference Conf = " << refConfNum << "\n";
  cout << "Current Conf = " << confNum << "\n";
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

    case '0': {
      color_option = 0;
    } break;

    case '1': {
      colorTable.setMinMax(0.5, 1.0);
      colorTable.setTableID(2);
      colorTable.Rebuild();
      color_option = 1;
    } break;

    case '2': {
      colorTable.setMinMax(0.5, 1.0);
      colorTable.setTableID(2);
      colorTable.Rebuild();
      color_option = 2;
    } break;

    case 'a': {
      alpha_particles = std::max(0.0f, alpha_particles - 0.05f);
    } break;
    case 'A': {
      alpha_particles = std::min(1.0f, alpha_particles + 0.05f);
    } break;

    case 'b': {
      alpha_ghosts = std::max(0.0f, alpha_ghosts - 0.05f);
    } break;
    case 'B': {
      alpha_ghosts = std::min(1.0f, alpha_ghosts + 0.05f);
    } break;

    case 'd': {
      show_displacements = 1 - show_displacements;
      if (show_displacements == 1 && show_fluctuations == 1) {
        show_fluctuations = 0;
      }
    } break;

    case 'f': {
      show_fluctuations = 1 - show_fluctuations;
      if (show_fluctuations == 1 && show_displacements == 1) {
        show_displacements = 0;
      }
    } break;

    case 'g': {
      show_ghosts = 1 - show_ghosts;
    } break;

    case 'i': {
      printInfo();
    } break;

    case 'n': {
      std::cout << "Go to file number ";
      int conNumTry;
      std::cin >> conNumTry;
      try_to_readConf(conNumTry, Conf, confNum);
    } break;

    case 'o': {
      showOrientations = 1 - showOrientations;
    } break;

    case 'q': {
      exit(0);
    } break;

    case 'r': {
      ref_fixed = 1 - ref_fixed;
      updateTextLine();
    } break;

    case 'S': {
      vScale *= 1.05;
    } break;

    case 's': {
      vScale *= 0.95;
      if (vScale < 0.0) vScale = 1.0;
    } break;

    case 'w': {
      ghost_width = std::max(0.0, ghost_width - 0.05);
    } break;
    case 'W': {
      ghost_width = std::min(0.5, ghost_width + 0.05);
    } break;

    case '-': {
      int deltaNum = confNum - refConfNum;

      if (ref_fixed == 0 && refConfNum > 0) {
        std::cout << "Reference Configuration: ";
        try_to_readConf(refConfNum - 1, RefConf, refConfNum);
      }
      if (confNum - refConfNum > deltaNum || (ref_fixed == 1 && refConfNum == 0)) {
        std::cout << "Current Configuration: ";
        try_to_readConf(confNum - 1, Conf, confNum);
      }

      updateTextLine();
    } break;

    case '+': {
      std::cout << "Current Configuration: ";
      bool fileRead = try_to_readConf(confNum + 1, Conf, confNum);
      if (fileRead == true && ref_fixed == 0) {
        std::cout << "Reference Configuration: ";
        try_to_readConf(refConfNum + 1, RefConf, refConfNum);
      }

      updateTextLine();
    } break;

    case '=': {
      fit_view();
      reshape(width, height);
    } break;
  };

  glutPostRedisplay();
}

void updateTextLine() {
  textZone.addLine("%s[%d](%0.4g s) -> [%d](%0.4g s), delta = [%d](%0.4g s)", (ref_fixed == 1) ? "*" : " ", refConfNum,
                   RefConf.t, confNum, Conf.t, confNum - refConfNum, Conf.t - RefConf.t);
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
          mouse_mode = PAN;
        else
          mouse_mode = ROTATION;
        break;
      case GLUT_MIDDLE_BUTTON:
        mouse_mode = ZOOM;
        break;
    }
  }
}

void motion(int x, int y) {

  if (mouse_mode == NOTHING) return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

    case ZOOM: {
      double ddy = (worldBox.max.y - worldBox.min.y) * dy;
      double ddx = (worldBox.max.x - worldBox.min.x) * dy;
      worldBox.min.x -= ddx;
      worldBox.max.x += ddx;
      worldBox.min.y -= ddy;
      worldBox.max.y += ddy;
    } break;

    case PAN: {
      double ddx = (worldBox.max.x - worldBox.min.x) * dx;
      double ddy = (worldBox.max.y - worldBox.min.y) * dy;
      worldBox.min.x -= ddx;
      worldBox.max.x -= ddx;
      worldBox.min.y += ddy;
      worldBox.max.y += ddy;
    } break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  display();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  drawBox();
  drawParticles();
  if (show_ghosts == 1) drawGhosts();
  if (show_displacements == 1) drawDisplacements();
  if (show_fluctuations == 1) drawFluctuations();

  textZone.draw();

  glFlush();
  glutSwapBuffers();
}

void fit_view() {
  //
  // 3 x ------- x 2
  //   |         |
  // 0 x ------- x 1
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = Conf.Cell.h.xy;
  double y1 = Conf.Cell.h.yy;
  double x2 = Conf.Cell.h.xy + Conf.Cell.h.xx;
  double y2 = Conf.Cell.h.yy + Conf.Cell.h.yx;
  double x3 = Conf.Cell.h.xx;
  double y3 = Conf.Cell.h.yx;
  worldBox.min.x = std::min(std::min(std::min(x0, x1), x2), x3);
  worldBox.min.y = std::min(std::min(std::min(y0, y1), y2), y3);
  worldBox.max.x = std::max(std::max(std::max(x0, x1), x2), x3);
  worldBox.max.y = std::max(std::max(std::max(y0, y1), y2), y3);
  reshape(width, height);
}

void reshape(int w, int h) {
  width = w;
  height = h;

  double left = worldBox.min.x;
  double right = worldBox.max.x;
  double bottom = worldBox.min.y;
  double top = worldBox.max.y;
  double worldW = right - left;
  double worldH = top - bottom;
  double dW = 0.1 * worldW;
  double dH = 0.1 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;

  if (worldW > worldH) {
    worldH = worldW * ((GLfloat)height / (GLfloat)width);
    top = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * ((GLfloat)width / (GLfloat)height);
    right = 0.5 * (left + right + worldW);
    left = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, bottom, top);

  glutPostRedisplay();
}

// Draw the periodic cell
void drawBox() {
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
  glLineWidth(1.0f);

  glBegin(GL_LINE_LOOP);
  glVertex2f(0.0f, 0.0f);
  glVertex2f(Conf.Cell.h.xy, Conf.Cell.h.yy);
  glVertex2f(Conf.Cell.h.xy + Conf.Cell.h.xx, Conf.Cell.h.yy + Conf.Cell.h.yx);
  glVertex2f(Conf.Cell.h.xx, Conf.Cell.h.yx);
  glEnd();
}

void setColor(int /*i*/, GLfloat alpha) {
  switch (color_option) {

    case 0: {
      glColor4f(0.8f, 0.8f, 0.9f, alpha);
    } break;

    case 1: {
      /*
      colorRGBA col;
      colorTable.getRGB(Conf.grains[i].zncc1, &col);
      glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
      */
    } break;

    case 2: {
      /*
      colorRGBA col;
      colorTable.getRGB(Conf.grains[i].zncc2, &col);
      glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
      */
    } break;

    default: {
      glColor4f(0.8f, 0.8f, 0.9f, alpha);
    } break;
  }
}

void add_ghost_pos(int i, double mn, double mx, std::vector<vec2r>& lst) {
  lst.clear();
  vec2r pos = Conf.Particles[i].pos;
  if (pos.x > mn && pos.x < mx && pos.y > mn && pos.y < mx) {
    return;
  }

  vec2r ghostPos;

  if (pos.x <= mn) {
    ghostPos.set(pos.x + 1.0, pos.y);
    lst.push_back(ghostPos);

    if (pos.y <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.y >= mx) {
      ghostPos.set(pos.x + 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.x >= mx) {
    ghostPos.set(pos.x - 1.0, pos.y);
    lst.push_back(ghostPos);

    if (pos.y <= mn) {
      ghostPos.set(pos.x - 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.y >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.y <= mn) {
    ghostPos.set(pos.x, pos.y + 1.0);
    lst.push_back(ghostPos);

    if (pos.x <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.x >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.y >= mx) {
    ghostPos.set(pos.x, pos.y - 1.0);
    lst.push_back(ghostPos);

    if (pos.x <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.x >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }
}

void drawParticles() {
  glLineWidth(1.0f);

  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos;
    double R = Conf.Particles[i].radius;

    setColor(i, alpha_particles);
    glBegin(GL_POLYGON);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
    }
    glEnd();

    glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);
    glBegin(GL_LINE_LOOP);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
    }
    glEnd();

    if (showOrientations) {
      double rot = Conf.Particles[i].rot;
      glBegin(GL_LINES);
      glVertex2f(pos.x, pos.y);
      glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
      glEnd();
    }
  }
}

void drawGhosts() {
  // if (mouse_mode != NOTHING && box.Particles.size() > 2000) return;

  std::vector<vec2r> lst_pos;  // list of reduced positions of ghost particles
  double mn = ghost_width;
  double mx = 1.0 - ghost_width;
  // GLColorRGBA color;
  glLineWidth(1.0f);
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    add_ghost_pos(i, mn, mx, lst_pos);
    double R = Conf.Particles[i].radius;
    for (size_t ig = 0; ig < lst_pos.size(); ig++) {

      vec2r pos = Conf.Cell.h * lst_pos[ig];

      setColor(i, alpha_ghosts);
      glBegin(GL_POLYGON);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();

      glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();

      if (showOrientations) {
        double rot = Conf.Particles[i].rot;
        glBegin(GL_LINES);
        glVertex2f(pos.x, pos.y);
        glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
        glEnd();
      }
    }
  }
}

void drawContacts() {}

void drawForces() {}

void drawDisplacements() {
  glLineWidth(1.5f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  glBegin(GL_LINES);
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r displ = pos - RefConf.Cell.h * RefConf.Particles[i].pos;

    glVertex2f(pos.x - displ.x, pos.y - displ.y);
    glVertex2f(pos.x, pos.y);
  }
  glEnd();
}

void drawFluctuations() {
  glLineWidth(1.5f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  vec2r vecFluctMean;
  std::vector<vec2r> allFluct(Conf.Particles.size());
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r delta = Conf.Particles[i].pos - RefConf.Particles[i].pos;
    delta.x -= floor(delta.x + 0.5);
    delta.y -= floor(delta.y + 0.5);
    allFluct[i] = Conf.Cell.h * delta;
    vecFluctMean.x += fabs(allFluct[i].x);
    vecFluctMean.y += fabs(allFluct[i].y);
  }
  vecFluctMean /= (double)(Conf.Particles.size());
  double normalisationFactor = 1.0 / norm(vecFluctMean);

  glBegin(GL_LINES);
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r fluct = allFluct[i] * normalisationFactor;
    fluct *= vScale;

    glVertex2f(pos.x - fluct.x, pos.y - fluct.y);
    glVertex2f(pos.x, pos.y);
  }
  glEnd();
}

bool try_to_readConf(int num, PBC& CF, int& OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << file_name << std::endl;
    OKNum = num;
    CF.loadConf(file_name);
    return true;
  } else {
    std::cout << file_name << " does not exist" << std::endl;
  }
  return false;
}

void menu(int num) {
  switch (num) {

    case 0:
      exit(0);
      break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  int submenu1 = glutCreateMenu(menu);  // Particle Colors
  glutAddMenuEntry("None", 100);
  glutAddMenuEntry("Velocity Magnitude", 101);
  glutAddMenuEntry("Sum of Normal Contact Forces", 102);

  int submenu2 = glutCreateMenu(menu);  // Force Colors
  glutAddMenuEntry("None", 200);
  glutAddMenuEntry("Magnitude", 201);

  glutCreateMenu(menu);  // Main menu
  glutAddSubMenu("Particle Colors", submenu1);
  glutAddSubMenu("Force Colors", submenu2);
  glutAddSubMenu("Velocity Colors", submenu2);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {

  if (argc == 1) {
    refConfNum = confNum = 0;
    ref_fixed = 1;
  } else if (argc == 2) {
    refConfNum = confNum = atoi(argv[1]);
    ref_fixed = 1;
  } else if (argc == 3) {
    refConfNum = atoi(argv[1]);
    confNum = atoi(argv[2]);
    ref_fixed = 0;
  }

  std::cout << "Reference Configuration: ";
  try_to_readConf(refConfNum, RefConf, refConfNum);
  std::cout << "Current Configuration: ";
  try_to_readConf(confNum, Conf, confNum);

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("CONF VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  mouse_mode = NOTHING;

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Other initialisations
  glText::init();
  // textZone.addLine("%s[%d] -> [%d]", (ref_fixed == 1) ? "*" : " ", refConfNum, confNum);
  updateTextLine();

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
