#include <GL/glut.h>
#include <GL/gl.h>
#include <unistd.h>
#include <stdio.h>
#include <GL/glu.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define CIRCLES 2
#define DIMS 2

const int step = 15;

double dt = 0.015, g = 9.80665, dtnew = 1, r = 0.5;
int n = 0, size = 10;
double colTime;
int cir, dir;

double pos[CIRCLES][DIMS], veloc[CIRCLES][DIMS];

void nextHit(void);
void wallBounce(void);
void drawShapes(void);
void drawCircle(int);
void update(const int);
void manageResize(int, int);
double skalsuc(const double[CIRCLES], const double[CIRCLES]);

int main(int argc, char *argv[]){
    srand(time(0));

    for (int i = 0; i < CIRCLES; i++){
        int alfa = rand() % 360;
        int v0 = rand() % size + 5;

        pos[i][1] = (i ? 1 : -1) * (size / 4) + (i ? -r - (r / 9) : r + (r / 9));
        pos[i][0] = (i ? 1 : -1) * (size / 3) + (i ? -r - (r / 9) : r + (r / 9));

        veloc[i][0] = (v0 * cos(alfa * ((M_PI / 180))));
        veloc[i][1] = (v0 * sin(alfa * ((M_PI / 180))));
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE);
    glutInitWindowSize(1000, 600);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Jan Agh - biliard");
    glutDisplayFunc(drawShapes);
    glClearColor(0.2, 0.8, 0.8, 0);
    glutReshapeFunc(manageResize);
    glutTimerFunc(step, update, 1);
    glutMainLoop();
    return 0;
}

double skalsuc(const double first[CIRCLES], const double second[CIRCLES]){
    double sum = 0.0;
    for (int cart = 0; cart < CIRCLES; cart++){
        sum += first[cart] * second[cart];
    }
    return sum;
}

void update(const int ihod){
    if(ihod == 1){
        nextHit();
        n = 1 + (int) (colTime / dt);
        dtnew = (colTime / n) < 1e-6 ? 1e-6 : colTime / n;
    }

    for (int i = 0; i < CIRCLES; i++){
        for (int j = 0; j < DIMS; j++){
            pos[i][j] += veloc[i][j] * dtnew;
        }
    }
    glutPostRedisplay();
        
    if(ihod == n){
        double rx = pos[1][0] - pos[0][0];
        double ry = pos[1][1] - pos[0][1];
        double vx = veloc[1][0] - veloc[0][0];
        double vy = veloc[1][1] - veloc[0][1];
        double distance = sqrt(rx * rx + ry * ry);

        double Ix = rx * (rx * vx + ry * vy) / (distance * distance);
        double Iy = ry * (rx * vx + ry * vy) / (distance * distance);
        
        if(dir > -1){
            wallBounce();
        }
        else if(dir == -1){
            if (distance <= 2 * r){
                veloc[0][0] += Ix;
                veloc[0][1] += Iy;
                veloc[1][0] += -Ix;
                veloc[1][1] += -Iy;
            }
        }
        glutTimerFunc(step, update, 1);
    }
    else{
        glutTimerFunc(step, update, ihod + 1);
    }
}

void wallBounce(void){
    if ((pos[cir][0] + r > size / 2 - (size / 10)) || (pos[cir][0] - r < -size / 2 + (size / 10))){
        veloc[cir][0] *= -1;
    }
    if ((pos[cir][1] + r > size / 3 - (size / 10)) || (pos[cir][1] - r < -(size / 3 - (size / 10)))){
        veloc[cir][1] *= -1;
    }
}

void drawShapes(void){
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(0, 1.0, 0.5);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glRectf(-size / 2 + (size / 10), -size / 3 + (size / 10), size / 2 - (size / 10), size / 3 - (size / 10));

    for (int i = 0; i < CIRCLES; i++){
        glLoadIdentity();
        glTranslatef(pos[i][0], pos[i][1], 0.0);
        glScalef(r, r, 1.0);
        drawCircle(i);
    }
    glutSwapBuffers();
}

void drawCircle(int id)
{
    int pocLucov = 360;
    int x = 0, y = 0;

    glBegin(GL_TRIANGLE_FAN); 
    glColor3f(id == 0 ? 1 : 0, 0, id == 0 ? 0 : 1);
    glVertex2f(x, y);
    
    for (int i = 0; i <= pocLucov; i++){ 
        double fin = (i * 2 * M_PI) / pocLucov;
        glVertex2f(x + sin(fin), y + cos(fin));
    }
    glEnd();
}

void nextHit(void){
    double strany[2] = {size / 2 - (size / 10.0), size / 3 - (size / 10.0)};
    double min = DBL_MAX;
    double minCandidate = 0;
    
    for (int i = 0; i < CIRCLES; i++){
        for (int j = 0; j < DIMS; j++){
            if (veloc[i][j] < 0){
                minCandidate = (pos[i][j] - (-1 * strany[j] + r)) / (-veloc[i][j]);
                if (minCandidate < min){
                    min = minCandidate, cir = i, dir = j;
                }
            }
            else if (veloc[i][j] > 0){
                minCandidate = ((strany[j] - (pos[i][j] + r))) / veloc[i][j];
                if (minCandidate < min){
                    min = minCandidate, cir = i, dir = j;
                }
            }
            else{
                min = 1e-6;
            }
        }
    }
    
    double bij = 0.0, rsq = 0.0, vsq = 0.0;
    for (int i = 0; i < DIMS; i++){
        bij += (pos[1][i] - pos[0][i]) * (veloc[1][i] - veloc[0][i]);
        rsq += (pos[1][i] - pos[0][i]) * (pos[1][i] - pos[0][i]);
        vsq += (veloc[1][i] - veloc[0][i]) * (veloc[1][i] - veloc[0][i]);
    }
    if (bij < 0.0 && bij * bij - vsq * (rsq - 4 * r * r) > 0.0){
        double time = (-bij - sqrt(bij * bij - (rsq - 4 * r * r) * vsq)) / vsq;
        if (time >= 0.0 && time < min){
            min = time, dir = -1;
        }
    }
    colTime = min;
}

void manageResize(int width, int height){
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (width == 0) width++;
    const float sideRatio = ((float)height / width);

    gluOrtho2D(-0.5 * size, 0.5 * size, -0.5 * size * sideRatio, 0.5 * size * sideRatio);
}