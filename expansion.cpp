#include <GL/glut.h>
#include <cstdlib>
#include <cmath>
#include <random>
#include <ctime>

#define N 1000

struct Vector {
    double x;
    double y;
    double z;
};

struct Particle {
    Vector p;
    Vector q;
    Vector F;
    double d;
};

std::mt19937 ran((int)time(0));
std::uniform_real_distribution<double> distr(0,1);

double Random() {                                //Random number in [0,1].
    return distr(ran);
}

Particle* p = new Particle [N+1];
int* order = new int [N+1];
Vector zero;
double R = 1e4;
double T = 0.001;
double L;
int period = 1;

inline double dotProduct(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double amplitude(Vector q) {
    return sqrt(dotProduct(q, q));
}

inline Vector addition(Vector a, Vector b) {
    return Vector {
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    };
}

inline Vector multiplication(double a, Vector v) {
    return Vector {
        a * v.x,
        a * v.y,
        a * v.z
    };
}

inline Vector direction(Vector r) {
    double amp = amplitude(r);
    return Vector {
        r.x / amp,
        r.y / amp,
        r.z / amp
    };
}

inline double distance(Vector q_1, Vector q_2) {
    return amplitude(addition(q_1, multiplication(-1.0, q_2)));
}

Vector eye = {
    0.0,
    -2.2e6,
    1.5e6
};

void Sort(int l, int r) {
    if(r <= l)
        return;
    int i = l;
    int j = r;
    while(i < j) {
        while(p[order[i]].d >= p[order[l]].d && i < r)
            i++;
        while(p[order[j]].d <= p[order[l]].d && j > l)
            j--;
        if(i >= j)
            break;
        int tmp = order[i];
        order[i] = order[j];
        order[j] = tmp;
    }
    int tmp = order[j];
    order[j] = order[l];
    order[l] = tmp;
    Sort(l, j - 1);
    Sort(j + 1, r);
}

void expansionDisplay() {
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0, 2.0, 1e5, 1e7);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye.x, eye.y, eye.z, 0.0, 0.5e6, 0.45e6, 0.0, 0.0, 1.0);
    
    glColor3f(0.5f, 0.5f, 0.5f);
    glBegin(GL_LINES);                  // Container
    glLineWidth(1e4f);
    
    glVertex3d(-1e6, 0.0, 0.0);
    glVertex3d(-1e6, 1e6, 0.0);
    glVertex3d(-1e6, 1e6, 0.0);
    glVertex3d(-1e6, 1e6, 1e6);
    glVertex3d(-1e6, 1e6, 1e6);
    glVertex3d(-1e6, 0.0, 1e6);
    glVertex3d(-1e6, 0.0, 1e6);
    glVertex3d(-1e6, 0.0, 0.0);
    
    glVertex3d(1e6, 0.0, 0.0);
    glVertex3d(1e6, 1e6, 0.0);
    glVertex3d(1e6, 1e6, 0.0);
    glVertex3d(1e6, 1e6, 1e6);
    glVertex3d(1e6, 1e6, 1e6);
    glVertex3d(1e6, 0.0, 1e6);
    glVertex3d(1e6, 0.0, 1e6);
    glVertex3d(1e6, 0.0, 0.0);
    
    glVertex3d(-1e6, 0.0, 0.0);
    glVertex3d(1e6, 0.0, 0.0);
    glVertex3d(-1e6, 1e6, 0.0);
    glVertex3d(1e6, 1e6, 0.0);
    glVertex3d(-1e6, 1e6, 1e6);
    glVertex3d(1e6, 1e6, 1e6);
    
    glEnd();
    
    Sort(1, N);							// Draw paricles in the back first.
    for(int i=1;i<=N;i++) {             // Particles
        double l = 1.0/3.0 + 2.0/3.0 * L * L / p[order[i]].d / p[order[i]].d;
        glColor3f(l, l, 0.0f);			// Light on a particle is 1/(r^2)
        glTranslated(p[order[i]].q.x, p[order[i]].q.y, p[order[i]].q.z);
        glutSolidSphere(R, 10, 10);
        Vector q = multiplication(-1.0, p[order[i]].q);
        glTranslated(q.x, q.y, q.z);
    }
    
    glColor3f(0.5f, 0.5f, 0.5f);
    glBegin(GL_LINES);                  // Container - one last edge
    glLineWidth(1e4f);
    
    glVertex3d(-1e6, 0.0, 1e6);
    glVertex3d(1e6, 0.0, 1e6);
    
    glEnd();
    
    glFlush();
    glutSwapBuffers();
}

void timeEvolution() {
    for(int t=1;t<=period;t++) {
        for(int i=1;i<=N;i++) {             // Moving
            p[i].q.x += p[i].p.x * T;
            if(p[i].q.x > 1e6) {
                p[i].q.x = 2e6 - p[i].q.x;
                p[i].p.x = 0.0 - p[i].p.x;
            } else if(p[i].q.x < -1e6) {
                p[i].q.x = -2e6 - p[i].q.x;
                p[i].p.x = 0.0 - p[i].p.x;
            }
            p[i].q.y += p[i].p.y * T;
            if(p[i].q.y > 1e6) {
                p[i].q.y = 2e6 - p[i].q.y;
                p[i].p.y = 0.0 - p[i].p.y;
            } else if(p[i].q.y < 0.0) {
                p[i].q.y = 0.0 - p[i].q.y;
                p[i].p.y = 0.0 - p[i].p.y;
            }
            p[i].q.z += p[i].p.z * T;
            if(p[i].q.z > 1e6) {
                p[i].q.z = 2e6 - p[i].q.z;
                p[i].p.z = 0.0 - p[i].p.z;
            } else if(p[i].q.z < 0.0) {
                p[i].q.z = 0.0 - p[i].q.z;
                p[i].p.z = 0.0 - p[i].p.z;
            }
            p[i].d = distance(eye, p[i].q);
            order[i] = i;
        }
        for(int i=1;i<N;i++) {				// Collisions
            for(int j=i+1;j<=N;j++) {
                if(distance(p[i].q, p[j].q) <= 2.0 * R) {
                    Vector r = direction(addition(p[j].q, multiplication(-1.0, p[i].q)));
                    double impulse = dotProduct(p[i].p, r) - dotProduct(p[j].p, r);
                    p[i].p = addition(p[i].p, multiplication(-1.0 * impulse, r));
                    p[j].p = addition(p[j].p, multiplication(impulse, r));
                }
            }
        }
    }
    expansionDisplay();
}

bool legal(Vector q, int index) {
    if(q.x < R-1e6 || q.x > 1e6-R ||q.y < R || q.y > 1e6-R || q.z < R || q.z > 1e6-R)
        return false;
    for(int i=1;i<index;i++)
        if(distance(q, p[i].q) < 2.0 * R)
            return false;
    return true;
}

void initiation() {
    for(int i=1;i<=N;i++) {
        do {
            p[i].q.x = Random() * (-1e6);
            p[i].q.y = Random() * 1e6;
            p[i].q.z = Random() * 1e6;
        } while(!legal(p[i].q, i));
        p[i].p.x = Random() * 2e6 - 1e6;
        p[i].p.y = Random() * 2e6 - 1e6;
        p[i].p.z = Random() * 2e6 - 1e6;
        p[i].d = distance(eye, p[i].q);
        order[i] = i;
    }
    L = distance(eye, {0.0, 0.0, 1e6});
}


int main(int argc, char** argv) {
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowPosition(60, 120);
    glutInitWindowSize(1400, 700);
    glutCreateWindow("Gas Expansion (Simulated with Particles)");
    initiation();
    glutDisplayFunc(&expansionDisplay);
    glutIdleFunc(&timeEvolution);
    glutMainLoop();
    
    return 0;
}
