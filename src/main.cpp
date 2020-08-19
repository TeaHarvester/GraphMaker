#include<iostream>
#include<GL/glut.h>
#include "graph.h"
#include "graphicobject.h"

void render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_TRIANGLES);
    glVertex3f(-0.5,-0.5,0.0);
	glVertex3f(0.5,0.0,0.0);
	glVertex3f(0.0,0.5,0.0);
	glEnd();

	glutSwapBuffers();
}

int main(int argc, char **argv)
{
    // initiate freeglut and create window

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(320, 320);
    glutCreateWindow("GraphMaker");

    // register callbacks

    glutDisplayFunc(render);

    // enter the processing loop

    //glutMainLoop();

    Graph testgraph; 
    testgraph.GenerateLFRGraph(20, 2, 10, 2, 3, 2, 0.8);
    testgraph.adjacency_matrix->Print();

    GraphicObject renderate(testgraph);

    return 0;
}

