#include<iostream>
#include<cmath>
#include<windows.h>
#include<GL/glew.h>
#include<GL/glut.h>
#include "graph.h"
#include "graphicobject.h"

void DrawCircle(float origin_x, float origin_y, float radius);
void Render();

const float pi = 3.141592f;
GraphicObject* gl_input;

int main(int argc, char **argv)
{
    // initiate freeglut and create window

	// glewInit();

    // if (GLEW_OK != glewInit())
    // {
    //     // GLEW failed!
    //     exit(1);
    // }
    Graph testgraph;
    testgraph.GenerateLFRGraph(100, 5, 50, 2, 5, 2, 0.8);
    GraphicObject g(testgraph);
    gl_input = &g;
    // assign global GraphicObject

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(320, 320);
    glutCreateWindow("GraphMaker");

    // Enable alpha blending
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );

    // Init();

    
    //testgraph.adjacency_matrix->Print();

    // register callbacks

    glutDisplayFunc(Render);

    // enter the processing loop

    glutMainLoop();

    return 0;
}

void Render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPointSize(5);
    glLineWidth(2);

    float*& VAO = gl_input->vertex_array;

    // render vertices as circles
    
    for (unsigned int i = 0; i < gl_input->n_vertices; ++i)
    {   
        unsigned int vertex_ptr = i * 7;
        glColor4f(VAO[vertex_ptr + 3], VAO[vertex_ptr + 4], VAO[vertex_ptr + 5], VAO[vertex_ptr + 6]);
    //    glVertex3f(VAO[vertex_ptr], VAO[vertex_ptr + 1], VAO[vertex_ptr + 2]);
        DrawCircle(VAO[vertex_ptr], VAO[vertex_ptr + 1], 0.02);
    }

    // render edges as lines

    glBegin(GL_LINES);
    for (unsigned int i = 0; i < gl_input->n_indices; ++i)
    {   
        unsigned int vertex_ptr = gl_input->index_array[i] * 7;
        glColor4f(VAO[vertex_ptr + 3], VAO[vertex_ptr + 4], VAO[vertex_ptr + 5], 0.3f);
        glVertex3f(VAO[vertex_ptr], VAO[vertex_ptr + 1], VAO[vertex_ptr + 2]);
    }
    glEnd();

	glutSwapBuffers();
}

void DrawCircle(float origin_x, float origin_y, float radius)
{
    glBegin(GL_LINE_LOOP);
    for (float i = 0; i < 2 * pi; i += 1 / 8.0f)
    {
        glVertex2f(radius * std::cos(i) + origin_x, radius * std::sin(i) + origin_y);
    }
    glEnd();
}


