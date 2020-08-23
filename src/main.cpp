#include<iostream>
#include<cmath>
#include<windows.h>
#include<GL/glew.h>
#include<GL/glut.h>
#include "graph.h"
#include "graphicobject.h"

void DrawCircle(float origin_x, float origin_y, float radius, bool filled);
void Render();

const float pi = 3.141592f;
GraphicObject* gl_input;

int main(int argc, char **argv)
{
	// glewInit();

    // if (GLEW_OK != glewInit())
    // {
    //     // GLEW failed!
    //     std::cout << "fail!";
    //     exit(1);
    // }

    // assign global GraphicObject
    Graph testgraph;
    testgraph.GenerateLFRGraph(100, 4, 40, 2.5, 7, 2.0, 0.2);
    testgraph.Louvain(testgraph, 0);
    //testgraph.GetMixingParameter(true);
    //testgraph.GetMixingParameter(false);
    GraphicObject g(testgraph);
    gl_input = &g;

    // initialise freeglut and open window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutCreateWindow("GraphMaker");

    // enable alpha blending
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );

    // register callbacks
    glutDisplayFunc(Render);

    // enter the processing loop
    glutMainLoop();

    return 0;
}

// a_max = pi r_max^2
// r_max = sqrt (a_max / pi)

void Render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLineWidth(2);
    float max_vertex_area = 0.06f;
    const unsigned int& max_degree = gl_input->source_graph->max_degree;
    const std::vector<unsigned int>& degree = *gl_input->source_graph->degree;

    float*& VAO = gl_input->vertex_array;

    // render vertices as circles
    
    for (unsigned int i = 0; i < gl_input->n_vertices; ++i)
    {   
        unsigned int vertex_ptr = i * 7;
        float radius = (((float)degree[i] / (float)max_degree) * max_vertex_area) + 0.01f;
        glColor4f(VAO[vertex_ptr + 3], VAO[vertex_ptr + 4], VAO[vertex_ptr + 5], 0.6f);
        DrawCircle(VAO[vertex_ptr], VAO[vertex_ptr + 1], radius, true);
        glColor4f(VAO[vertex_ptr + 3], VAO[vertex_ptr + 4], VAO[vertex_ptr + 5], 1.0f);
        DrawCircle(VAO[vertex_ptr], VAO[vertex_ptr + 1], radius, false);
    }

    // render edges as lines
    glLineWidth(1);
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

void DrawCircle(float origin_x, float origin_y, float radius, bool filled)
{
    glBegin(filled ? GL_TRIANGLE_FAN : GL_LINE_LOOP);
    for (float i = 0; i < 2 * pi; i += 1 / 20.0f)
    {
        glVertex2f(radius * std::cos(i) + origin_x, radius * std::sin(i) + origin_y);
    }
    glEnd();
}


