#include<iostream>
#include<windows.h>
#include<GL/glew.h>
#include<GL/glut.h>
#include "graph.h"
#include "graphicobject.h"
// #include "renderer.h"


// static unsigned int shaderProgram;

// int Init();
void Render();

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
    testgraph.GenerateLFRGraph(40, 5, 10, 2, 3, 2, 0.8);
    GraphicObject g(testgraph);
    gl_input = &g;
    // assign global GraphicObject

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(320, 320);
    glutCreateWindow("GraphMaker");

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

    float*& VAO = gl_input->vertex_array;

    glBegin(GL_POINTS);
    for (unsigned int i = 0; i < gl_input->n_vertices; ++i)
    {   
        unsigned int vertex = i * 7;
        glColor4f(VAO[vertex + 3], VAO[vertex + 4], VAO[vertex + 5], VAO[vertex + 6]);
        glVertex3f(VAO[vertex], VAO[vertex + 1], VAO[vertex + 2]);
    }

    // glVertex3f(-0.5,-0.5,0.0);
	// glVertex3f(0.5,0.0,0.0);
	// glVertex3f(0.0,0.5,0.0);

	glEnd();

	glutSwapBuffers();
}

// int Init()
// {
// 	const char* vertexShaderSource;
// 	const char* fragmentShaderSource;

//     vertexShaderSource = "#version 330 core\n"
// 		"layout (location = 0) in vec3 aPos;\n"
// 		"layout (location = 1) in vec4 aCol;\n"
// 		"out vec4 colour;\n"
// 		"uniform mat4 model;\n"
// 		"uniform mat4 view;\n"
// 		"uniform mat4 projection;\n"
// 		"void main()\n"
// 		"{\n"
// 		"	gl_Position = projection * view * model * vec4(aPos, 1.0);\n"
// 		"	colour = aCol;\n"
// 		"}\0";

// 	fragmentShaderSource = "#version 330 core\n"
// 		"in vec4 colour;\n"
// 		"out vec4 FragColour;\n"
// 		"void main()\n"
// 		"{\n"
// 		"	FragColour = colour;\n"
// 		"}\n\0";

//     glEnable(GL_DEPTH_TEST);

//     unsigned int vertexShader;
//     vertexShader = glCreateShader(GL_VERTEX_SHADER);
//     glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
// 	glCompileShader(vertexShader);

// 	unsigned int fragmentShader;
// 	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
// 	glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
// 	glCompileShader(fragmentShader);

// 	int success;
// 	char infoLog[512];
// 	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);

// 	if (!success)
// 	{
// 		glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
// 		std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
// 	}

// 	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);

// 	if (!success)
// 	{
// 		glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
// 		std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
// 	}

// 	shaderProgram = glCreateProgram();
// 	glAttachShader(shaderProgram, vertexShader);
// 	glAttachShader(shaderProgram, fragmentShader);
// 	glLinkProgram(shaderProgram);

// 	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
// 	if (!success) {
// 		glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
// 		std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
// 	}

// 	glDeleteShader(vertexShader);
// 	glDeleteShader(fragmentShader);

// 	return 0;
// }

// void Render()
// {
//     GLuint VAO;
// 	GLuint VBO; 
// 	GLuint EBO;

// 	// create buffer and vertex array object
// 	glGenVertexArrays(1, &VAO);
// 	glGenBuffers(1, &VBO);
// 	glGenBuffers(1, &EBO);

// 	// Initialisation for objects that infrequently change

// 	glBindVertexArray(VAO);
// 	glBindBuffer(GL_ARRAY_BUFFER, VBO);
// 	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * g->n_vertices * 7, g->vertices, GL_STATIC_DRAW);

// 	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
// 	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * g->n_indices * 2, g->indices, GL_STATIC_DRAW);

// 	// set vertex attribute pointers

// 	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(0));
// 	glEnableVertexAttribArray(0);

// 	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(3 * sizeof(float)));
// 	glEnableVertexAttribArray(1);

// 	glBindBuffer(GL_ARRAY_BUFFER, 0);
// 	glBindVertexArray(0);

//     // input
//     // rendering commands
//     glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
//     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//     glUseProgram(shaderProgram);
//     glBindVertexArray(VAO);
//     glPointSize(10.0f);

//     glDrawArrays(GL_POINTS, 0, g->n_vertices);
//     glDrawElements(GL_LINES, g->n_indices, GL_UNSIGNED_INT, 0);

//     glBindVertexArray(0);

//     // check and call events, swap the buffers
//     glutSwapBuffers();
//     glDeleteVertexArrays(1, &VAO);
//     glDeleteBuffers(1, &VBO);
//     glDeleteBuffers(1, &EBO);
// }