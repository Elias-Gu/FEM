#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <GLFW/glfw3.h>

#include "constants.h"
#include "helpers.h"
#include "solver.h"


/* Declarations */
void initGLContext();
GLFWwindow* initGLFWContext();
Solver* simulation;



/* -----------------------------------------------------------------------
|						  FINITE ELEMENT METHOD							 |
----------------------------------------------------------------------- */


void FiniteElementMethod()
{
	
}



/* -----------------------------------------------------------------------
|								MAIN									 |
----------------------------------------------------------------------- */


int main()
{	
	simulation = new Solver();
	Output output; 

	int step = 0;
	simulation->Init();
	output.OutputNodeCoordinates(simulation->nodes_coo, simulation->val,step);
	while (simulation->tn < 8)//step < 2)
	{
		step++;
		simulation->Step();
		output.OutputNodeCoordinates(simulation->nodes_coo, simulation->val, step);
	}

	std::cin.get();
	return 0;
}

//int main()
//{
//	Renderer renderer;
//	Solver *Simulation = new Solver();
//
//	GLFWwindow* window = initGLFWContext();
//	initGLContext();
//	while (!glfwWindowShouldClose(window))
//	{
//		glClear(GL_COLOR_BUFFER_BIT);
//
//		renderer.DrawMesh(Simulation->nodes_coo, Simulation->mesh);
//
//		glfwSwapBuffers(window);
//		glfwPollEvents();
//	}
//
//	glfwTerminate();
//
//	return 0;
//}



/* -----------------------------------------------------------------------
|						     	  OPENGL 								 |
----------------------------------------------------------------------- */


/* Define window. */
GLFWwindow* initGLFWContext()
{
	if (!glfwInit())
		exit(EXIT_FAILURE);

	/* Create Window*/
	GLFWwindow* window = glfwCreateWindow(X_WINDOW, Y_WINDOW, "Simulation", NULL, NULL);

	if (!window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(window);
	return window;
}


/* Define view parameters */
void initGLContext()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, X_GRID, 0, Y_GRID, -1, 1);					// original
	glViewport(0, 0, (GLsizei)X_WINDOW, (GLsizei)Y_WINDOW);	// transfo

	glClearColor(.7f, .7f, .7f, .0f);
	glClear(GL_COLOR_BUFFER_BIT);
}