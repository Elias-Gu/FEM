#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <GLFW/glfw3.h>

#include "constants.h"
#include "helpers.h"
#include "heat_solver.h"
#include "poisson_solver.h"


/* Declarations */
void initGLContext();
GLFWwindow* initGLFWContext();
HeatSolver* simulation;
PoissonSolver* simulation2;



/* -----------------------------------------------------------------------
|						  FINITE ELEMENT METHOD							 |
----------------------------------------------------------------------- */


void FiniteElementMethod()
{
	
}



/* -----------------------------------------------------------------------
|								MAIN									 |
----------------------------------------------------------------------- */


//int main()
//{	
//	simulation = new HeatSolver();
//	simulation2 = new PoissonSolver();
//	simulation2->FEMSolver();
//	Output output; 
//
//	int step = 0;
//	simulation->Init();
//	output.OutputNodeCoordinates(simulation->nodes_coo, simulation->val,step);
//	while (simulation->tn < 8)//step < 2)
//	{
//		step++;
//		simulation->Step();
//		output.OutputNodeCoordinates(simulation->nodes_coo, simulation->val, step);
//	}
//
//	std::cin.get();
//	return 0;
//}

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

int main()
{
	Renderer renderer;
	int Nn = n * n;
	int Ne = 2 * (n - 1) * (n - 1);

	std::vector<Vector2d> nodes_coo;
	nodes_coo.resize(Nn);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++)
		{
			nodes_coo[i * n + j] = Vector2d(i * dx, j * dx);
			//nodes_coo[i * n + j][0] += 0.25;
			//nodes_coo[i * n + j][1] += 0.25;
		}
	}

	// Nodes connectivity
	std::vector<Vector3i> mesh;
	mesh.resize(Ne);
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++)
		{
			int q = (n - 1) * i + j;
			mesh[2 * q] = Eigen::Vector3i((i + 1) * n + j, (i + 1) * n + j + 1, i * n + j);
			mesh[2 * q + 1] = Eigen::Vector3i(i * n + j + 1, i * n + j, (i + 1) * n + j + 1);
		}
	}

	GLFWwindow* window = initGLFWContext();
	initGLContext();
	while (!glfwWindowShouldClose(window))
	{
		glClear(GL_COLOR_BUFFER_BIT);
	
		renderer.DrawMesh(nodes_coo, mesh);
	
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	
	glfwTerminate();
	
	return 0;

}



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