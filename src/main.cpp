#include <iostream>

#include <GLFW/glfw3.h>


/* Declarations */
void initGLContext();
GLFWwindow* initGLFWContext();
// TODO: Clear this
const double dx = 0.001;
const int X_GRID = int(1 / dx);
const int Y_GRID = int(1 / dx);
const int X_WINDOW = 1400;
const int Y_WINDOW = X_WINDOW * Y_GRID / X_GRID;


/* -----------------------------------------------------------------------
|								MAIN									 |
----------------------------------------------------------------------- */


int main()
{
	GLFWwindow* window = initGLFWContext();
	initGLContext();
	while (!glfwWindowShouldClose(window))
	{
		glClear(GL_COLOR_BUFFER_BIT);

		// TODO: Draw call here
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

	glClearColor(.55f, .55f, .55f, .0f);
	glClear(GL_COLOR_BUFFER_BIT);
}