#include "renderer.h"

void Renderer::DrawConnectivity(const std::vector<Eigen::Vector2d>& nodes_coo, const std::vector<Eigen::Vector3i>& connectivity)
{
	glLineWidth(4);
	glColor3f(0.6f, 0.6f, 0.6f);

	glBegin(GL_LINES);
	for (size_t i = 0; i < connectivity.size(); i++)
	{
		glVertex2f(float(nodes_coo[connectivity[i][0]][0] * X_GRID), float(nodes_coo[connectivity[i][0]][1] * Y_GRID));
		glVertex2f(float(nodes_coo[connectivity[i][1]][0] * X_GRID), float(nodes_coo[connectivity[i][1]][1] * Y_GRID));
		glVertex2f(float(nodes_coo[connectivity[i][2]][0] * X_GRID), float(nodes_coo[connectivity[i][2]][1] * Y_GRID));
		glVertex2f(float(nodes_coo[connectivity[i][0]][0] * X_GRID), float(nodes_coo[connectivity[i][0]][1] * Y_GRID));
	}
	glEnd();
}


void Renderer::DrawPoints(const std::vector<Eigen::Vector2d> &nodes_coo)
{
	glPointSize(5.0f);
	glColor3f(0.3f, 0.3f, 0.3f);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);
	for (size_t i = 0; i < nodes_coo.size(); i++)
	{
		glVertex2f(float(nodes_coo[i][0] * X_GRID), float(nodes_coo[i][1] * Y_GRID));
	}
	glEnd();
}


void Renderer::DrawMesh(const std::vector<Eigen::Vector2d>& nodes_coo, const std::vector<Eigen::Vector3i>& connectivity)
{
	DrawConnectivity(nodes_coo, connectivity);
	DrawPoints(nodes_coo);
}