#include "loads.h"

/* -----------------------------------------------------------------------
|							 INTERNAL FORCES							 |
----------------------------------------------------------------------- */


double Loads::InternalForce(const Vector2d& coo)
{
	double force_f = 4 * double(PI) * sin(2 * double(PI) * coo[0]) * cos(2 * double(PI) * coo[1]);

	return force_f;
}


double Loads::InternalForce(const Vector2d& coo, const double tn)
{
	double force_f = 4 * double(PI) * sin(double(PI) * coo[0]) * sin(double(PI) * coo[1]);

	return force_f;
}



/* -----------------------------------------------------------------------
|							NEUMANN FORCES								 |
----------------------------------------------------------------------- */


std::vector<Vector2i> Loads::NeumannEdges(const std::vector<Vector2d>& coo)
{
	std::vector<Vector2i> neumann_edges;
	for (int k = 0; k < n - 1; k++)									// Top edge
		neumann_edges.push_back(Vector2i((k + 1) * n - 1, (k + 2) * n - 1));
	for (int k = 0; k < n - 1; k++)									// Right edge
		neumann_edges.push_back(Vector2i(n * n - (k + 1), n * n - (k + 2)));
	for (int k = 0; k < n - 1; k++)									// Bottom edge
		neumann_edges.push_back(Vector2i((n - 1 - k) * n, (n - 2 - k) * n));

	return neumann_edges;
}


double Loads::NeumannForce(const Vector2d& coo, const Vector2d& normal)
{
	Vector2d force_g = Vector2d(double(PI)/5.0 * cos(2 * double(PI) * coo[0]), sin(2 * double(PI) * coo[1]));

	return force_g.dot(normal);
}


double Loads::NeumannForce(const Vector2d& coo, const Vector2d& normal, const double tn)
{
	Vector2d force_g;

	if (tn < 2)
		force_g = Vector2d(2 * double(PI) * sin(2 * double(PI) * coo[1]), 2 * double(PI) * sin(2 * double(PI) * coo[1]));
	else
		force_g = Vector2d(0, 0);

	if (coo[0] == 1)
		return 4.0 * force_g.dot(normal);
	else
		return 0.0;
}



/* -----------------------------------------------------------------------
|							DIRICHLET FORCES							 |
----------------------------------------------------------------------- */


std::vector<bool> Loads::DirichletNodes(const std::vector<Vector2d>& coo)
{
	std::vector<bool> dirichlet_nodes(coo.size(), false);

	for (size_t i = 0; i < coo.size(); i++)
		if (i < n)
			dirichlet_nodes[i] = true;

	return dirichlet_nodes;
}


double Loads::DirichletValue(const Vector2d& coo)
{
	double value_d = cos(2 * double(PI) * coo[0]) * cos(2 * double(PI) * coo[1]);

	return value_d; 
}


double Loads::DirichletValue(const Vector2d& coo, const double tn)
{
	return 4.0 * cos(2 * double(PI) * tn);
}


double Loads::DtDirichletValue(const Vector2d& coo, const double tn)
{
	return -8.0 * double(PI) * sin(2 * double(PI) * tn);
}