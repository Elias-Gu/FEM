#pragma once

#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <GLFW/glfw3.h>

#include "constants.h"


struct Renderer
{
	Renderer() {}
	~Renderer() {}
	
	void DrawMesh(const std::vector<Eigen::Vector2d>& nodes_coo, const std::vector<Eigen::Vector3i>& connectivity);
	void DrawPoints(const std::vector<Eigen::Vector2d>& nodes_coo);
	void DrawConnectivity(const std::vector<Eigen::Vector2d>& nodes_coo, const std::vector<Eigen::Vector3i>& connectivity);

	void OutputNodeCoordinates(const std::vector<Eigen::Vector2d>& nodes_coo, const Eigen::VectorXd& sol);
	void OutputNodeCoordinates(const std::vector<Eigen::Vector2d>& nodes_coo, const Eigen::VectorXd& sol, const int step);
	void OutputMesh(const std::vector<Eigen::Vector3i>& mesh);
};