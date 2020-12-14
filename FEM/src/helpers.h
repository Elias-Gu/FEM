#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "constants.h"


struct Output
{
	Output() {}
	~Output() {}


	void OutputNodeCoordinates(const std::vector<Eigen::Vector2d>& nodes_coo, const Eigen::VectorXd& sol)
	{
		std::ofstream output;
		output.open("out/nodes_coo.txt");
		for (size_t i = 0; i < nodes_coo.size(); i++)
		{
			std::string str = std::to_string(nodes_coo[i][0]) + " , " +
				std::to_string(nodes_coo[i][1]) + " , " +
				std::to_string(sol[i]);
			output << str << "\n";
		}
		output.close();
	}


	void OutputNodeCoordinates(const std::vector<Eigen::Vector2d>& nodes_coo, const Eigen::VectorXd& sol, const int step)
	{
		std::ofstream output;
		output.open("out/steps/nodes_coo_" + std::to_string(step) + ".txt");
		for (size_t i = 0; i < nodes_coo.size(); i++)
		{
			std::string str = std::to_string(nodes_coo[i][0]) + " , " +
				std::to_string(nodes_coo[i][1]) + " , " +
				std::to_string(sol[i]);
			output << str << "\n";
		}
		output.close();
	}


	void OutputMesh(const std::vector<Eigen::Vector3i>& mesh)
	{
		std::ofstream output;
		output.open("out/mesh.txt");
		for (size_t i = 0; i < mesh.size(); i++)
		{
			std::string str = std::to_string(mesh[i][0]) + " , " +
				std::to_string(mesh[i][1]) + " , " +
				std::to_string(mesh[i][2]);
			output << str << "\n";
		}
		output.close();
	}
};


struct ScopedTimer
{
	std::string task_name;
	double* time;
	bool verbose;
	std::chrono::high_resolution_clock::time_point start;
	

	ScopedTimer(const std::string& _task_name = "", double* _time = nullptr, const bool _verbose = true)
		: task_name(_task_name),
		time(_time),
		verbose(_verbose),
		start(std::chrono::high_resolution_clock::now())
	{}


	~ScopedTimer()
	{
		auto stop = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = stop - start;
		double time_elapsed = (double)(elapsed.count() * 1000.0);

		if (time)
			*time = time_elapsed;

		if (verbose)
			std::cout << task_name << " ran in " << time_elapsed << " ms." << std::endl;
	}
};