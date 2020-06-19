#pragma once

using namespace Eigen;

static const int n = 20;

static const double dx = 1.0 / (n - 1);
static const int X_GRID = int(1 / dx);
static const int Y_GRID = int(1 / dx);
static const int X_WINDOW = 1400;
static const int Y_WINDOW = X_WINDOW * Y_GRID / X_GRID;


#define PI 3.1415927410125732421875