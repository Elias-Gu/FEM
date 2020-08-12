#pragma once

using namespace Eigen;

static const int n = 40;
static const double dx = 1.0 / (n - 1);
static const double dt = 0.01;

static const int X_GRID = int(1 / dx);
static const int Y_GRID = int(1 / dx);
static const int X_WINDOW = 1400;
static const int Y_WINDOW = X_WINDOW * Y_GRID / X_GRID;


constexpr auto PI = 3.1415927410125732421875;