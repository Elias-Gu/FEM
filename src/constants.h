#pragma once

using namespace Eigen;

static const double dx = 0.01;
static const int X_GRID = int(1 / dx);
static const int Y_GRID = int(1 / dx);
static const int X_WINDOW = 1400;
static const int Y_WINDOW = X_WINDOW * Y_GRID / X_GRID;

static const int n = 20;