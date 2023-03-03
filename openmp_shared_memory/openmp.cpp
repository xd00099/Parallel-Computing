#include "common.h"
#include <omp.h>
#include <vector>
#include <iostream>
#include <cmath>

// Put any static global variables here that you will use throughout the simulation.

double bin_size = cutoff * 2.1;
int grid_size;
const int neighbors[4][2] = {{0,1},{1,0},{1,-1},{1,1}};
std::vector<particle_t*> **grid;
// set locks for the bins
omp_lock_t **locks;

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
#pragma omp atomic
    particle.ax += coef * dx;
#pragma omp atomic
    particle.ay += coef * dy;
#pragma omp atomic
    neighbor.ax -= coef * dx;
#pragma omp atomic
    neighbor.ay -= coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    double p_x_prev = p.x;
    double p_y_prev = p.y;
    int prev_bin_x = p.x / bin_size;
    int prev_bin_y = p.y / bin_size;
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }

    // move to bins here
    int cur_bin_x = p.x / bin_size;
    int cur_bin_y = p.y / bin_size;

    // doesn't need to move
    if ((prev_bin_x == cur_bin_x) && (prev_bin_y == cur_bin_y)) {
        return;
    }

    // set lock before edit bins
    omp_set_lock(&locks[prev_bin_x][prev_bin_y]);

    // delete the previous point
    std::vector<particle_t*> &prev_grid = grid[prev_bin_x][prev_bin_y];
    for (int i = 0; i < prev_grid.size(); i++) {
        if (prev_grid[i] == &p){
            prev_grid.erase(prev_grid.begin() + i);
            break;
        }
    }
    omp_unset_lock(&locks[prev_bin_x][prev_bin_y]);

    // add to the new grid
    omp_set_lock(&locks[cur_bin_x][cur_bin_y]);
    grid[cur_bin_x][cur_bin_y].push_back(&p);
    omp_unset_lock(&locks[cur_bin_x][cur_bin_y]);

}


void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here

	grid_size = floor(size / bin_size) + 1;
    grid = new std::vector<particle_t*>*[grid_size];
    locks = new omp_lock_t *[grid_size];
    for (int i = 0; i < grid_size; i++) {
        grid[i] = new std::vector<particle_t*>[grid_size];
        locks[i] = new omp_lock_t[grid_size];

        for (int j = 0; j < grid_size; j++) {
            omp_init_lock(&locks[i][j]);
        }
    }

    // Initialize particles:
    for (int i = 0; i < num_parts; ++i) {
        // Add particle to grid:
        int x = parts[i].x / bin_size;
        int y = parts[i].y / bin_size;
        grid[x][y].push_back(&parts[i]);
    }
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // reset accelerations
	
	#pragma omp for
    for(int i = 0; i < num_parts; ++i){
        parts[i].ax = 0;
        parts[i].ay = 0;
    }
	
    // Compute Forces
    // traverse thru each grid
	#pragma omp for collapse(2)
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            
            // traverse thru each particle in the grid
            for (int k = 0; k < grid[i][j].size(); ++k) {
                particle_t *this_p = grid[i][j][k];

                // applying forces to each other in the same grid
                for (int l = k+1; l < grid[i][j].size(); ++l) {
                    apply_force(*this_p, *grid[i][j][l]);
                }

                // applying forces from neighboring grids
                for (int a = 0; a < 4; a++) {
                    int neighbor_i = i + neighbors[a][0];
                    int neighbor_j = j + neighbors[a][1];
                    if (neighbor_i < 0 || neighbor_i >= grid_size || neighbor_j < 0 || neighbor_j >= grid_size) {
                        continue;
                    }
                    for (int b = 0; b < grid[neighbor_i][neighbor_j].size(); ++b) {
                        apply_force(*this_p, *grid[neighbor_i][neighbor_j][b]);
                    }
                }
            }
        }
    }
    
    // Move Particles
	#pragma omp for
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
}
