#include "common.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

double bin_size = cutoff * 2;
int grid_size;
const int neighbors[4][2] = {{0,1},{1,0},{1,-1},{1,1}};
std::vector<particle_t*> **grid;

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
    particle.ax += coef * dx;
    particle.ay += coef * dy;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
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
}


void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here

    // Define a 2D array in which each cell is a list of particles:
    // int grid_size = 1 / cutoff;
    // // Define global grid:
    // std::vector<particle_t*> grid[grid_size][grid_size];

    grid_size = floor(size / bin_size) + 1;
    grid = new std::vector<particle_t*>*[grid_size];
    for (int i = 0; i < grid_size; i++) {
        grid[i] = new std::vector<particle_t*>[grid_size];
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
    for(int i = 0; i < num_parts; ++i){
        parts[i].ax = 0;
        parts[i].ay = 0;
    }
    
    // Compute Forces
    // traverse thru each grid
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            
            // traverse thru each particle in the grid
            for (int k = 0; k < grid[i][j].size(); ++k) {
                particle_t *this_p = grid[i][j][k];

                // applying forces to each other in the same grid
                for (int l = k+1; l < grid[i][j].size(); ++l) {
                    apply_force(*this_p, *grid[i][j][l]);
                }

                // commented: applying to all 8 neighbors 
                // for (int a = i - 1; a <= i + 1; ++a) {
                //     for (int b = j - 1; b <= j + 1; ++b) {
                //         if (a < 0 || a >= grid_size || b < 0 || b >= grid_size || (a == i && b == j)) {
                //                 continue;
                //             }
                //         for (int c = 0; c < grid[a][b].size(); ++c) {
                //             apply_force(*this_p, *grid[a][b][c]);
                //         }
                //     }
                // }

                // applying forces from neighboring grids: only 4 neighbors is needed
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
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
    // Erase vectors:
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            grid[i][j].clear();
        }
    }
    for (int i = 0; i < num_parts; ++i) {
        // Add particle to grid:
        int x = parts[i].x / bin_size;
        int y = parts[i].y / bin_size;
        grid[x][y].push_back(&parts[i]);
    }
}

// Below is the original O(n^2) code
// #include "common.h"
// #include <cmath>

// // Apply the force from neighbor to particle
// void apply_force(particle_t& particle, particle_t& neighbor) {
//     // Calculate Distance
//     double dx = neighbor.x - particle.x;
//     double dy = neighbor.y - particle.y;
//     double r2 = dx * dx + dy * dy;

//     // Check if the two particles should interact
//     if (r2 > cutoff * cutoff)
//         return;

//     r2 = fmax(r2, min_r * min_r);
//     double r = sqrt(r2);

//     // Very simple short-range repulsive force
//     double coef = (1 - cutoff / r) / r2 / mass;
//     particle.ax += coef * dx;
//     particle.ay += coef * dy;
// }

// // Integrate the ODE
// void move(particle_t& p, double size) {
//     // Slightly simplified Velocity Verlet integration
//     // Conserves energy better than explicit Euler method
//     p.vx += p.ax * dt;
//     p.vy += p.ay * dt;
//     p.x += p.vx * dt;
//     p.y += p.vy * dt;

//     // Bounce from walls
//     while (p.x < 0 || p.x > size) {
//         p.x = p.x < 0 ? -p.x : 2 * size - p.x;
//         p.vx = -p.vx;
//     }

//     while (p.y < 0 || p.y > size) {
//         p.y = p.y < 0 ? -p.y : 2 * size - p.y;
//         p.vy = -p.vy;
//     }
// }


// void init_simulation(particle_t* parts, int num_parts, double size) {
// 	// You can use this space to initialize static, global data objects
//     // that you may need. This function will be called once before the
//     // algorithm begins. Do not do any particle simulation here
// }

// void simulate_one_step(particle_t* parts, int num_parts, double size) {
//     // Compute Forces
//     for (int i = 0; i < num_parts; ++i) {
//         parts[i].ax = parts[i].ay = 0;
//         for (int j = 0; j < num_parts; ++j) {
//             apply_force(parts[i], parts[j]);
//         }
//     }

//     // Move Particles
//     for (int i = 0; i < num_parts; ++i) {
//         move(parts[i], size);
//     }
// }