#ifndef AMR_LATTICE_H
#define AMR_LATTICE_H

// Directions: 1 = right, 2 = left, 3 = up, 4 = down.
// Sites are numbered 1 through size*size.
int lattice_neighbour(int site, int direction, int size);

#endif