#ifndef AMR_LATTICE_H
#define AMR_LATTICE_H

// Directions: 1 = right, 2 = left, 3 = up, 4 = down.
// Sites are numbered 1 through size*size.
int lattice_neighbour(int site, int direction, int size);

// Store empty neighbours in direction order: right, left, up, down.
// Returns the number stored. empty must have room for four entries.
int lattice_empty_neighbours(int site, int size, const int *occupied, int empty[4]);

#endif