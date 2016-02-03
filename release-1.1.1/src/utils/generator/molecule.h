/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef MOLECULE_H
#define MOLECULE_H

/** molecule type
 *
 *  The molecule_t can store the data for a single molecule.
 */
typedef struct molecule_t {
    double r[3]; /**< coordinate of center of mass (x,y,z) */
    long cid; /**< component id */
} molecule_t;

#endif /* MOLECULE_H */
