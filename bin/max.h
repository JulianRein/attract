#ifndef MAX_H
#define MAX_H /* to make sure that we include only once... */

const int MAXATOM = 10000;
const int MAXRES = 3000;
const int TOTMAXATOM = 100000;
const int TOTMAXRES = 10000;
const int MAXLIG = 100;
const int MAXMODE = 10;
const int MAXMOLPAIR = 2000000;
const int MAXDOF = 1070; //100 ligands with modes and index modes, use 2 ligands for index modes only
const int MAXATOMTYPES = 99;
const int MAXATOMTYPES_EXT=127; // extended atomtypes, for gaa-atoms with step potential
const int MAXSELECTION = 5000; //maximum size of selection; NOTE: a static array of MAXSELECTION*MAXSELECTION Coors+doubles is kept in memory!
const int MAXRESTRAINTS = 10000;
const int MAXENS = 3000; //maximum ensemble size
const int MAXLENINDEXMODE = 10; //maximum number of nonzero entries in index modes
const int MAXINDEXMODE = 1000; // maximum number of index modes for flexible interface
const int MAXGRIDNEIGHBOUR = 200000000; //maximum size of total neighbour grid; 10 bytes per item

const int MAXCOPIES=31; //ALSO CHANGE IN read_pdb.f, max number of allowed copies per residue, used by cartsate.ncop -> more leads to overwrite!
const int MAXMOVER = 4; // monte carlo mover number; rigid_body & ensemble & both & side-chain switching now

typedef char char4[4];
typedef double dof2[MAXLIG];
typedef int idof2[MAXLIG];
typedef double coors2[MAXLIG][3];
typedef double modes2[MAXLIG][MAXMODE+MAXINDEXMODE];

typedef double Parameters[MAXATOMTYPES][MAXATOMTYPES];
typedef int iParameters[MAXATOMTYPES][MAXATOMTYPES];

#endif
