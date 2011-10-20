#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

#include "max.h"

char buf[100000];

extern "C" FILE *read_dof_init_(const char *f_, int nlig, int &line, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len) {
  srand((unsigned int)time(0));
  char f[1000];
  memcpy(f, f_, f_len);
  f[f_len] = 0;  
    
  FILE *fil = new FILE;
  fil = fopen(f, "r");
  line = 1;
  
  //read pivot header
  bool done = 0;
  bool ok = 1;
  while (!done) {
    done = 1;
    if(!fgets(buf,100000,fil))ok = 0;
    else if (strlen(buf) < 8 || strncmp(buf, "#pivot ", 7)) {
      if (strlen(buf) < 2 || buf[0] != '#' || buf[1] != '#') ok = 0;
      else done = 0;
    }
    if (!ok) {  
      fprintf(stderr, "Reading error in %s, line %d: expected pivot point or '#pivot auto'\n", f, line);
      exit(0);
    }
    if (!done) continue;
    if (!(strncmp(buf+7,"auto", 4))) {
      auto_pivot = 1;
    }
    else {
      auto_pivot = 0;
      bool ok = 1;
      for (int n = 0; n < nlig; n++) {
	int nr; float px, py, pz;
	if (ok) {
	  if (strncmp(buf, "#pivot ", 7)) ok = 0;
	  else if (sscanf(buf+7, "%d %f %f %f", &nr, &px, &py, &pz) < 4) ok = 0;
	  else if (nr != n+1) ok = 0;
	}
	if (!ok) {
          fprintf(stderr, "Reading error in %s, line %d: expected pivot point\n", f, line);
          exit(0);      
	}
	pivot[0][n] = px;
	pivot[1][n] = py;
	pivot[2][n] = pz;
	if (n < nlig-1) {
          ok = fgets(buf,100000,fil);
          line++;  
	}
      }    
    }
    if (done) break;    
  }
  
  bool cr_true = 0, cr_false = 0;
  if(fgets(buf,100000,fil)) {
    if (strlen(buf) >= 24 && !strncmp(buf, "#centered receptor: true", 24)) 
      cr_true = 1;
    if (strlen(buf) >= 25 && !strncmp(buf, "#centered receptor: false", 25)) 
      cr_false = 1;
  }    
  if (!(cr_true||cr_false)) {
    fprintf(stderr, "Reading error in %s, line %d: expected '#centered receptor: true' or '#centered receptor: false'\n", f, line);
    exit(0);
  }
  if (cr_true) centered_receptor = 1; 
  else centered_receptor = 0;

  bool cl_true = 0, cl_false = 0;
  if(fgets(buf,100000,fil)) {
    if (strlen(buf) >= 23 && !strncmp(buf, "#centered ligands: true", 23)) 
      cl_true = 1;
    if (strlen(buf) >= 24 && !strncmp(buf, "#centered ligands: false", 24)) 
      cl_false = 1;
    if (strlen(buf) >= 22 && !strncmp(buf, "#centered ligand: true", 22)) 
      cl_true = 1;
    if (strlen(buf) >= 23 && !strncmp(buf, "#centered ligand: false", 23)) 
      cl_false = 1;
  }    
  if (!(cl_true||cl_false)) {
    fprintf(stderr, "Reading error in %s, line %d: expected '#centered ligands: true' or '#centered ligands: false'\n", f, line);
    exit(0);
  }
  if (cl_true) centered_ligands = 1; 
  else centered_ligands = 0;
  
  return fil;
}

extern "C" int read_dof_(FILE *fil, int &line, int &nstruc, const char *f_, idof2 &ens, dof2 &phi, dof2 &ssi, dof2 &rot, dof2 &xa, dof2 &ya, dof2 &za, modes2 &dlig, const int &nlig, const int *nhm, const int *nrens0, int &seed, char *&label, int f_len) {
  int nrens00[MAXLIG];
  memset(nrens00,0,MAXLIG*sizeof(int));
  
  const int *nrens;
  if (nrens0 == NULL) {
    nrens = nrens00;
  }
  else nrens = nrens0;

  char f[1000];
  memcpy(f, f_, f_len);
  f[f_len] = 0;  

  int mode = 0;
  //modes: 
  // 0 = expecting comment
  // 1 = expecting data
  bool next = 1;
  int clig = 0;
  float fields[10000];
  int cseed = 0;
  char *currlabel = NULL;
  while (1) {
    if (feof(fil)) return 1;
    line++;
    if (next) {
      if(!fgets(buf,100000,fil)) {
        if (mode == 0) continue;
        fprintf(stderr, "Reading error in %s, line %d\n", f, line);
        exit(0);
      }
    }
    else next = 1;
    if (mode == 0) {
      if (buf[0] != '#') {
        fprintf(stderr, "Reading error in %s, line %d: expecting #\n", f, line);
	exit(1);
      }
      if (atoi(&buf[1]) == nstruc+1) {
        nstruc++;
        clig = 0;
        cseed = rand() % 1000000 + 1;
	label = NULL;
	mode = 1;
	continue;
      }
    }
    if (mode == 1) {
      char *seedstr = "### SEED ";
      if (strlen(buf) > strlen(seedstr) 
       && !strncmp(buf,seedstr,strlen(seedstr))) {	
	cseed = atoi(&buf[strlen(seedstr)]);
      }
      seed = cseed;
      if (strlen(buf) > strlen(seedstr) 
       && !strncmp(buf,seedstr,strlen(seedstr))) {	
	continue;
      }
      if (strlen(buf) >= 3 && buf[0] == '#' && buf[1] == '#' && buf[2] == '#') {
        if (currlabel == NULL) { 
	  currlabel = new char[10000];
	  currlabel[0] = 0;               
	}
	strcat(currlabel, buf);
      }
      if (buf[0] == '#' && buf[1] == '#') continue;
      if (buf[0] == '#') {
        fprintf(stderr, "Reading error in %s, line %d: expecting data\n", f, line);
	exit(1);
      }
      char chars[] = {32,10,13,0};
      int firstfield;
      char *field = strtok(buf,chars);
      int nf = 0;
      while (field != NULL) {
	if (nf == 10000) {
   	  fprintf(stderr, "Reading error in %s, line %d: too many values on a line\n", f, line);
   	  exit(1);
	}
	if (nf == 0) firstfield = atoi(field);
	fields[nf] = atof(field);
	//printf("%s\n", field);
	field = strtok(NULL,chars);
	nf++;
      }
      int has_ens = 0;
      if (nrens[clig]) has_ens = 1;
      if (nf < has_ens + 6 || nf > has_ens + 6 + nhm[clig]) {
   	  fprintf(stderr, "Reading error in %s, line %d: read %d values, expected %d values\n", f, line, nf, has_ens + 6 + nhm[clig]);
	  exit(1);
      }
      memset(ens, 0, MAXLIG*sizeof(int));
      if (has_ens) {
        ens[clig] = firstfield; //fields[0]
	if (firstfield <= 0 || firstfield > nrens[clig]) {
   	  fprintf(stderr, "Reading error in %s, line %d: ensemble copy must be in range 1-%d, is %d\n", f, line, nrens[clig], firstfield);
	  exit(1);
	}
      }
      phi[clig] = fields[has_ens+0];
      ssi[clig] = fields[has_ens+1];
      rot[clig] = fields[has_ens+2];
      xa[clig] = fields[has_ens+3];
      ya[clig] = fields[has_ens+4];
      za[clig] = fields[has_ens+5];
      for (int n=0;n<nhm[clig];n++) {
        if (has_ens+6+n < nf) {
          dlig[clig][n] = fields[has_ens+6+n];
	}
	else 
	  dlig[clig][n] = 0;
      }
      
      clig++;
      if (clig == nlig) {        
        mode = 0;
	if (currlabel != NULL) {
	  label = new char[strlen(currlabel)+1];
	  strcpy(label, currlabel);
	}
	return 0;
      }
    }
  }
}

