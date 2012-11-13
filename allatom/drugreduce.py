break_threshold = 3.5 #3.5 A before a N/C-terminal patch is applied
threshsq = break_threshold * break_threshold

import sys, os

import parse_cns_top 


def run(pdblines, topfile, transfile, outf, mapf, patches):
  global rescounter, atomcounter

  rescounter = 1
  atomcounter = 1  
  residues, presidues = parse_cns_top.parse_stream(open(topfile))

  code_to_type = {}
  for l in open(transfile):
    ll = l.split()
    type = int(ll[0])
    for code in ll[3:]:
      if code.startswith("#"): break
      assert code not in code_to_type, code
      code_to_type[code] = type

  class PDBres:
    def __init__(self, resid, resname):
      self.resid = resid
      self.resname = resname
      self.coords = {} 
      self.ter = False

  pdbres = []

  def write_res(res, resprev=None, resnext=None):
    global rescounter, atomcounter

    r = residues[res.resname.lower()].copy()
    if res.resid in patches:
      p = patches[res.resid]
      if isinstance(p, str): p = [p]
      for pp in p:
        r.patch(presidues[pp.lower()])

    done = set()
    for a in r.atomorder: 
      atom = r.atoms[a]
      if a.startswith("h") and atom.charge == 0: continue
      aa = a.upper()
      x = " XXXXXXX"
      y = x; z = x
      if aa.strip() in res.coords:
        x,y,z = ("%8.3f" % v for v in res.coords[aa.strip()])
      xyz = x + y + z
      type = code_to_type[atom.type.upper()]
      a0 = aa
      if len(a0) < 4:
        a0 = " " + a0 + "   "[len(a0):]
      print >> outf, "ATOM %6d %4s %s %5d    %s %4d %7.3f 0 1.00" % \
       (atomcounter, a0, res.resname, rescounter, xyz, type, atom.charge )
      atomcounter += 1
    print >> mapf, res.resid.split()[-1], rescounter
    rescounter += 1


  curr_res = None

  for l in pdblines:
    if l.startswith("TER"):
      if curr_res is not None: curr_res.ter = True
      continue
    if l.startswith("ATOM") or l.startswith("HETATM"):
      atomcode = l[12:16].strip()
      assert l[16] == " ", l
      resname = l[17:20]   
      resid = l[21:27]
      x = float(l[30:38])
      y = float(l[38:46])
      z = float(l[46:54])
      if curr_res is None \
       or (resid != curr_res.resid) or resname != curr_res.resname:
        curr_res = PDBres(resid, resname)
        pdbres.append(curr_res)
      curr_res.coords[atomcode] = (x,y,z)

  for n in range(len(pdbres)):
    res = pdbres[n]
    resprev = None
    if n > 0: resprev = pdbres[n-1]
    resnext = None
    if n < len(pdbres)-1: resnext = pdbres[n+1]
    write_res(res, resprev, resnext)
  print >> outf, "END"
  
if __name__ == "__main__":
  pdb = sys.argv[1]
  pdblines = open(pdb).readlines()
  transfile = sys.argv[2]
  topfile = sys.argv[3]

  outfile = os.path.splitext(pdb)[0] + "-aa.pdb"
  outf = open(outfile, "w")

  run(pdblines, topfile, transfile, outf, sys.stdout, {})
