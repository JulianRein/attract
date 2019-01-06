#      program redatom
import sys     
#  run by:python redatom.py input.pdb >out.pdb

if '--atomtypefile' in sys.argv:
  atomtypes = sys.argv[sys.argv.index('--atomtypefile')+1]
else:
  print 'please give file with all the atomtypes for the model'
  sys.exit()
  
build_atom = open(atomtypes, 'r')
atomTypes=build_atom.readlines()

if '--output' in sys.argv:
  outname=sys.argv[sys.argv.index('--output')+1]
else:
  print 'give outputname'
  sys.exit()

if '--pdb' in sys.argv:
  eingabe=sys.argv[sys.argv.index('--pdb')+1]
else: 
  print 'please give pdb to convert'
  sys.exit()
  
  
  
if '--addNH' in sys.argv:
  f=open(eingabe, 'r') 
  i=0
  iflater=0
  pdbLines=[]
  for line in f:
      if (line[:4]=='ATOM') and line[13:16]!='OXT' and line[13:14]!='H' and line[12:13]!='H':
	  pdbLines.append(line)
	  i+=1

  f.close()
  # automatic addition of hydrogens
  # treat first N separately
  natom=i
  nnew=0
  for a in range(natom):
      if pdbLines[a][13:16]=='N  ' and pdbLines[a][17:20]!='PRO':
	  nnew+=1
  for a in range(natom+nnew):
      if pdbLines[a][13:16]=='N  ' and pdbLines[a][17:20]!='PRO':
	  pdbLines.insert(a+1,pdbLines[a][:13]+'%-3s' % 'H'+pdbLines[a][16:30]+'%8s' % round(float(pdbLines[a][30:38])*
			1.66-float(pdbLines[a+1][30:38])*0.66,3)+'%8s' % round(float(pdbLines[a][38:46])*1.66-float(pdbLines[a+1][38:46])*0.66,3)+
			'%8s' % round(float(pdbLines[a][46:54])*1.66-float(pdbLines[a+1][46:54])*0.66,3)+pdbLines[a][54:69]+'\n')

else:
  fobj=open(eingabe, 'r')
  lines=fobj.readlines()
  i=0
  pdbLines=[]
  for line in lines:
    if line[:4]=='ATOM':
      pdbLines.append(line)
      i+=1
  nnew=0
  natom=i
  fobj.close()
#
#     generierung der pdbcade tabelle und der integercode tabellen
#     fuer jede Aminosaeure;
#
if '--atomtypesonly' in sys.argv:
  outobj=open(outname,'w+')
  t=0
  atm=natom
  test=0
  for i in range(natom):
      for aType in range(len(atomTypes)):
	  if atomTypes[aType][:4]==pdbLines[i][17:21] and atomTypes[aType][6:10]==pdbLines[i][12:16]:
	      test+=1
	      outobj.write(pdbLines[i][:55]+atomTypes[aType][12:16]+pdbLines[i][59:67]+" "+"0"+" "+"1.00"+"\n")

      if pdbLines[i][:3]=='TER' or pdbLines[i][:3]=='END':
	  outobj.write(pdbLines[i])
	  t=0
	  atm=0
	  test=0

  outobj.close()
else:  
  outobj=open(outname,'w+')
  t=0
  atm=0
  test=0
  c_sidechainAtoms=0
  #go through every atom and test every atomtype ->write out for match...
  for atom in range(natom+nnew):
      #atm+=1 not working here with new reduce/addHydrogens, too many h-atoms that are removed

      if atom<(natom+nnew-1):
	  if pdbLines[atom][21:26]!=pdbLines[atom-1][21:26]: #new res
	      t+=1
      for aType in range(len(atomTypes)):
	  if atomTypes[aType][:4]==pdbLines[atom][17:21] and atomTypes[aType][6:10]==pdbLines[atom][12:16]:
	      atm+=1
	      test+=1
	      c_copy=0
	      if pdbLines[atom][16]==" ": c_sidechainAtoms=0 #useless
	      elif pdbLines[atom][16]=="A": c_sidechainAtoms+=1 #useless
	      else:
	        c_copy=ord(pdbLines[atom][16])-ord("A")
	        #if pdbLines[atom][16]!=pdbLines[atom-1][16]: #new copy, reset counter... ; now useless, will be incremented even if copy
	          #atm-=c_sidechainAtoms
	          #test-=c_sidechainAtoms

	      outobj.write(pdbLines[atom][:6]+'%5s' % atm+pdbLines[atom][11:21]+'%5s' % t+pdbLines[atom][26:55]+atomTypes[aType][12:16]+atomTypes[aType][17:25]+'%2s' % c_copy+" 1.00"+"\n") #upgraded "new res-detection"
	      #outobj.write(pdbLines[atom][:6]+'%5s' % atm+pdbLines[atom][11:21]+pdbLines[atom][21:26]+pdbLines[atom][26:55]+atomTypes[aType][12:16]+atomTypes[aType][17:25]+'%2s' % c_copy+" 1.00"+"\n") #use old resNr instead of newly computed one, does not rely on "N" first instead of maybe NH

      if pdbLines[atom][:3]=='TER' or pdbLines[atom][:3]=='END':
	  outobj.write(pdbLines[atom])
	  t=0
	  atm=0
	  test=0

  outobj.close()
if test!=atm:
  print 'something went wrong, wrote: ', test, atm
	
