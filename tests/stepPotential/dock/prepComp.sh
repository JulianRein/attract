#transform pdbs + make grid
 #for all in CompNameFile: ->maybe for all later -> for all... not here, only args for each
#   make working dir with subfolder
#   link poses as COMPLEXstartposes.dat
#   bound structs to ligandr, receptorr
#   make grids for receptorr
set -e -x

rec=$1
lig=$2
outDir=$3
fAtomTypes=$4
dScripts=$5 
useH=$6
redo=$7

cd $outDir

if [ $useH -eq 1 ] ; then
    if   [ ! -f recH.pdb ] || [ $redo -eq 1 ]; then
        pdb2pqr --ff=amber $rec recH.pdb 
        sed -i 's/ H / HN/g' recH.pdb #only "HN" will be recognized by redatom.py as the H to N in backbone
        pdb2pqr --ff=amber $lig ligH.pdb 
        sed -i 's/ H / HN/g' ligH.pdb
    fi
    toTransformRec=recH.pdb
    toTransformLig=ligH.pdb
    else
    toTransformRec=$rec
    toTransformLig=$lig
fi

python $dScripts/redatom.py --atomtypefile $fAtomTypes --pdb $toTransformRec --output recGAA.pdb
sort -k 1,1 -k1.24,1.26n -k 1.17,1.17 recGAA.pdb -o recGAA.pdb #sort ATOM,ResNr,copyId ->maybe add AtomType gaa...


python $dScripts/redatom.py --atomtypefile $fAtomTypes --pdb $toTransformLig --output ligGAA.pdb
sort -k 1,1 -k1.24,1.26n -k 1.17,1.17 ligGAA.pdb -o ligGAA.pdb

mv recGAA.pdb receptorr.pdb
mv ligGAA.pdb ligandr.pdb

#rmsd
ln -sf $rec receptor-aa.pdb
ln -sf $lig ligand-aa.pdb
#native-scoring -> TODO: for unbound-docking always from bound state...
ln -sf receptorr.pdb receptorr-nat.pdb
ln -sf ligandr.pdb ligandr-nat.pdb
ln -sf $dScripts/native.dat native.dat

#reference structures
#change later if unbound docking...
ln -sf receptor-aa.pdb receptor-refe.pdb
ln -sf ligand-aa.pdb ligand-refe.pdb
python $ATTRACTDIR/../allatom/aareduce.py receptor-refe.pdb receptor-refeHeavy.pdb --heavy --pdb2pqr > /dev/null
python $ATTRACTDIR/../allatom/aareduce.py ligand-refe.pdb ligand-refeHeavy.pdb --heavy --pdb2pqr > /dev/null
ln -sf receptor-refeHeavy.pdb receptor-heavy.pdb
ln -sf ligand-refeHeavy.pdb ligand-heavy.pdb

if [ ! -f receptorgrid.grid ] || [ $redo -eq 1 ] ; then
awk '{print substr($0,58,2)}' ligandr.pdb | sort -nu > receptorgrid.alphabet
$ATTRACTDIR/make-grid-omp receptorr.pdb $ATTRACTDIR/../attract.par 7 7 receptorgrid.grid  --alphabet receptorgrid.alphabet
fi 


start=systsearch.dat

if [ ! -f $start ] ; then
echo '**************************************************************'
echo 'Generate starting structures...'
echo '**************************************************************'

# python $ATTRACTDIR/../tools/randsearch.py 2 5000 receptor-aa.pdb --radius 80 --fix-receptor >systsearch.dat #does not take receptor-size into account
cat $ATTRACTDIR/../rotation.dat > rotation.dat
$ATTRACTDIR/translate receptorr.pdb ligandr.pdb > translate.dat
$ATTRACTDIR/systsearch > $start
# 
# 
#add pivot of pose-files to startPositions
#echo "`head -n 2 ../training/start.dat``tail -n +2 systsearch.dat`" >systsearch.dat #in one line
# mv $start systsearch.dat.tmp
# head -n 2 ../training/start.dat >$start
# tail -n +2 systsearch.dat.tmp >>$start
fi
