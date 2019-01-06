#!/bin/bash
set -u -e 
$ATTRACTDIR/shm-clean


echo "setting params"
stepParam=${1}


#name of the run
name=${2:-attract_bench}

#docking parameters
#params="$ATTRACTDIR/../attract.par receptorr.pdb ligandr.pdb --fix-receptor"
params="$ATTRACTDIR/../attract.par receptorr.pdb ligandr.pdb --fix-receptor"
scoreparams="$ATTRACTDIR/../attract.par receptorr.pdb ligandr.pdb --score --fix-receptor"


#grid parameters
gridparams=" --grid 1 receptorgrid.grid"
start=systsearch.dat #init later


echo '**************************************************************'
echo 'Docking'
echo '**************************************************************'
rm -f out_${name}*[.dat,.lrmsd]
rm -rf result.dat result.pdb result.lrmsd result.irmsd result.fnat >& /dev/null


#$ATTRACTDIR/attract $start $params $gridparams --mc --mcscalerot 30 --mcscalecenter 3 --mcmax 1000 --gravity 2 --rstk 0.0001 --mcprobs 70 0 0 0 --step-pot $stepParam > out_$name.dat #vorher rstk 0.01
parals="--np 4 --chunks 4"  #chuncs can be more than cores! -> does 1 chunk at a time (or np chunks)
python $ATTRACTDIR/../protocols/attract.py $start $params $gridparams $parals --mc --mcscalerot 30 --mcscalecenter 3 --mcmax 1100 --gravity 2 --rstk 0.001 --mcprobs 70 0 0 0 --step-pot $stepParam --output out_$name.dat

# echo '**************************************************************'
# echo 'Scoring Native'
# echo '**************************************************************'
rm -f native.score
nativescoreparams="$ATTRACTDIR/../attract.par receptorr-nat.pdb ligandr-nat.pdb --score --fix-receptor" 

#touch native.score
$ATTRACTDIR/attract native.dat $nativescoreparams $gridparams  --step-pot $stepParam > native.score #for $gridparams: atomzahl in ordung bringen:


mv -f out_$name.dat out_$name-scored.dat

echo '**************************************************************'
echo 'Sort structures'
echo '**************************************************************'
python $ATTRACTTOOLS/sort.py out_$name-scored.dat > out_$name-sorted.dat

echo '**************************************************************'
echo 'remove redundant'
echo '**************************************************************'
$ATTRACTDIR/deredundant out_$name-sorted.dat 2 | python $ATTRACTTOOLS/fill-deredundant.py /dev/stdin out_$name-sorted.dat > out_$name-sorted-dr.dat #--lim 5 bei dered um alles mit rmsd 5 zu vereinen
#mv -f out_$name-sorted.dat out_$name-sorted-dr.dat

echo '**************************************************************'
echo 'Soft-link the final results'
echo '**************************************************************'
ln -f -s out_$name-sorted-dr.dat result.dat

echo '**************************************************************'
echo 'collect top 50 structures:'
echo '**************************************************************'
$ATTRACTTOOLS/top out_$name-sorted-dr.dat 50 > out_$name-top50.dat
# $ATTRACTDIR/collect out_$name-top50.dat receptor-sc.pdb ligand-sc.pdb > out_$name-top50.pdb #hier richtige sietenketten für end-pdb wählen, aus receptor-sc.pdb
$ATTRACTDIR/collect out_$name-top50.dat receptorr.pdb ligandr.pdb > out_$name-top50.pdb #hier richtige sietenketten für end-pdb wählen, aus receptor-sc.pdb

ln -f -s out_$name-top50.pdb result.pdb



echo '**************************************************************'
echo 'calculate backbone ligand RMSD'
echo '**************************************************************'
python $ATTRACTDIR/lrmsd.py out_$name-sorted-dr.dat ligand-aa.pdb ligand-refe.pdb --receptor receptor-aa.pdb > out_$name-sorted-dr.lrmsd
ln -s out_$name-sorted-dr.lrmsd result.lrmsd

echo '**************************************************************'
echo 'calculate backbone interface RMSD'
echo '**************************************************************'
python $ATTRACTDIR/irmsd.py out_$name-sorted-dr.dat receptor-heavy.pdb receptor-refeHeavy.pdb ligand-heavy.pdb ligand-refeHeavy.pdb  > out_$name-sorted-dr.irmsd
ln -s out_$name-sorted-dr.irmsd result.irmsd

echo '**************************************************************'
echo 'calculate fraction of native contacts'
echo '**************************************************************'
python $ATTRACTDIR/fnat.py out_$name-sorted-dr.dat 5 receptor-heavy.pdb receptor-refeHeavy.pdb ligand-heavy.pdb ligand-refeHeavy.pdb > out_$name-sorted-dr.fnat
ln -s out_$name-sorted-dr.fnat result.fnat

$ATTRACTDIR/shm-clean


