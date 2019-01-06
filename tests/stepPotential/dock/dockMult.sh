


# scripts in workingdir, train and docking files in COMPLExsubfolder
#
# 
# for all in ComplexNameFile:
#   merge InputPoses with dockingResults (none at first)
#   make trainGrids
# train ->result in main-Folder
# for all in ComplexNameFiles:
#   sample.sh
#   copy best for next run
#args or defaults

#benchark according to benchamrk 5.0
set -u -e 
fParams=${1}
compDir=${2:-~/transfer/benchmark5_attract}  #files from attract-benchmark5-filder, poses have no unbound...
outDir=${4:-~/testNoH_dock}
fCompNames=${5:-`pwd`/complexesTest.txt}
fAtomTypes=${6:-`pwd`/atomtypes_gaa.dat}
scriptDir=`pwd`




useH=0 #recompile attract!! DO prepComp!
redo=1



# compDir=${1:-~/Uni/BA/benchmarks/benchmark5_attract}  #files from attract-benchmark5-filder, poses have no unbound...
# poseDir=${2:-~/Uni/BA/stepPotentials/Attract_benchmark_wDecoys}
# outDir=${3:-~/Uni/BA/benchmarks/test_multComplexes/test}
# fCompNames=${4:-~/Uni/BA/benchmarks/test_multComplexes/complexes.txt}
# fAtomTypes=${5:-~/Uni/BA/benchmarks/test_multComplexes/atomtypes_gaa.dat}
# scriptDir=`pwd`




###################prep
while read compName; do
echo "-----------preparing " $compName

#for smapling
  dSample=$outDir/$compName/sampling
  mkdir -p $dSample
  # delete 2x "-refe" for unbound docking
  if [ 1 -eq 1 ]; then ### move and change to disable parts of the protocol
  ./prepComp.sh $compDir/$compName/receptor-refe.pdb $compDir/$compName/ligand-refe.pdb $dSample $fAtomTypes $scriptDir $useH $redo
  fi

done <$fCompNames
##########################





    while read compName; do
    echo "-----------docking " $compName
    name=DOCK5_`basename $fParams .par`
        cd $outDir/$compName/sampling
        time $scriptDir/sample.sh $fParams $name

        python $scriptDir/capristars.py result.lrmsd result.irmsd result.fnat
        cp result.lrmsd $name.lrmsd;cp result.irmsd $name.irmsd;cp result.fnat $name.fnat;cp result.capstars $name.capstars; cp result.dat $name.dat
	cp native.score ${name}_native.score
    done <$fCompNames
