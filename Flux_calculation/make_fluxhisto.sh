####################################
###### setup your needed products here, e.g. geant4 etc...
####################################

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v08_48_00 -q e19:prof

####################################
#### This is where you copy all of your executable/necessary files to the worker node 
#### ( If applicable )
####################################

###### this is where you copy your executable - I have a simple hello.out code here.
# ifdh cp /pnfs/uboone/scratch/users/kreslo/beamextract beamextract
 echo "Start coping files"
 ifdh cp $1 fluxreader.root
 echo "copied flux root file"
 ifdh cp /pnfs/uboone/scratch/users/tmettler/v08_48_00/plot plot
 echo "copied plot executable"

####### 
####### ifdh cp does not preserve permissions, so need to add executable. #########
#######
 chmod u+x plot
#######
####### launch executable
#######

# ./beamextract
 echo "Start running over flux file"
#root -l -b -q 'CoincScale_V8.C+g("'$1'",1)'  >/dev/null
 ./plot fluxreader.root
#######
####### Copy results back 
#######

echo "Finished process, start copy back"
#ifdh cp *.root /pnfs/uboone/scratch/users/kreslo/datamay/
#ifdh mkdir ${SCRATCH_DIR}/${GRID_USER}/output_${CLUSTER}.${PROCESS}
echo "Made folder, start copy"
ifdh cp *.root $(dirname "$1")/
echo "ifdh cp *.root ${SCRATCH_DIR}/${GRID_USER}/output_${CLUSTER}.${PROCESS}/"
echo "Finished everything"

