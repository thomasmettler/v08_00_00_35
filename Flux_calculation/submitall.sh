#for f in /pnfs/uboone/data/uboone/raw/online/crt_seb/v2_0/00/00/00/00/ProdRun201705??_*-crt0?.1.crtdaq
#for f in /pnfs/uboone/data/uboone/raw/online/crt_seb/v2_0/00/00/00/00/ProdRun201704??_*-crt0?.1.crtdaq
COUNTER=0
echo "File in $1"

for f in $(cat $1)
do
  echo "Processing Nr: $COUNTER, file: $f"
  #echo $(basename $f)
  let COUNTER=COUNTER+1
  # take action on each file. $f store current file name
#jobsub_submit -G uboone --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL7 --expected-lifetime='1d' --disk=30GB --role=Analysis "file://$(pwd)/make_fluxhisto.sh $(basename $f)"
jobsub_submit -G uboone --role=Analysis -N 1 --resource-provides="usage_model=DEDICATED,OPPORTUNISTIC" --OS="SL7" --expected-lifetime='1d' --disk=30GB file://$(pwd)/make_fluxhisto.sh $f
#jobsub_submit -G uboone --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 --expected-lifetime='1d' --role=Analysis file://exbeam.sh $(basename $f)
#jobsub_submit -G uboone --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC  --expected-lifetime='3d' --disk=130GB --role=Analysis file://exbeam2.sh $(basename $f)

done

#ProdRun20170515_041007-crt01.1.crtdaq
