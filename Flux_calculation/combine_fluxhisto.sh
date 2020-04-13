


COUNTER=0
echo "File in $1 and output in $2"

for f in $(cat $1)
do
  echo "Processing $f file..."
  if [ $COUNTER -gt 0 ]
  then
    echo "Adding file nr: $COUNTER"
    #echo $(basename $f)
    mv $2 $2_tmp.root
    hadd $2 $2_tmp.root $f
  else
    echo "Copying file nr: $f"
    cp $f $2
  fi
  let COUNTER=COUNTER+1
done

rm $2_tmp.root