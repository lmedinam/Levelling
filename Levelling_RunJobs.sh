HOMEDIR=/afs/cern.ch/user/l/lmedinam
PROGRAMDIR=$HOMEDIR/Level_2

# CHECK LIKES 39 and 61 for the Luminosity Version!

cd
cd $PROGRAMDIR
echo ""
echo -n "START DIRECTORY: "
pwd

FILE="Levelling_List.txt"
ARRAY=( $(<"$FILE") )

echo ""

for (( i=0; i<${#ARRAY[@]}; i++ )); do

  echo ""
  echo "======================================================================"
  echo "Scenario:" ${ARRAY[i]}
  echo "======================================================================"
  echo ""

  # GENERIC NAME

  NAME=${ARRAY[i]}

  # CREATE A FOLDER FOR THE RESULTS BY MAKING A COPY OF THE ORIGINAL PROGRAM DIRECTORY:

  echo ""
  echo -n "DIRECTORY WITH RESULTS: "

  RESULTDIR=$HOMEDIR/Level_2_results_$NAME
  
  mkdir $RESULTDIR
  
  cp Levelling_Beam.py $RESULTDIR
  cp Levelling_Config.py $RESULTDIR
  cp Levelling_Densities.py $RESULTDIR
  cp Levelling_jobtau.madx $RESULTDIR
  #cp Levelling_jobtau-old.madx $RESULTDIR
	#cp Levelling_jobtau_helhc.madx $RESULTDIR
  cp Levelling_Luminosity_v062.py $RESULTDIR
  cp Levelling_Others.py $RESULTDIR
  cp Levelling_Run.py $RESULTDIR
  cp HLLHC-sb150.seq $RESULTDIR

  cd $RESULTDIR 

  pwd

  echo ""
  echo "FILES IN DIRECTORY:"
  echo ""
  
  ls -l
  
  echo ""

bsub -q 8nh  -J "LEV[$((i+1))]" << EOF
#bsub -q 1nd  -J "LEV2[$((i+1))]" << EOF
#bsub -q 2nd -J "LEV3[$((i+1))]" << EOF
echo ""
pwd
cp -r $RESULTDIR/* .
/afs/cern.ch/user/l/lmedinam/anaconda/bin/python Levelling_Run.py -n ${ARRAY[i]}
ls -l
cp -r * $RESULTDIR/
cd $RESULTDIR
rm db5
rm Levelling*
rm slhc
rm twiss*
#rm 1*out
cd
EOF

done

cd $HOMEDIR

echo ""
