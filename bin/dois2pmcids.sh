#!/bin/bash
#@author: madeline

#User input arguments are --dois (filepath to txt file of dois to look for) and --out_file (filepath to final output file, a .csv)
# eg. ./bin/dois2pmcids.sh -o out_dois.csv -i test_dois.txt

if [ "$1" = "--dois" -o "$1" = "-i" ]; then
    DOIS=$2
elif [ "$1" = "--out_file" -o "$1" = "-o" ]; then
    OUT=$2
else
    echo "Both --dois (-i) and --out_file (-o) are required."
fi

if [ "$3" = "--dois" -o "$3" = "-i" ]; then
    DOIS=$4
elif [ "$3" = "--out_file" -o "$3" = "-o" ]; then
    OUT=$4
else
    echo "Both --dois (-i) and --out_file (-o) are required."

fi



# create a csv that's just the header, for the outputs of curl for this mutation to be concatenated to
echo '"PMID","PMCID","DOI","Version","MID","IsCurrent","IsLive","ReleaseDate","Msg"' > $OUT

#read dois into a bash array
IFS=$'\n' read -d '' -r -a dois_arr < $DOIS

# query group_size records at a time and append csv chunks to the main csv
group_size=50

for((i=0; i < ${#dois_arr[@]}; i+=group_size))
do
  dois_arr_chunk="${dois_arr[@]:i:group_size}"
  id_string=$(echo $dois_arr_chunk | sed 's/ /,/g')
  curl -X GET "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?format=csv&email=m.iseminger@alumni.ubc.ca&tool=VarRec&idtype=doi&ids=${id_string}" | sed 1d >> $OUT
done
