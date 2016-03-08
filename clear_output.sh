#! /bin/bash

## Prompt the user for confirmation that they actually want to
## obliterate all of the output

echo "This will remove all snooppy output (ALL of it)."
echo "Are you REALLY SURE you want to proceed? (yes/NO)?"
read YNResponse

if [ "$YNResponse" == "yes" ]
then
  rm -f time* *.dat dump* output* nohup* ./data/*
else
  echo "Taking no action."
fi
