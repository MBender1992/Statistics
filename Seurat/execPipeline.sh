#!/bin/bash

# bash script to sequentially execute R scripts for scRNAseq

# directory of R scripts
SCRIPT_DIR="/run/user/1001/gvfs/smb-share:server=192.168.35.240,share=zellbiologie/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/Github/Statistics/Seurat/R"

# check if directory exists
if [ ! -d "$SCRIPT_DIR" ]; then
  echo "Error: Directory $SCRIPT_DIR does not exist."
  exit 1
fi

# define all scripts which start with a number
SCRIPTS=$(find "$SCRIPT_DIR" -maxdepth 1 -type f -name "[0-9][1-9]*.R" | sort)

echo "Found the following scripts:"
echo "$SCRIPTS"  # This will display the list of scripts found

# loop through scripts
for SCRIPT in $SCRIPTS
do
  echo "Running $SCRIPT ..."
  Rscript "$SCRIPT" || { echo "Error running $SCRIPT"; exit 1; }
  if [ $? -ne 0 ]; then
    echo "Error running $SCRIPT. Exiting."
    exit 1
  fi
done

echo "All scripts executed successfully!"



