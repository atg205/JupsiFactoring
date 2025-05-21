#! /usr/bin/bash
N="493"
mkdir "energies/$N"
for filepath in IsingTerms/$N/*; do 
    filename=$(basename "$filepath")
	BLOCKSIZE=16 ../../xubo/xubo_ising IsingTerms/$N/$filename  | tee energies/$N/$filename.energies;
	 

done
