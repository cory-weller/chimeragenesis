#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-3:00:00
#SBATCH --account sadhumj
#SBATCH -p norm

mkdir -p uniprot_db

function get_related_structures() {
  cluster=${1}
  pdb_id=${2}
  species=${3}
  wget -q -O - "http://www.uniprot.org/uniprot/${pdb_id}.txt" > uniprot_db/${pdb_id}.uniprot && sleep 0.35
}

while read cluster bitscore species inparanoidscore pdb_id pidentity; do
  get_related_structures ${cluster} ${pdb_id} ${species}
done < sqltable.H.sapiens-S.cerevisiae
