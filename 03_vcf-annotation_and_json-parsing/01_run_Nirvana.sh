##### 01_run_Nirvana.sh #####
# Alysson Henrique Urbanski & Tiago Minuzzi Freire da Fontoura Gomes @ Fiocruz Paraná, 2026 Jul.
# Nirvana annotation of aPRoVAR multisample VCF  

#!/bin/bash

INPUT="$1"
OUTPUT="${INPUT%.vcf.gz}"

/opt/dragen/4.3.13/share/nirvana/Nirvana -c /staging/NirvanaData/Cache/GRCh38/ \
					 -r /staging/NirvanaData/References/Homo_sapiens.GRCh38.Nirvana.dat \
					 --sd /staging/NirvanaData/SupplementaryAnnotation/GRCh38/ \
					 -i $INPUT \
					 -o $OUTPUT \
					 --vcf-info F_MISSING,AF,AF_EXCL
