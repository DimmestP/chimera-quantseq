mkdir indexed_genome
for file in *.fa
do
	hisat2-build $file "./indexed_genome/${file%.*}" 
done
