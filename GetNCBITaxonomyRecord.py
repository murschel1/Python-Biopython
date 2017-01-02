############################################################
####                                                    ####
#### Biopython Entrez Taxonomy search automation script ####
####                                                    ####
############################################################

infile=open("C:\PostDoc\Geomicro_server\BioPythonScripts\Species_file.txt",'r')

from Bio import Entrez
Entrez.email = "urschm@rpi.edu"

#Open BLAST log file (tab-delimited text file)
TaxonomyOutFile = open("C:\PostDoc\Geomicro_server\BioPythonScripts\Taxonomy_out.txt","w")


for line in infile:
    handle = Entrez.esearch(db="Taxonomy", term=line)
    record = Entrez.read(handle)
    handle.close()

    handle2 = Entrez.efetch(db="Taxonomy", id=record["IdList"][0], retmode="xml")
    records = Entrez.read(handle2)
    records[0]['OtherNames']['Synonym']
    handle2.close()

    #Write header to BLAST log file
    #TaxonomyOutFile.write(records[0]['OtherNames']['Synonym'])
    print records[0]['OtherNames']['Synonym']

            

TaxonomyOutFile.close()
infile.close()
