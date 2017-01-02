###########################################
####                                   ####
#### Biopython BLAST automation script ####
####                                   ####
###########################################

from Bio.Blast import NCBIWWW
fasta_string = open("C:\PostDoc\Geomicro_server\BioPythonScripts\BLASTAutomation\FASTA_sequences_input.txt").read()
result_handle = NCBIWWW.qblast("blastn", "rRNA_typestrains/prokaryotic_16S_ribosomal_RNA", fasta_string,hitlist_size=1000)

from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(result_handle)

from Bio import Entrez
Entrez.email = "urschm@rpi.edu"

#Open BLAST log file (tab-delimited text file)
BLASTLogFile = open("C:\PostDoc\Geomicro_server\BioPythonScripts\BLASTAutomation\BLAST_Log_File.txt","w")

#Build string for BLAST log file header (column names)
BLASTHeader = "Database\t"
BLASTHeader = BLASTHeader+"Sequence\t"
BLASTHeader = BLASTHeader+"Taxonomy\t"
BLASTHeader = BLASTHeader+"Accession number\t"
BLASTHeader = BLASTHeader+"ID number\t"
BLASTHeader = BLASTHeader+"Alignment Length\t"
BLASTHeader = BLASTHeader+"Query sequence length\t"
BLASTHeader = BLASTHeader+"AlignQueryLength\t"
BLASTHeader = BLASTHeader+"Query start\t"
BLASTHeader = BLASTHeader+"Query end\t"
BLASTHeader = BLASTHeader+"Subject start\t"
BLASTHeader = BLASTHeader+"Subject end\t"
BLASTHeader = BLASTHeader+"Score\t"
BLASTHeader = BLASTHeader+"Bits\t"
BLASTHeader = BLASTHeader+"E-value\t"
BLASTHeader = BLASTHeader+"Coverage\t"
BLASTHeader = BLASTHeader+"Num_alignments\t"
BLASTHeader = BLASTHeader+"Identities\t"
BLASTHeader = BLASTHeader+"% Identities\t"
BLASTHeader = BLASTHeader+"Positives\t"
BLASTHeader = BLASTHeader+"% Positives\t"
BLASTHeader = BLASTHeader+"Gaps\t"
BLASTHeader = BLASTHeader+"% Gaps\t"
BLASTHeader = BLASTHeader+"Strand\t"
BLASTHeader = BLASTHeader+"Frame"

#Write header to BLAST log file
BLASTLogFile.write(BLASTHeader+"\n")

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            #Write BLAST results to log file

            #Calculate query length for alignment (including gaps)
            AlignQueryLength = float(hsp.query_end) - float(hsp.query_start)
            AlignQueryLength = AlignQueryLength + float(hsp.gaps)
            
            #Calculate percent identity
            IDPercent = float(hsp.identities)/float(AlignQueryLength)
            IDPercent = IDPercent * 100

            #Calculate percent positives
            PositivePercent = float(hsp.positives)/float(AlignQueryLength)
            PositivePercent = PositivePercent * 100

            #Calculate percent gaps
            GapsPercent = float(hsp.gaps)/float(AlignQueryLength)
            GapsPercent = GapsPercent * 100

            #Calculate coverage
            CoveragePercent = (float(hsp.query_end) - float(hsp.query_start))/float(blast_record.query_letters)
            CoveragePercent = CoveragePercent * 100

            #Get accession number of subject sequence
            Accession = alignment.title.split("|",4)[3]

            #Get id number of subject sequence
            IDNum = alignment.title.split("|",4)[1]

            #gb_ID_string = gb_ID_string+IDNum+","

            #Build next row for BLAST log file             
            BLASTLogString = str(blast_record.database)+"\t"
            BLASTLogString = BLASTLogString+str(alignment.title)+"\t"

            #Get taxonomy for sequence
            TaxHandle = Entrez.efetch(db="nucleotide", id=IDNum, retmode="xml")
            TaxRecord = Entrez.read(TaxHandle)
            TaxHandle.close()
            BLASTLogString = BLASTLogString+TaxRecord[0]["GBSeq_taxonomy"]+"\t"
            
            BLASTLogString = BLASTLogString+str(Accession)+"\t"
            BLASTLogString = BLASTLogString+str(IDNum)+"\t"
            BLASTLogString = BLASTLogString+str(alignment.length)+"\t"
            BLASTLogString = BLASTLogString+str(blast_record.query_letters)+"\t"
            BLASTLogString = BLASTLogString+str(float(AlignQueryLength))+"\t"
            BLASTLogString = BLASTLogString+str(hsp.query_start)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.query_end)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.sbjct_start)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.sbjct_end)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.score)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.bits)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.expect)+"\t"
            BLASTLogString = BLASTLogString+str(float(CoveragePercent))+"\t"
            BLASTLogString = BLASTLogString+str(hsp.num_alignments)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.identities)+"\t"
            BLASTLogString = BLASTLogString+str(float(IDPercent))+"\t"
            BLASTLogString = BLASTLogString+str(hsp.positives)+"\t"
            BLASTLogString = BLASTLogString+str(PositivePercent)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.gaps)+"\t"
            BLASTLogString = BLASTLogString+str(GapsPercent)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.strand)+"\t"
            BLASTLogString = BLASTLogString+str(hsp.frame)
            
            #Write next row to BLAST log file
            BLASTLogFile.write(BLASTLogString+"\n")
            

BLASTLogFile.close()

