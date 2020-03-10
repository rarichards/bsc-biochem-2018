#The following script was used to extract ribosomal genes and write them into one FASTA file per genome


#Procedures must be placed at the beginning of the script, but they are not called until later
#Complement procedure to reverse sequence and exchange bases

	proc Complement s {

	global revs

	regsub -all {a} $s {1} s
	regsub -all {t} $s {a} s
	regsub -all {1} $s {t} s

	regsub -all {c} $s {2} s
	regsub -all {g} $s {c} s
	regsub -all {2} $s {g} s

	set revs [string reverse $s]

	return}



#start

set flist [glob -nocomplain -directory filtered_genomes *.embl]


#set up a counter

	set num_genomes [llength $flist]
	set cnt 1



foreach genome $flist {

	regexp {genomes/(.*?)\.embl} $genome all accn

	puts "Analysing $cnt of $num_genomes : $accn"


	#Don't add to an existing file

	if {[file exists "$accn_rib.fasta"] ==1} {file delete "$accn_rib.fasta"} 

	set out_genes [open "$accn_rib.fasta" a]
	
	
	
	#Identify where to look for sequence and extract it
	
	set in [open $genome r]
	set contents [read $in]
	close $in

	set posxx [string first "XX\nSQ" $contents]
	set pos_other [string first "other;" $contents $posxx]
	set posnl [string first \n $contents $pos_other]

	set raw_seq [string range $contents $posnl end]


	#Clean up and remove characters which are not letters
	
	regsub -all {[^A-Za-z]} $raw_seq {} clean_seq


	#Label everything else in the file as 'annotation'
	
	set annot [string range $contents 0 $posxx] 



	unset raw_seq
	unset contents



	#Split the annotation into genes using "FT   " which preceeds each gene
	
	regsub -all {FT   CDS             } $annot "£FT   CDS             " annot
	
	set gene_list [split $annot £]

	set gene_list [lrange $gene_list 1 end]
	

	set num_genes [llength $gene_list]

	set gene_cnt 1
	
	
	set num_CDS 0
	set num_RNA 0


	

	foreach gn $gene_list {


		#Identify gene position
	
		regexp {[^0-9]([0-9]+?)\.\.([0-9]+?)[^0-9]} $gn all strt nd

		
		#Identify gene product
		
		if {[regexp {/product="(.*?)"} $gn]} {
		regexp {/product="(.*?)"} $gn all product
		} 
		
		if {[regexp {ribosome|ribosomal} $product]} {
		} {
		
		#only do the following if the gene description includes "ribosome" or "ribosomal"
		
		
		#Identify the locus tag - If a gene does not have a locus tag,use the gene count

		if {[regexp {/locus_tag="(.*?)"} $gn] } {
		regexp {/locus_tag="(.*?)"} $gn all loctag
		} else {set loctag $gene_cnt}



		#Extract nucleotide sequence and clean up

		set seq [string range $clean_seq [expr $strt -1] [expr $nd -1]]
	
		set seq [string tolower $seq]
		set seql [string length $seq]



		#Account for complementary strands using complement procedure
	
		regexp {(.*?)\n} $gn all firstline
	
		if {[regexp {complement} $firstline]} {
			set comp 1
		} else {set comp 0}
	
	
		if {$comp ==1} {
			Complement $seq
			set seq $revs		
		} 
	


		#Quality control tests
		
		set fail 0
			
		while {$fail ==0} {
		

			#1. divisible by three
			
				if {[expr $seql % 3] >0} {set fail 1}
				

			#2 includes letters other than atgc
				if {[regexp {[^a|c|t|g]} $seq]} {set fail 2}
					

			#3. ends in codon thats not a stop codon
				set last_codon [string range $seq end-2 end]
				if {[regexp {[^tga|taa|tag]} $last_codon]} {set fail 3}
				
			#4. contains internal stop codon
				for {set i 0} {$i < [expr $seql - 5]} {incr i +3} {
				set cdn [string range $seq $i [expr $i +2]]
				if {[regexp {tga|taa|tag} $cdn]} {set fail 4}}		
			
			#5. starts with a start codon - almost all bacterial genomes start codons end tg 
				set first_codon [string range $seq 1 2]
				if {![regexp {tg} $first_codon]} {set fail 5}
				
				
			#6. length >=74nt to give meaningful data by avoiding 5' 3' overlap
				if {$seql <75} {set fail 6}
				
				if {$fail >0} {
				puts "$loctag : fail"
					
				}

			#if it passes the 4 fail tests then it must be reassigned -1
				if {$fail ==0} {set fail -1}
				
			}

			

			if {$fail == -1} {
			
			#Output to file if passed quality control
			puts $out_genes ">$loctag\n$seq"
			
			}

	}	 
	#Increase gene counter
	incr gene_cnt	

	}
	
	#Close output file
	close $out_genes
	
	#Increase genome counter
	incr cnt

}


