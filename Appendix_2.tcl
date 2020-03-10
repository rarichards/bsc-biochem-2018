#APPENDIX 2

#The following script was used to extract bacterial genome sequences from embl files. 
#In a version of this script, the output file names were changed (indicated where appropriate) for use with archaeal genomes.


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




#First, open the folder containing the embl genome files

set flist [glob -nocomplain -directory filtered_genomes *.embl]

set num_genomes [llength $flist]



#Set up a genome counter
	
set cnt 1



#Make output file

set out_gnlst [open bac_genomes.csv w]

#In the archaeal version of this file the line read: 	set out_gnlst [open arc_genomes.csv w]

puts $out_gnlst "Accn,Genus,FullName,GenGC,cdsGC3,fveGC3,thrGC3,srnaGC"



#Open foreach loop - the following to be carried out on all genomes

foreach genome $flist {

	regexp {genomes/(.*?)\.embl} $genome all accn

	puts "Analysing $cnt of $num_genomes : $accn"


	
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



	#Find the species' name
	
	set posx [string first XX\nOS $contents]
	set posos [string first OS $contents $posx]
	set posn [string first \n $contents $posos]

	set speciesline [string range $contents [expr $posos + 2] $posn]
	set speciesline [string trim $speciesline]
	set parts [split $speciesline]
	set genus [lindex $parts 0]


	unset raw_seq
	unset contents



	#Split the annotation into genes using "FT   " which preceeds each gene
	
	regsub -all {FT   CDS             } $annot "£FT   CDS             " annot
	regsub -all {FT   tRNA            } $annot "£FT   tRNA            " annot
	regsub -all {FT   rRNA            } $annot "£FT   rRNA            " annot
	
	set gene_list [split $annot £]

	set gene_list [lrange $gene_list 1 end]
	set num_genes [llength $gene_list]



	#Set/Reset gene counters

	set gene_cnt 1
	set num_CDS 0
	set num_RNA 0

	
	
	#Set/Reset parameters
	
	set total_gc3 0.0
	set total5_gc3 0.0
	set total3_gc3 0.0
	set totalrna_gc 0.0



	#Calculate genomic GC content
	
	set gnme_l [string length $clean_seq]
	
	set gnme_gc [regexp -all {[c|g]} $clean_seq]
	set gnme_gc [expr $gnme_gc * 1.0] 
	set gnme_gc [expr {$gnme_gc/$gnme_l} * 100]



	#Open foreach loop - the following to be carried out on all genes
	
	foreach gn $gene_list {


		#Assign gene type

		if {[regexp {FT   CDS             } $gn] ==1 } {
			set gn_type "CDS"
			incr num_CDS
	
		} elseif {[regexp {FT   rRNA            } $gn] ==1 } {
			set gn_type "rRNA"
			incr num_RNA
	
		} elseif {[regexp {FT   tRNA            } $gn] ==1 } {
			set gn_type "tRNA"
			incr num_RNA
	
		} else {
			set gn_type "error"
			puts "Gene Type Error"
		}
		
	

		#Identify gene position
	
		regexp {[^0-9]([0-9]+?)\.\.([0-9]+?)[^0-9]} $gn all strt nd



		#Identify the locus tag - If a gene does not have a locus tag,use the gene count

		if {[regexp {/locus_tag="(.*?)"} $gn] } {
		regexp {/locus_tag="(.*?)"} $gn all loctag
		} else {set loctag $gene_cnt}



		#Set the gene name

		if {[regexp {/gene="(.*?)"} $gn]} {
		regexp {/gene="(.*?)"} $gn all genename
		} else {set genename $gene_cnt}
	
	
	
		#Extract the gene nucleotide sequence and clean up

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
	

		
		#The following only needs to be carried out on coding DNA sequences, not RNAs
		
		if {$gn_type == "CDS"} {



			#Quality control tests
			
			set fail 0
			while {$fail ==0} {
		

				#1. Divisible by three (whole number of codons)
					if {[expr $seql % 3] >0} {set fail 1}
				

				#2 Only includes letters a, t, g or c
					if {[regexp {[^a|c|t|g]} $seq]} {set fail 2}
					

				#3. Ends in a stop codon
					set last_codon [string range $seq end-2 end]
					if {[regexp {[^tga|taa|tag]} $last_codon]} {set fail 3}
				
				#4. Does not contains internal stop codons
					for {set i 0} {$i < [expr $seql - 5]} {incr i +3} {
					set cdn [string range $seq $i [expr $i +2]]
					if {[regexp {tga|taa|tag} $cdn]} {set fail 4}}		
			
				#5. Starts with a start codon - almost all bacterial genomes start codons end 'tg'
					set first_codon [string range $seq 1 2]
					if {![regexp {tg} $first_codon]} {set fail 5}
				
				
				#6. Gene length is >=74nt, to avoid 5' 3' overlap
					if {$seql <75} {set fail 6}
				
					if {$fail >0} {
					puts "$loctag : fail"
					
					}

				#Reassign fail code if it passes
					if {$fail ==0} {set fail -1}
				
			}



			#The following is carried out only on genes which pass quality control

			if {$fail == -1} {


		
			#Calculate the GC3 content of whole gene.
	
			set cdns [expr $seql / 3]
			set gc3 0
	
			for {set i 2} {$i <= $seql} {incr i +3} {
				set pos3 [string range $seq $i $i]
				if {$pos3 == "c"} {set gc3 [expr $gc3 + 1.0]}
				if {$pos3 == "g"} {set gc3 [expr $gc3 + 1.0]}
			}

			set gc3_pcnt [expr {$gc3/$cdns} * 100]
			set total_gc3 [expr $total_gc3 + $gc3_pcnt]
	
		
		
						
		
			#Extract 5 prime sequence
			#If it is the first gene in the genome, may need to find preceding bases from the end of the genome
		
			set fiveprm [string range $clean_seq [expr $strt -5] [expr $strt +35]]
			
			if {$strt == 1} {
				set pre_seq [string range $clean_seq end-3 end]
				set fiveprm "$pre_seq$fiveprm"
			}
		
			if {$strt == 2} {
				set pre_seq [string range $clean_seq end-2 end]
				set fiveprm "$pre_seq$fiveprm"
			}
		
			if {$strt == 3} {
				set pre_seq [string range $clean_seq end-1 end]
				set fiveprm "$pre_seq$fiveprm"
			}
		
			if {$strt == 4} {
				set pre_seq [string range $clean_seq end end]
				set fiveprm "$pre_seq$fiveprm"
			}
	
	
			#Extract 3 prime sequence
			#If it is the last gene in the genome, may need to find last bases from the start of the genome
		
			set threeprm [string range $clean_seq [expr $nd -36] [expr $nd +4]]
		
			set clean_seql [string length $clean_seq]
			set finalpos [expr $clean_seql - $nd]
		
			if {$finalpos == 0} {
				set postseq [string range $clean_seq 0 3]
				set threeprm "$threeprm$postseq"
			}
		
			if {$finalpos == 1} {
				set postseq [string range $clean_seq 0 2]
				set threeprm "$threeprm$postseq"
			}
		
			if {$finalpos == 2} {
				set postseq [string range $clean_seq 0 1]
				set threeprm "$threeprm$postseq"
			}
		
			if {$finalpos == 3} {
				set postseq [string range $clean_seq 0 0]
				set threeprm "$threeprm$postseq"
			}
	
	
			set fiveprm [string tolower $fiveprm] 
			set threeprm [string tolower $threeprm]	
			
		
	
			#Account for complementarity in 5' and 3' sequences, using complement procedure

			if {$comp ==1} {
		
				Complement $fiveprm
				set threeprmA $revs
		
				Complement $threeprm
				set fiveprm $revs
		
				set threeprm $threeprmA
			
			} 
			
			
			
			# Calculate 5' GC3 content.

			set fivel [string length $fiveprm]
			set fivecdns [expr $fivel / 3]
			set fivegc3 0
		
			for {set i 0} {$i <= $fivel} {incr i +3} {
				set pos3 [string range $fiveprm $i $i]
				if {$pos3 == "c"} {set fivegc3 [expr $fivegc3 + 1.0]}
				if {$pos3 == "g"} {set fivegc3 [expr $fivegc3 + 1.0]}
			}
			
			set fivegc3_pcnt [expr {$fivegc3/$fivecdns} * 100]
			set total5_gc3 [expr $total5_gc3 + $fivegc3_pcnt]	
		


			#Calculate 3' GC3 content.
	
			set threel [string length $threeprm]
			set threecdns [expr $threel / 3]
			set threegc3 0
		
			for {set i 0} {$i <= $threel} {incr i +3} {
				set pos3 [string range $threeprm $i $i]
				if {$pos3 == "c"} {set threegc3 [expr $threegc3 + 1.0]}
				if {$pos3 == "g"} {set threegc3 [expr $threegc3 + 1.0]}
			}
			
			set threegc3_pcnt [expr {$threegc3/$threecdns} * 100]
			set total3_gc3 [expr $total3_gc3 + $threegc3_pcnt]	
			
			}
	
		
		#End of actions on CDS
		} else {
	
	
	
			#Calculate RNA GC content.
	
			set gc [regexp -all {[c|g]} $seq]
			set gc [expr $gc * 1.0]
			set rna_gc [expr {$gc/$seql} * 100]
		
			set totalrna_gc [expr $totalrna_gc + $rna_gc]
			
		}
	
	
	#Increase gene counter
	incr gene_cnt	

	
	#End of gene loop
	}
	
	
	
	#Calculate genome means
	
	set mean_gc3 [expr $total_gc3 / $num_CDS]
	set mean5_gc3 [expr $total5_gc3 / $num_CDS]
	set mean3_gc3 [expr $total3_gc3 / $num_CDS]
	if {$totalrna_gc >0} {set meanrna_gc [expr $totalrna_gc / $num_RNA]} else {set meanrna_gc "NaN"}

	set d53 [expr $mean5_gc3 - $mean3_gc3]



	#Output parameters to file
	if {$meanrna_gc != "NaN"} {
		puts $out_gnlst "$accn,$genus,$speciesline,$gnme_gc,$mean_gc3,$mean5_gc3,$mean3_gc3,$meanrna_gc"
 	}
	
	
	#Increase genome counter
	incr cnt
	
	
	
#End of genome loop
}


#Close the output file
close $out_gnlst
