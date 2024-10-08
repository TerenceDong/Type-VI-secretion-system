gunzip E.coli_genome.fasta.gz
grep '>' E.coli_genome.fasta | wc -l # 3780 总数

# select the DNA sequence containing prophage BP-933W -------------------------------------------
makeblastdb -in E.coli_genome.fasta -dbtype nucl -out E.coli_genome_db
blastn -query BP-933W-changed.fasta -db E.coli_genome_db -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -evalue 1e-10 -out results1.tsv
awk -F '\t' '$14>30 {print $2}' results1.tsv | sort | uniq > matched_headers1.txt
cat matched_headers1.txt | wc -l # 492
grep -A 1 -F -f matched_headers1.txt E.coli_genome.fasta > E.coli_genome_filter1.fasta

# predict the genes in E.coli_genome_filter1.fasta and translate it into AA -----------------------------------
mkdir -p Ecoli_split
num_records=$(grep -c "^>" E.coli_genome_filter1.fasta)
echo $num_records
k=123 
records_per_part=$(( (num_records + k - 1) / k )) 
awk -v k=$k -v rpp=$records_per_part '
  BEGIN { 
    file_index = 1; 
    file_name = "A" file_index ".fasta"; 
  }
  /^>/ {
    if (count == rpp) { 
      file_index++; 
      file_name = "A" file_index ".fasta"; 
      count = 0; 
    } 
    count++;
  }
  {
    print > file_name;
  }
' E.coli_genome_filter1.fasta
mv A*.fasta Ecoli_split
cd Ecoli_split
ls *.fasta | sed 's/\.fasta$//g' > Ecoli_list.txt

mkdir -p prodigal
parallel -j 50 --xapply \
"prodigal -i {1}.fasta \
-d prodigal/{1}.fa \
-o prodigal/{1}.gff \
-p meta -f gff" \
::: `tail -n+1 Ecoli_list.txt`

cat prodigal/*.fa > prodigal/A_all.fa

transeq -sequence prodigal/A_all.fa \
-outseq prodigal/A_all_protein.fa -trim Y 

# determine the protein sequence containing th first 10bp of EHEC_tssM-changed.fasta (1154bp) --------------------------------------

makeblastdb -in A_all_protein.fa -dbtype prot -out A_all_protein
blastp -query EHEC_tssM_ancestor-changed.fasta -db A_all_protein -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -evalue 1e-10 -out blast_results_1.tsv
awk -F '\t' '$8==1 && $14>30 && $3>70 {print $2}' blast_results_1.tsv | sed 's/_.*$//g' | sort | uniq > matched_headers2.txt 
awk -F '\t' '$8==1 && ($14<=30 || $3<=70) {print $2}' blast_results_1.tsv | sed 's/_.*$//g' | sort | uniq > matched_headers3.txt 
awk -F '\t' '$8==1 {print $2}' blast_results_1.tsv | sed 's/_.*$//g' | sort | uniq > matched_headers4.txt 

comm -23 matched_headers1.txt matched_headers4.txt > Ecoli_notcontaining_pre.txt 
cat Ecoli_notcontaining_pre.txt matched_headers3.txt > Ecoli_notcontaining.txt 
# manually add 13 strains
# CP104645.1
# CP104647.1
# CP104649.1
# CP120944.1
# CP124740.1
# CP124819.1
# CP126904.1
# CP126906.1
# CP129257.1
# CP129269.1
# CP129345.1
# CP132297.1
# CP132298.1

grep -A 1 -F -f Ecoli_notcontaining.txt E.coli_genome.fasta > Ecoli_notcontaining.fa

# identify the genes of Ecoli_notcontaining.fa --------------------------------------------------------------

mkdir -p prodigal_2
num_records=$(grep -c "^>" Ecoli_notcontaining.fa)
echo $num_records
k=179  # Set the number of parts you want to split into
records_per_part=$(( (num_records + k - 1) / k ))  
mv B*.fasta prodigal_2
cd prodigal_2
ls *.fasta | sed 's/\.fasta$//g' > Ecoli_list.txt

parallel -j 50 --xapply \
"prodigal -i {1}.fasta \
-d {1}.fa \
-o {1}.gff \
-p meta -f gff" \
::: `tail -n+1 Ecoli_list.txt`

cat B*.fa > Ecoli_notcontaining_prodigal.fa
transeq -sequence Ecoli_notcontaining_prodigal.fa \
-outseq Ecoli_notcontaining_all_protein.fa -trim Y 

# determine the species list containing EHEC_tssM sequence --------------------------------------------------------------

makeblastdb -in Ecoli_notcontaining_all_protein.fa -dbtype prot -out Ecoli_notcontaining_all_protein
blastp -query EHEC_tssM.fasta -db Ecoli_notcontaining_all_protein -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -evalue 1e-10 -out blast_results_2.tsv
awk -F '\t' '$8==1 && $14>30 && $3>70 {print $2}' blast_results_2.tsv | sed 's/_.*$//g' | sort | uniq > matched_headers5.txt 
grep -F -f matched_headers5.txt E.coli_genome.fasta > matched_headers5_fullname.txt # N=295
# identify the species list without EHEC_tssM sequence 
awk -F '\t' '$8==1 {print $2}' blast_results_2.tsv | sed 's/_.*$//g' | sort | uniq > matched_headers6.txt 
sort Ecoli_notcontaining.txt -o Ecoli_notcontaining.txt
comm -23 Ecoli_notcontaining.txt matched_headers6.txt > matched_headers7.txt 
grep -F -f matched_headers7.txt E.coli_genome.fasta > matched_headers7_fullname.txt

# extract the strain names, ID and aligned sequences (N=169, 33, 290) ---------------------------

# N=169--------------------------------------------------------------------------------------------

awk -F '\t' '$8==1' ../blast_results_1.tsv > blast_N169.tsv

awk '/^>/ {if (seq) print seq; print; seq=""} !/^>/ {seq=seq""$0} END {if (seq) print seq}' Ecoli_split/prodigal/A_all_protein.fa \
> A_all_protein_changed.fa

bash sequence_extraction.sh A_all_protein_changed.fa blast_N169.tsv N169_sequence.tsv
awk -F'\t' '{ $1 = gensub(/_.*/, "", "g", $1); print $0 }' OFS='\t' N169_sequence.tsv > N169_sequence_tmp.tsv

grep -F -f ../matched_headers4.txt ~/TRAINING/sun/references/E.coli_genome.fasta > N169_strainname.tsv
# To split the content of a TSV file into two columns based on the first blank space
awk '{first=$1; $1=""; gsub(/ /, "_"); print first, substr($0, 2)}' N169_strainname.tsv > N169_strainname_tmp.tsv
sed -i "s/>//g" N169_strainname_tmp.tsv

# merge two tsv files
sort -k1,1 N169_sequence_tmp.tsv > N169_sequence_sorted.tsv
sort -k1,1 N169_strainname_tmp.tsv > N169_strainname_tmp_sorted.tsv
awk 'NR==FNR{a[$1]=$2; next} {print $1, a[$1], $2}' N169_strainname_tmp_sorted.tsv N169_sequence_sorted.tsv  | tr ' ' '\t' > N169_final.tsv


# N=323 prior to subsequent steps--------------------------------------------------------------------------------------------

awk -F '\t' '$14>30 {print $2}' ../blast_results_1.tsv | sed 's/_.*$//g' | sort | uniq > N459.txt 
comm -23 N459.txt ../matched_headers4.txt > N290.txt

comm -23 matched_headers1.txt N459.txt > N33.txt


# N=290--------------------------------------------------------------------------------------------

grep -F -f N290.txt ../blast_results_1.tsv > blast_N290.tsv # N=294 including replicate strains

# awk -F '\t' '{print $2}' blast_N290.tsv | sed 's/_.*$//g' | sort | uniq > no.txt # N=290
# rm -f no.txt

bash sequence_extraction.sh A_all_protein_changed.fa blast_N290.tsv N290_sequence.tsv
awk -F'\t' '{ $1 = gensub(/_.*/, "", "g", $1); print $0 }' OFS='\t' N290_sequence.tsv > N290_sequence_tmp.tsv

grep -F -f N290.txt ~/TRAINING/sun/references/E.coli_genome.fasta > N290_strainname.tsv
# To split the content of a TSV file into two columns based on the first blank space
awk '{first=$1; $1=""; gsub(/ /, "_"); print first, substr($0, 2)}' N290_strainname.tsv > N290_strainname_tmp.tsv
sed -i "s/>//g" N290_strainname_tmp.tsv

# merge two tsv files
sort -k1,1 N290_sequence_tmp.tsv > N290_sequence_sorted.tsv
sort -k1,1 N290_strainname_tmp.tsv > N290_strainname_tmp_sorted.tsv
awk 'NR==FNR{a[$1]=$2; next} {print $1, a[$1], $2}' N290_strainname_tmp_sorted.tsv N290_sequence_sorted.tsv  | tr ' ' '\t' > N290_final.tsv

# N=33--------------------------------------------------------------------------------------------

grep -F -f N33.txt ~/TRAINING/sun/references/E.coli_genome.fasta > N33_strainname.tsv
# To split the content of a TSV file into two columns based on the first blank space
awk '{first=$1; $1=""; gsub(/ /, "_"); print first, substr($0, 2)}' N33_strainname.tsv > N33_strainname_final.tsv
sed -i "s/>//g" N33_strainname_final.tsv























