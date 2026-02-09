#!/bin/bash

module load bcftools minimap2
module load StdEnv/2020 jellyfish/2.3.0
# This file is to keep common snippets used but are not included in the pipeline
###########################################
# Work on NCLCN-261 assembly
# Work with RepeatMasker file

  module load repeatmasker
  grep -P '^[[:space:]]*[0-9]+' NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.haplotype1.fa.out \
  | grep -E 'LINE/L1|SINE/Alu|SVA' \
  > NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.haplotype1.youngTEs.out

  awk 'BEGIN{OFS="\t"}
  {
    # convert to 0-based BED start:  start = $6 - 1
    # contig = $5, end = $7
    # then tack on all original columns as annotation
    printf("%s\t%d\t%s", $5, $6-1, $7)
    for(i=1;i<=NF;i++) printf("\t%s", $i)
    print ""
  }' NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.haplotype2.youngTEs.out \
  > to_lift.hap2.bed

  k8 /home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/tools/paftools/paftools.js liftover asm2chm13.paf to_lift.hap2.bed > hap22chm13.fa.out

  awk 'BEGIN{FS=OFS="\t"}
    # First file: unified liftover results
    NR==FNR {
      key = $4                    # e.g. "haplotype1-0000001_0_2340"
      chrMap[key]   = $1          # new chr
      startMap[key] = $2          # new 0-based start
      endMap[key]   = $3          # new end
      next
    }
    # Second file: original to_lift.bed
    {
      key = $1 "_" $2 "_" $3      # build the same key from cols 1–3
      if(key in chrMap) {
        $1 = chrMap[key]
        $2 = startMap[key]
        $3 = endMap[key]
      } else {
        # optional: warn if you have a bed entry with no mapping
        print "WARNING: no mapping for " key > "/dev/stderr"
      }
      print
    }
  ' hap12chm13.fa.out to_lift.hap1.bed > hap1.chm13.bed




  awk -F'\t' ' 
  # First pass: count up how many times each ID appears
  NR==FNR { cnt[$4]++; next }  
  # Second pass: on the same file again, print only if count>1
  cnt[$4]>1
  ' hap12chm13.fa.out hap12chm13.fa.out

# After running Genpipes, rename scaffolds to chr based on chromAlias.txt
  result_dir="$wd/results/iPSC_K27"
  markname="iPSC_K27"
  module load mugqic/genpipes/4.4.5
  module load mugqic/python/3.10.4
  mkdir -p $result_dir
    cd $result_dir || exit
    directory=linear_chm13_w_trimmomatic
    mkdir -p $directory
    chipseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
                    $MUGQIC_PIPELINES_HOME/pipelines/common_ini/narval.ini t2t.ini \
                    -r ${markname}.readset.tsv -d ${markname}.design.txt \
                    -o $directory -f > "chipseqScript_${markname}.txt" 
    bash "chipseqScript_${markname}.txt" 
    cd $wd || exit
  awk '
    BEGIN{ FS=OFS="\t" }
    # read chromAlias.txt into map: key=old name, val=new UCSC name
    NR==FNR {
      map[$1] = $2
      next
    }
    # then for each line of your BED:
    {
      # if there is a mapping for $1, replace it; otherwise leave it
      $1 = ($1 in map ? map[$1] : $1)
      print
    }
  ' Genome/NCLCN-261/NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.chromAlias.txt   results/iPSC_K27/linear_NCLCN-261/peak_call/iPSC_K27/iPSC_K27/iPSC_K27.iPSC_K27_peaks.narrowPeak.bed > results/iPSC_K27/linear_NCLCN-261/peak_call/iPSC_K27/iPSC_K27/iPSC_K27.iPSC_K27_peaks.narrowPeak.ucsc.bed
# liftOver
  module load minimap2
  cd /home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/Genome/NCLCN-261 || exit
  CHM13=/home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/Graph_genome_data/chm13v2.0.fa
  peakresult=results/iPSC_K27/linear_NCLCN-261_hap2/peak_call/iPSC_K27/iPSC_K27/
  out=$wd/$peakresult/iPSC_K27.hap12chm13.bed
  ASM=$wd/Genome/NCLCN-261/NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.haplotype2.fa
  # minimap2 -t64 -x asm20 -c $CHM13 $ASM > $wd/Genome/NCLCN-261/hap22chm13.paf
  bed=$wd/$peakresult/iPSC_K27.iPSC_K27_peaks.narrowPeak.bed
  k8 $wd/tools/paftools/paftools.js liftover \
    $wd/Genome/NCLCN-261/hap12chm13.paf $bed \
    > $out

   awk 'BEGIN{FS=OFS="\t"}
    # First file: unified liftover results
    NR==FNR {
      key = $4                    # e.g. "haplotype1-0000001_0_2340"
      chrMap[key]   = $1          # new chr
      startMap[key] = $2          # new 0-based start
      endMap[key]   = $3          # new end
      next
    }
    # Second file: original to_lift.bed
    {
      key = $1 "_" $2 "_" $3      # build the same key from cols 1–3
      if(key in chrMap) {
        $1 = chrMap[key]
        $2 = startMap[key]
        $3 = endMap[key]
      } else {
        # optional: warn if you have a bed entry with no mapping
        print "WARNING: no mapping for " key > "/dev/stderr"
      }
      print
    }
  ' Genome/NCLCN-261/hap12chm13.fa.out Genome/NCLCN-261/to_lift.hap1.bed > Genome/NCLCN-261/hap1.chm13.bed

  awk 'BEGIN{FS=OFS="\t"}
     # select lines where the contig field contains “haplotype” (novel insertions)
     $1 ~ /haplotype/ {
       # print chr, start, end, subfamily (col 13 in your example)
       print $1, $2, $3, $13
     }
  ' Genome/NCLCN-261/hap1.chm13.bed > Genome/NCLCN-261/hap1.novelinsert.bed

# 1) Count by family (Alu, L1, SVA), sorted descending
  awk 'BEGIN{FS=OFS="\t"}
     {
       if($4 ~ /^Alu/) fam="Alu"
       else if($4 ~ /^L1/) fam="L1"
       else if($4 ~ /^SVA/) fam="SVA"
       else next
       famCount[fam]++
     }
     END {
       for(f in famCount) print f, famCount[f]
     }' Genome/NCLCN-261/hap1.novelinsert.bed \
  | sort -k2,2nr > Genome/NCLCN-261/hap1.family_counts.tsv

# 2) Count by subfamily, sorted descending
  {
  # Print header
  printf "subfamily\tcount\n"
  # Count by subfamily for novel insertions, prefix with family order for sorting
  awk 'BEGIN{FS=OFS="\t"}
       $1 ~ /haplotype/ { cnt[$4]++ }
       END {
         for(sf in cnt) {
           if   (sf ~ /^Alu/) fam=1
           else if(sf ~ /^L1/)  fam=2
           else if(sf ~ /^SVA/) fam=3
           else                   fam=4
           print sf, cnt[sf], fam
         }
       }' Genome/NCLCN-261/hap1.novelinsert.bed |
  # Sort by family (Alu→L1→SVA) then count descending
  sort -k3,3n -k2,2nr |
  # Drop the family‐order column, keeping only subfamily and count
  awk 'BEGIN{FS=OFS="\t"} { print $1, $2 }'
  } > Genome/NCLCN-261/hap1.sub_family_counts.tsv

  grep -v '^haplotype' Genome/NCLCN-261/hap1.chm13.bed \
    | awk -F'\t' '$19 != "*"' \
    | bedtools intersect -v -b - -a Genome/t2t_L1_Alu_SVA.bed \
    > Genome/NCLCN-261/hap1.excluded.chm13.bed


  grep -v '^haplotype' Genome/NCLCN-261/hap1.chm13.bed \
    | bedtools intersect -v -a - -b Genome/t2t_L1_Alu_SVA.bed \
    > Genome/NCLCN-261/hap1.excluded.chm13.bed

  awk 'BEGIN{
  while((getline<"Genome/NCLCN-261/chromAlias.from_mashmap.txt")>0){
    alias[$1]=$2
  }
  }
  {
  # replace query and/or target names if present in alias
  if($1 in alias) $1=alias[$1]
  if($6 in alias) $6=alias[$6]
  print
  }' OFS="\t" Genome/NCLCN-261/asm2chm13.paf > Genome/NCLCN-261/chr_asm2chm13.paf

  grep -v $'\ttp:A:S\b' Genome/NCLCN-261/chr_asm2chm13.paf > Genome/NCLCN-261/chr_asm2chm13_primsup.paf

  gunzip -c Genome/NCLCN-261/NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.dip.vcf.gz > Genome/NCLCN-261/NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.dip.vcf
  k8 /home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/tools/paftools/paftools.js liftover Genome/NCLCN-261/chr_asm2chm13_primsup.paf Genome/NCLCN-261/NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.dip.vcf \
      > Genome/NCLCN-261/NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.dip.chm13.vcf
  bcftools query -l NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.dip.vcf.gz 
  bcftools view -h NCLCN-261_verkko221_HiFi_ONTUL_POREC_round2.dip.vcf.gz  | grep '^##contig' 
  grep '^chr20\s' Genome/NCLCN-261/chr_asm2chm13_primsup.paf | awk '{ if($3 <= 12345 && $4 >= 12345) print }' | head
  cut -f1 chr_asm2chm13_primsup.paf | sort | uniq 

  
###########################################
  json=Pangenomes/vg_giraffe/vcfbub/graphs_chunks/chr1.json
  head -c 50 Pangenomes/vg_giraffe/vcfbub/graphs/_0_GRCh38#0#chr1#0_0_248177868.json
  grep -oP 'node.{0,50}' $json | head -n 1
  comm -23 Pangenomes/vg_giraffe/vcfbub/graphs/_0_GRCh38#0#chr1#0_0_248177868_ids.txt Pangenomes/vg_giraffe/vcfbub/graphs_chunks/chr1_ids.txt > Pangenomes/vg_giraffe/vcfbub/graphs_chunks/check_chr1_ids.txt

  json=Pangenomes/vg_giraffe/vcfbub/graphs/_0_GRCh38#0#chr1#0_0_248177868.json
  grep -oP '{"id": *"[^"]+", *"sequence": *"[^"]+"}' $json | \
      sed -E 's/.*"id": *"([^"]+)", *"sequence": *"([^"]+)".*/\1\t\2/' > Pangenomes/vg_giraffe/vcfbub/graphs/_0_GRCh38nodes_seq.tsv

  awk 'NR==FNR {ids[$1]; next} { if ($1 in ids) print }' \
    Pangenomes/vg_giraffe/vcfbub/graphs/unique_chr1_ids_GRCh38vschunk.txt \
    all_nodes.tsv > \
    Pangenomes/vg_giraffe/vcfbub/graphs/unique_chr1_ids_GRCh38vschunk.seq.txt

  awk -F'\t' '$1 ~ /^haplotype/ {count++} END{print count}' /home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/Genome/NCLCN-261/to_lift.bed
###########################################
# Look at unmapped reads
  module load bedtools 
  module load samtools
  samtools view -b -f 4 $wd/results/146/vg_giraffe_vcfbub/treatment_alignments.bam > $wd/results/146/vg_giraffe_vcfbub/unmapped_reads_graph.bam

  # samtools view -N unmapped_reads_graph.txt $wd/results/146/linear_chm13/alignment/146/ZNF146/146.ZNF146.sorted.dup.bam | awk '{print $1, $3, $4, $5, $6}' OFS='\t' > unmapped_in_graph_linear_summary.tsv
  samtools view -h -N unmapped_reads_graph.txt $wd/results/146/linear_chm13/alignment/146/ZNF146/146.ZNF146.sorted.dup.bam | awk 'BEGIN{OFS="\t"} /^@/ || ($5 == 60)' | samtools view -bS - > linear_mapq60_unmapped_in_graph.bam
  samtools sort -o linear_mapq60_unmapped_in_graph.sorted.bam linear_mapq60_unmapped_in_graph.bam
  samtools index linear_mapq60_unmapped_in_graph.sorted.bam
  bedtools bamtobed -i linear_mapq60_unmapped_in_graph.sorted.bam > reads.bed

  linear_bed=$wd/results/146/linear_chm13/peak_call/146/ZNF146/146.ZNF146_peaks.narrowPeak.bed
  graph_bed=$wd/results/146/vg_giraffe_vcfbub/callpeaks/all_linear_peaks.bed
# look at lost peaks that may be supported by unmapped reads in graph
  awk '$5 >= 50' $linear_bed > linear_50.bed
  bedtools intersect -v -a linear_50.bed -b $graph_bed > lost_peaks_linear.bed
  bedtools bamtobed -i linear_mapq60_unmapped_in_graph.sorted.bam > reads_unmapped_in_graph.bed
  bedtools intersect -u -a lost_peaks_linear.bed -b reads_unmapped_in_graph.bed > lost_peaks_supported_by_reads.bed

# Repeat content of unmapped reads
  BAM_ln="results/146/linear_chm13/alignment/146/ZNF146/146.ZNF146.primary.bam"
  BAM_CHM="results/146/vg_giraffe_chm13_hardhitcap2500_hitcap100_scorefraction0.8/treatment_alignments.bam"
  BAM_VCF="results/146/vg_giraffe_vcfbub_hardhitcap2500_hitcap100_scorefraction0.8/treatment_alignments.bam"
  OUTDIR="part_results/unmapped_reads_146_hardhitcap2500_hitcap100_scorefraction0.8/unmapped_compare"
  mkdir -p "$OUTDIR"
  # VCFBUB
  samtools view -@8 -f 4 -F 2304 "$BAM_ln" | awk '{print toupper($10)}' | sort -u > "$OUTDIR/ln.seq"
  samtools view -@8 -f 4 -F 2304 "$BAM_VCF" | awk '{print toupper($10)}' | sort -u > "$OUTDIR/vcf.seq"
  samtools view -@8 -f 4 -F 2304 "$BAM_CHM" | awk '{print toupper($10)}' | sort -u > "$OUTDIR/chm.seq"

samtools view -@8 -f 4 -F 2304 "$BAM_VCF" \
| awk 'BEGIN{OFS="\t"} {for(i=1;i<=11;i++){printf (i==1?$i:OFS $i)}; printf "\n"}' \
> $OUTDIR/vcf.names
  comm -23 "$OUTDIR/vcf.seq" "$OUTDIR/chm.seq" > "$OUTDIR/vcf.only.seq"

  awk 'BEGIN{OFS=""}
     {id++; printf(">vcf_only_%06d count=%s\n%s\n", id, ($1? $1:1), $2)}' \
     "$OUTDIR/vcf.only.seq" > "$OUTDIR/vcf.only.fa"

  mkdir -p "$OUTDIR/rm_out"
  RepeatMasker -pa 8 -species human -s -a -dir "$OUTDIR/rm_out" "$OUTDIR/vcf.only.unique.fa"

  # 2) Intersect & set-difference on names
  comm -12 "$OUTDIR/vcf.seq" "$OUTDIR/chm.seq" > "$OUTDIR/common.seq"
  comm -23 "$OUTDIR/vcf.seq" "$OUTDIR/chm.seq" > "$OUTDIR/vcf.only.seq"
  comm -13 "$OUTDIR/vcf.seq" "$OUTDIR/chm.seq" > "$OUTDIR/chm.only.seq"

  # # 3) Print counts
  # echo "VCFBUB unmapped: $(wc -l < "$OUTDIR/vcf.names")"
  # echo "CHM13  unmapped: $(wc -l < "$OUTDIR/chm.names")"
  # echo "Common unmapped IDs:          $(wc -l < "$OUTDIR/common.names")"
  # echo "VCFBUB-unique unmapped IDs:   $(wc -l < "$OUTDIR/vcf.unique.names")"
  # echo "CHM13-unique  unmapped IDs:   $(wc -l < "$OUTDIR/chm.unique.names")"

# extract mapq of unmapped reads in one another
  out="$OUTDIR/mapq_unmappedln_mappedCHM.txt"

  samtools view -F 2308 "$BAM_CHM" \
  | awk -v seqfile="$OUTDIR/ln.only.seq" -v out="$out" '
  function rc(s,   i,c,r){ r=""; s=toupper(s);
    for(i=length(s); i>0; i--){
      c=substr(s,i,1);
      if(c=="A") c="T"; else if(c=="T") c="A";
      else if(c=="C") c="G"; else if(c=="G") c="C";
      else c="N"; r=r c
    } return r
  }
  BEGIN{
    # Load seqfile into two sets: forward (F) and reverse-complement (R)
    while((getline s < seqfile)>0){
      gsub(/[ \t\r\n]/,"",s); if(s=="") continue
      s=toupper(s); F[s]=1; R[rc(s)]=1
    }
    close(seqfile)
  }
  {
    seq = toupper($10); mq = $5
    if (seq in F) { print mq > out; fw++ }
    else if (seq in R) { print mq > out; rcnt++ }
  }
  END{
    print "forward_matches=" fw > "/dev/stderr"
    print "rc_matches=" rcnt > "/dev/stderr"
    print "total_matches=" fw+rcnt > "/dev/stderr"
  }'



# Count how many linear-unmapped read NAMES appear anywhere in CHM.bam
  awk '
    FNR==NR {seen[$1]=1; total++; next}
    ($1 in seen) {hits++; delete seen[$1]}
    END {printf("QNAME hits=%d  total=%d  missing=%d\n", hits, total, total-hits)}
  ' <(samtools view -@8 -f 4 -F 2304 $BAM_ln) \
    <(samtools view -@8 $BAM_CHM)

  fq=(results/146/linear_chm13/trim/146/ZNF146/146Rep1.trim.pair1.fastq.gz results/146/linear_chm13/trim/146/ZNF146/146Rep1.trim.pair2.fastq.gz \
      results/146/linear_chm13/trim/146/ZNF146/146Rep2.trim.pair1.fastq.gz results/146/linear_chm13/trim/146/ZNF146/146Rep2.trim.pair2.fastq.gz)
  awk '
    FNR==1 {file++; r=0} {r++}

    # File 1: LN unmapped -> collect UNIQUE sequences (line 2 of each FASTQ record)
    file==1 && (r%4==2) {
      s=toupper($0); if(!(s in S)){ S[s]=1; total++ }; next
    }

    # Files 2..N: all raw FASTQs -> mark presence if sequence seen
    (r%4==2) {
      s=toupper($0); if(s in S) HAVE[s]=1; next
    }

    END {
      for(s in S) if(s in HAVE) hits++
      printf("LN-unmapped seqs present in RAW FASTQs: hits=%d missing=%d total=%d\n",
            hits, total-hits, total)
    }
  ' \
    <(samtools view -F 2308 "$BAM_CHM") \
    <(zcat -f "${fq[@]}")

# Troubleshooting difference in read sequence linear and graph bam (why can't just compare by matching)
  seqfile="$OUTDIR/ln.only.seq"
  NON="$OUTDIR/ln_CHMnonmatch.seq"

  ############################################
  # 1) Collect unique linear sequences not in CHM bam (neither FWD nor RC)
  ############################################
  samtools view -F 4 "$BAM_CHM" \
  | awk -v seqfile="$seqfile" -v out="$NON" '
      function rc(s,    i,c,r){
          r=""
          for(i=length(s); i>0; i--){
              c=substr(s,i,1)
              if(c=="A") c="T"
              else if(c=="T") c="A"
              else if(c=="C") c="G"
              else if(c=="G") c="C"
              else c="N"
              r = r c
          }
          return r
      }
      BEGIN{
          # load all sequences from ln.only.seq into set S[]
          while ((getline s < seqfile) > 0){
              gsub(/[ \t\r\n]/,"",s)
              if(s!=""){ S[toupper(s)]=1 }
          }
          close(seqfile)
      }
      {
          seq = toupper($10)
          if (seq in S) {
              delete S[seq]
          } else {
              rcseq = rc(seq)
              if (rcseq in S) delete S[rcseq]
          }
      }
      END{
          # print whatever is left in S (never seen in BAM)
          for (s in S) print s > out
      }
  '

  # 2 Extract ALL records from linear that carry any of those sequences

  samtools view -h -f 4 -F 2304 "$BAM_ln" \
  | awk -v list="$NON" '
  BEGIN{
    while ((getline s < list) > 0){ S[s]=1 }
    close(list)
  }
  /^@/ { print; next }
  {
    seq=toupper($10);
    if(seq in S) print;
  }
  ' | samtools view -b -o "$OUTDIR/ln_CHMnonmatch.bam" -

  # 1) Get read names from nonmatch_BAM
  samtools sort -@8 -o "$OUTDIR/treatment_alignments.sorted.bam" "$BAM_CHM"
  samtools view -@8 -h <( samtools view -@8 "$OUTDIR/ln_CHMnonmatch.bam" | cut -f1 | sort -u ) "$OUTDIR/treatment_alignments.sorted.bam" \
    | samtools sort -@8 -o $OUTDIR/ln_CHMnonmatchinCHM.bam -

  # Index the result
  samtools index $OUTDIR/ln_CHMnonmatchinCHM.bam

# Look at 1.1M read sequences become unmapped from CHM to VCF
  # Extract polymorphic L1 sequences
    VCF="Graph_genome_data/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz"
    OUTDIR="part_results/unmapped_reads_146_hardhitcap2500_hitcap100_scorefraction0.8/unmapped_compare"
    FA="$OUTDIR/full_annotation.fa"
    TSV="$OUTDIR/full_annotation.tsv"
    bcftools query -f '%CHROM\t%POS\t%ID\t%SVLEN\t%ALT\t%INFO/match_lengths\t%INFO/repeat_ids\n' "$VCF" \
    | awk -F'\t' -v OFS='\t' -v fa="$FA" -v tsv="$TSV" '
    function is_bases(s){ return s ~ /^[ACGTNacgtn]+$/ }
    function wrap60(s,  i,out){ out=""; s=toupper(s); gsub(/[ \t\r\n]/,"",s);
      for(i=1;i<=length(s);i+=60) out=out substr(s,i,60) "\n"; return out }

    BEGIN{ print "chrom","pos","id","svlen","repeat_id","match_len" > tsv }

    {
      chrom=$1; pos=$2; id=$3; svlen=$4
      alt=$5;  ml=$6;  rids=$7

      # If multiallelic, take first ALT allele
      nalt=split(alt, A, /,/); alt=A[1]

      # choose sequence from ALT only
      seq=""; if (is_bases(alt)) seq=alt
      gsub(/[ \t\r\n]/,"",seq); seq=toupper(seq)

      # choose repeat_id with largest match_length
      n1=split(ml,  L, /,/); n2=split(rids, R, /,/); n=(n1<n2?n1:n2)
      best=-1; rid="NA"
      for(i=1;i<=n;i++){ len=L[i]+0; if(len>best){ best=len; rid=R[i] } }

      print chrom,pos,id,svlen,rid,best >> tsv
      if(seq!=""){ printf(">%s|rid=%s|len=%s\n%s", id, rid, best, wrap60(seq)) >> fa }
    }
    '
    K=29
    : > "$FA"

    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/TYPE\t%INFO/LEN\n' "$VCF" \
    | awk -F'\t' -v fa="$FA" -v K="$K" '
    function is_bases(s){ return s ~ /^[ACGTNacgtn]+$/ }
    function wrap60(s,  i,out){ gsub(/[ \t\r\n]/,"",s); s=toupper(s);
      out=""; for(i=1;i<=length(s);i+=60) out=out substr(s,i,60) "\n"; return out }

    {
      chr=$1; pos=$2; id=$3; ref=$4; alt=$5; typ=$6; len=$7
      if(id=="." || id=="") id = chr ":" pos

      nalt=split(alt, A, /,/)
      ntyp=split(typ, T, /,/)
      nlen=split(len, L, /,/)

      for(i=1;i<=nalt;i++){
        a=A[i]
        if(!is_bases(a)) continue          # skip symbolic <INS>, etc.
        if(length(a) < K) continue         # too short to contribute K-mers
        t = (i<=ntyp ? T[i] : "NA")
        ln= (i<=nlen ? L[i] : length(a))
        printf(">var=%s|chr=%s|pos=%s|alle=%d|type=%s|len=%s\n%s",
              id, chr, pos, i, t, ln, wrap60(a)) >> fa
      }
    }
    '
  # Convert to FASTA
    nl -ba "$OUTDIR/vcf.only.seq" | awk '{printf(">q%s\n%s\n",$1,toupper($2))}' > $OUTDIR/vcf.only.fa
  
  # Check L1 content by alignment with minimap2 for "meaningful" hits (example: aligned length ≥ 60 and identity ≥ 80%)
    PAF="$OUTDIR/vcf.only.polyL1_hits.paf"
    TSV_ALN="$OUTDIR/vcf.only.polyL1_hits.tsv"

    minimap2 -t 8 -x sr -k15 -w5 -N 100 "$FA" "$OUTDIR/vcf.only.fa" > "$PAF"

    # PAF: col10 = nmatch, col11 = alnLen, col12 = mapq
    awk -v OFS='\t' '
      function rid_from_name(s,    m){ return match(s,/rid=([^|]+)/,m)?m[1]:"NA" }
      function var_from_name(s,    m){ return match(s,/var=([^|]+)/,m)?m[1]:"NA" }
      {
        q=$1; t=$6; nmatch=$10+0; alen=$11+0; mapq=$12+0;
        if(alen>=60 && nmatch/alen>=0.70){
          print q, var_from_name(t), rid_from_name(t), mapq, alen, nmatch/alen
        }
      }
    ' "$PAF" > "$TSV_ALN"

    # Optional quick counts per repeat_id
    awk '{c[$3]++} END{for(k in c) print k, c[k]}' "$TSV_ALN" | sort -k2,2nr > "$OUT/repeatid_counts.aln.txt"

  # Check k-mer enrichment of reads in VCF variants
    # k-mer DB from the VCF allele FASTA
    JFK="$OUTDIR/full_annotation_vcf_k29.jf"
    jellyfish count -m 29 -C -s 1G -t 8 -o "$JFK" "$OUTDIR/full_annotation.fa"
    jellyfish dump -c "$JFK" > "$OUTDIR/full_annotation_vcf_k29.counts"   # columns: kmer  count

    # per-read containment and heavy-k fraction
    K=29                          # must match jellyfish -m
    HITCAP=100                    # vg giraffe --hit-cap
    HARDCAP=2500                  # vg giraffe --hard-hit-cap
    SCOREF=0.8                    # vg giraffe --score-fraction

    TSV_KMER="$OUTDIR/full_anotation_read_kmer_metrics.tsv"

    awk -v K="$K" -v HITCAP="$HITCAP" -v HARDCAP="$HARDCAP" -v SCOREF="$SCOREF" '
      function rc(s,  i,c,r){ r=""; for(i=length(s); i; i--){ c=substr(s,i,1);
        if(c=="A")c="T"; else if(c=="T")c="A"; else if(c=="C")c="G"; else if(c=="G")c="C"; else c="N"; r=r c } return r }
      function canon(s,  t){ t=rc(s); return (s<t)?s:t }
      function pct(x,y){ return (y>0)? int((x*100.0/y)+0.5) : 0 }

      NR==FNR { C[$1]=$2+0; next }   # load k-mer counts

      # emit one row for the previous read (if any)
      function emit(){
        if(tot==0) return

        # present = k-mers with c>0
        contain_pct = pct(hit, tot)

        # split present by caps
        low  = present_low
        trim = present_trim
        drop = present_drop

        # fractions (still useful for reading; keep them)
        cap100_pct_present  = pct(trim + drop, hit)   # ≥100 over present
        cap2500_pct_present = pct(drop, hit)          # ≥2500 over present

        # counts (your request): survivors after hard cap, and survivors that pass 0.8 proxy
        survive_cnt = hit - drop                      # c < 2500
        retained_cnt = 0
        if (survive_cnt > 0) {
          cthr = cmin_survive / SCOREF                # keep if c <= cmin_survive/0.8
          for (cc in HC_survive) { c = cc + 0; if (c <= cthr) retained_cnt += HC_survive[cc] }
        }

        # print header once
        if (!header_printed){
          print "read_id",
                "tot_kmers","present_kmers",
                "containment_pct","cap100_pct_present","cap2500_pct_present",
                "survive_hardcap_count","score80_survive_count",
                "max_kmer_count","min_present_kmer_count","min_survivor_kmer_count",
                "eff_hits_100","eff_hits_2500"
          header_printed=1
        }

        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
          name,
          tot, hit,
          contain_pct, cap100_pct_present, cap2500_pct_present,
          survive_cnt, retained_cnt,
          cmax, (hit?cmin:0), (survive_cnt?cmin_survive:0),
          eff100, eff2500
      }

      # start a new read
      /^>/{
        if(NR>FNR) emit()
        name=substr($0,2)
        tot=hit=present_low=present_trim=present_drop=0
        cmax=0; cmin=0; cmin_survive=0
        eff100=eff2500=0
        delete HC_survive
        next
      }

      # sequence lines
      {
        s=toupper($0)
        for(i=1; i+K-1<=length(s); i++){
          k = canon(substr(s,i,K))
          c = (k in C)? C[k] : 0
          tot++
          if(c>0){
            hit++
            if(c>cmax) cmax=c
            if(cmin==0 || c<cmin) cmin=c

            if(c >= HARDCAP) present_drop++
            else {
              if(c >= HITCAP) present_trim++; else present_low++
              if(cmin_survive==0 || c<cmin_survive) cmin_survive=c
              HC_survive[c]++
            }

            eff100  += (c < HITCAP  ? c : HITCAP)
            eff2500 += (c < HARDCAP ? c : HARDCAP)
          }
        }
      }

      END{ emit() }
    ' "$OUTDIR/full_annotation_vcf_k29.counts" "$OUTDIR/vcf.only.fa" > "$TSV_KMER"
       # Columns:
          # 1 read_id
          # 2 containment_pct
          # 3 cap100_pct_present
          # 4 cap2500_pct_present
          # 5 score80_survive_pct
          # 6 max_kmer_count
          # 7 min_present_kmer_count
          # 8 eff_hits_100
          # 9 eff_hits_2500
          #10 tot_kmers
          #11 present_kmers

# Analyse 1.1M results
TSV="read_kmer_metrics.tsv"
# thresholds you want
CONT=50      # containment >= 50
SURV=10      # survive_hardcap_count <= 10

awk -F'\t' -v CONT="$CONT" -v SURV="$SURV" '
  NR==1{next}  # skip header
  {
    id=$1; pres=$3+0; cont=$4+0; surv=$7+0
    total++

    if (pres==0) zero++

    if (cont >= CONT) c_hit++

    # Only count survive<=SURV for reads with any present k-mers
    if (pres>0 && surv <= SURV) s_hit++

    if (cont >= CONT && pres>0 && surv <= SURV) both_hit++
  }
  END{
    printf("total\t%d\nzero_present\t%d\ncontainment_ge_%d\t%d\nsurvive_le_%d\t%d\nboth\t%d\n",
           total, zero, CONT, c_hit, SURV, s_hit, both_hit)
  }
' "$TSV"
###########################################
# Compare reads location between linear and graph
    LINEAR_BAM="$wd/results/146/linear_chm13/alignment/146/ZNF146/146.ZNF146.sorted.dup.bam"
    GRAPH_BAM="$wd/results/146/vg_giraffe_chm13/treatment_alignments.bam"
    result_dir="$wd/part_results/compare_reads_146"
    cd $result_dir || exit
    module load samtools


  # 1) Extract tuples from the linear BAM

  (
    # mate 1 (“fw”), primary only (-F 2304), MAPQ ≥ 60
    samtools view -f 64 -F 2304 -q 60 $LINEAR_BAM \
      | awk '{ print $1 "\tfw\t" $3 "\t" $4 }'
    # mate 2 (“rv”)
    samtools view -f 128 -F 2304 -q 60 $LINEAR_BAM \
      | awk '{ print $1 "\trv\t" $3 "\t" $4 }'
  ) | sort -k1,1 -k2,2 -k3,3 -k4,4 -u > linear.info
  # 2) Do the same for the graph BAM

    samtools view -F 2304 -q 60 $GRAPH_BAM \
      | awk '{ print $1 "\t" $3 "\t" $4 }'| sort  -k1,1 -k2,2 -k3,3 -u > graph.info


  # 3) Take the lines common to both
  join -t $'\t' \
    <( awk -F'\t' '{ key=$1 "|" $3; print key "\t" $0 }' linear.info | sort ) \
    <( awk -F'\t' '{ key=$1 "|" $2; print key "\t" $0 }' graph.info  | sort ) |
  awk -F'\t' '
      ($4 == $7) {                       # keep only same-chromosome matches
          # fields after join:
          #   $1 = "QNAME|CHR"
          #   $2 = QNAME      (linear)
          #   $3 = mate       ("fw"/"rv")
          #   $4 = CHR        (linear)
          #   $5 = POS_linear
          #   $6 = QNAME      (graph)
          #   $7 = CHR_graph
          #   $8 = POS_graph
          delta = $5 - $8
          print $2, $3, $4, $5, delta
      }
  ' OFS='\t' > common.info
  wc -l common.info

  awk -F'\t' '
  {
      abs = ($5 < 0 ? -$5 : $5)        # |Δ| of current line

      if (prev_q == $1 && prev_m == $2) {
          # We have a pair (two lines in a row with same read + mate)
          if (abs < prev_abs) {
              print $0                 # current one is better
          } else {
              print prev_line          # previous one is better
          }
          keep_prev = 0                # pair handled; don’t re-print
      } else {
          # New QNAME/mate; first write out the stored line (if any)
          if (keep_prev) print prev_line
          prev_line = $0
          prev_q    = $1
          prev_m    = $2
          prev_abs  = abs
          keep_prev = 1
      }
  }
  END {
      if (keep_prev) print prev_line   # flush last unpaired line
  }'  common.info > common.fixed

  awk '
  { cnt[$1]++ }                       # bump counter for each read name
  END {
      for (q in cnt) {
          if (cnt[q] == 1) n1++
          else if (cnt[q] == 2) n2++
          else if (cnt[q] == 3) n3++
          else if (cnt[q] == 4) n4++
      }
      printf "Read names appearing exactly 1× : %d\n", n1+0
      printf "Read names appearing exactly 2× : %d\n", n2+0
      printf "Read names appearing exactly 3× : %d\n", n3+0
      printf "Read names appearing exactly 4× : %d\n", n4+0
  }' linear.info

  # Count how many 0 in column 5
  awk -F'\t' '$5 == 0 { n++ } END { print n }' common.fixed

  # Extract non-zero Δ rows and sort them descending by column 5
  awk -F'\t' '$5 != 0'  common.fixed |
      sort -t$'\t' -k5,5nr  > diff_loc.info
  
  # check soft-clipped reads
  awk '$5 == -1 {print $1}' diff_loc.info > deltaminus1_readnames.txt
  module load samtools
  samtools view -N deltaminus1_readnames.txt /home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/results/146/linear_chm13/alignment/146/ZNF146/146.ZNF146.sorted.dup.bam > deltaminus1_linear.sam
  awk '$6 ~ /S/' deltaminus1_linear.sam | cut -f1,6 | less
  samtools view /home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/results/K27_FLU/vg_giraffe_vcfbub_17_7/treatment_alignments.bam "A00266:289:HMC3LDSXX:1:2171:31955:27492" | cut -f6
############################################
# Bin peaks to see if high score peaks in linear detected in graph
  cd part_results/bin_peaks_146/vcfbub
  module load bedtools
  vcfbub_bed=$wd/results/146/vg_giraffe_vcfbub/callpeaks/all_linear_peaks.bed
  chm13_bed=$wd/results/146/vg_giraffe_chm13/callpeaks/all_linear_peaks.bed
  # Add bin labels from 1 (highest scores) to 10 (lowest scores)
  awk '{
      score=$9;  # assuming score is in column 5
      if (score > 75) bin="bin1";
      else if (score > 40 && score <= 75) bin="bin2";
      else if (score > 10 && score <= 40) bin="bin3";
      else bin="bin4";
      print $0 "\t" bin;
  }' $linear_bed > linear_peaks_binned.bed


  for bin in bin1 bin2 bin3 bin4; do
      grep "${bin}$" filtered_linear_peaks.bed > ${bin}.bed
      bedtools intersect -wa -wb -a ${bin}.bed -b $graph_peak > tmp_${bin}_intersect_raw.bed
      awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$NF}' tmp_${bin}_intersect_raw.bed > graph_${bin}_intersect.bed
      echo "${bin}: $(wc -l < graph_${bin}_intersect.bed) intersected / $(wc -l < ${bin}.bed) total"
  done
  awk '{print $1}' /home/hoangnhi/projects/rrg-bourqueg-ad/hoangnhi/projects/results/iPSC_K27/linear_NCLCN-261/peak_call/iPSC_K27/iPSC_K27/iPSC_K27.iPSC_K27_peaks.narrowPeak.bed | sort -u


  bedtools intersect -v -a graph_bin1.bed -b bin1.bed > lost_peaks_vcfbubvschm13_bin1.bed
  awk '$9 < 900' lost_peaks_vcfbubvschm13.bed > filtered_lost_peaks_vcfbubvschm13.bed
  sort -k9,9nr filtered_lost_peaks_vcfbubvschm13.bed > sorted_lost_peaks_vcfbubvschm13.bed

# Check repeat contents
  bedtools intersect -a $wd/part_results/bin_peaks_146/vcfbub/graph_bin1_intersect_smaller50.bed -b $wd/Genome/t2t.fa.out.bed -wa -wb > $wd/part_results/bin_peaks_146/vcfbub/graph_bin1_intersect_smaller50_TE.bed

############################################
# Trim first 15 bases of reads
  module load seqtk
  dir=data/iPSC
  seqtk trimfq -b 8 $dir/Treatment.pair1.chr6.42M_45M.fastq.gz \
    | gzip -c > $dir/Treatment.pair1.chr6.42M_45M.15.fastq.gz
  seqtk trimfq -b 8 $dir/Treatment.pair2.chr6.42M_45M.fastq.gz \
    | gzip -c > $dir/Treatment.pair2.chr6.42M_45M.15.fastq.gz
  seqtk trimfq -b 8 $dir/p_N_iPSC_Input_r1_S18_R1_001.fastq.gz \
    | gzip -c > $dir/p_N_iPSC_Input_r1_S18_R1_001.15.fastq.gz
  seqtk trimfq -b 8 $dir/p_N_iPSC_Input_r1_S18_R2_001.fastq.gz \
    | gzip -c > $dir/p_N_iPSC_Input_r1_S18_R2_001.15.fastq.gz

################################################
# For Source code: Find line of code in all files of the directory
  grep -RIn --exclude="peaks_to_linear_matches.txt" \
    "class" "$wd/tools/graph_peak_caller/" \
    >> "$wd/tools/graph_peak_caller/peaks_to_linear_matches.txt"
  grep -RIn --exclude="qvalues_matches.txt" "qvalue" \
  "$wd/tools/graph_peak_caller/" \
  > "$wd/tools/graph_peak_caller/qvalues_matches.txt"

###################################################
# Troubleshooting peaks_to_linear
  find . -maxdepth 1 \
  -type f \
  -name "chr*_max_paths.intervalcollection" \
  ! -name "chr*_all_max_paths.intervalcollection" \
  -exec cat {} + \
  > combined.nonmax_paths.intervalcol
  mv combined.nonmax_paths.intervalcol max_peaks.intervalcollection

  if cmp -s <(sort max_peaks.intervalcollection) <(sort all_peaks.intervalcollection); then
    echo "Files are identical"
  else
    echo "Files differ"
  fi


  # Check if max_peaks.intervalcollection part of all_max_peaks.intervalcollection
  jq -c '{start, end, region_paths}' max_peaks.intervalcollection \
    | sort \
    > max_norm.json

  jq -c '{start, end, region_paths}' all_max_peaks.intervalcollection \
    | sort \
    > all_norm.json

  comm -12 max_norm.json all_norm.json  > in_both.json
  comm -23 max_norm.json all_norm.json  > only_in_max.json
  comm -13 max_norm.json all_norm.json  > only_in_all.json


###################################################
# Determine polymorphic peaks
  graph=Pangenomes/vg_giraffe/vcfbub/graphs_backbone
  

  perl -nE 'while (/"id": *"([^"]*)"/g) { say $1 }' Pangenomes/vg_giraffe/vcfbub/graphs/chr1.json | sed 's/"id": *"//;s/"$//' > Pangenomes/vg_giraffe/vcfbub/graphs/chr1_ids.txt


  comm -13 <(sort -n Pangenomes/vg_giraffe/vcfbub/graphs/chr1_ids.txt | uniq) <(sort -n Pangenomes/vg_giraffe/vcfbub/graphs_withpath/_0_GRCh38#0#chr1#0_0_248177868_ids.txt | uniq) > Pangenomes/vg_giraffe/vcfbub/polymorphic_ids.txt
###################################################
# Access change in MAPQ=0 reads when moving from linear to graph (Excel)
# Input BAM files
  LINEAR_BAM="$wd/results/K27_FLU/linear_chm13/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.sorted.dup.bam"
  GRAPH_BAM="$wd/results/K27_FLU/vg_giraffe_chm13/treatment_alignments.bam"
  result_dir="$wd/part_results/compare_reads_K27_FLU"
  cd $result_dir || exit
  module load samtools

  # Function to extract read names for a MAPQ range and primary alignments
  extract_reads() {
      bam=$1
      label=$2
      range=$3

      if [[ $range == "mapq0" ]]; then
          samtools view -F 2304 "$bam" | awk '$5 == 0 { print $1 }' | sort > "$label.txt"
      elif [[ $range == "mapq1_29" ]]; then
          samtools view -F 2304 "$bam" | awk '$5 > 0 && $5 < 30 { print $1 }' | sort > "$label.txt"
      elif [[ $range == "mapq30_59" ]]; then
          samtools view -F 2304 "$bam" | awk '$5 >= 30 && $5 < 60 { print $1 }' | sort > "$label.txt"
      elif [[ $range == "mapq60" ]]; then
          samtools view -F 2304 "$bam" | awk '$5 == 60 { print $1 }' | sort > "$label.txt"
      fi
  }

  # MAPQ bins
  extract_reads "$linear_bam" linear_mapq0 mapq0
  extract_reads "$graph_bam"  graph_mapq0  mapq0

  extract_reads "$linear_bam" linear_mapq1_29 mapq1_29
  extract_reads "$graph_bam"  graph_mapq1_29  mapq1_29

  extract_reads "$linear_bam" linear_mapq30_59 mapq30_59
  extract_reads "$graph_bam"  graph_mapq30_59  mapq30_59

  extract_reads "$linear_bam" linear_mapq60 mapq60
  extract_reads "$graph_bam"  graph_mapq60  mapq60

  # Function to compare read name files and report counts
  compare_bins() {
      label=$1
      comm -12 linear_${label}.txt graph_${label}.txt > common_${label}.txt
      comm -23 linear_${label}.txt graph_${label}.txt > unique_linear_${label}.txt
      comm -13 linear_${label}.txt graph_${label}.txt > unique_graph_${label}.txt

      echo "==== MAPQ bin: $label ===="
      echo "Common reads     : $(wc -l < common_${label}.txt)"
      echo "Unique to linear : $(wc -l < unique_linear_${label}.txt)"
      echo "Unique to graph  : $(wc -l < unique_graph_${label}.txt)"
      echo ""
  }

  # Compare all bins
  compare_bins mapq0
  compare_bins mapq1_29
  compare_bins mapq30_59
  compare_bins mapq60


###########################################
# Continue to compare their location
    LINEAR_BAM="$wd/results/K27_FLU/linear_chm13/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.sorted.dup.bam"
    GRAPH_BAM="$wd/results/K27_FLU/vg_giraffe_vcfbub_17_7/treatment_alignments.bam"
    result_dir="$wd/part_results/compare_reads_K27_FLU"
    cd $result_dir || exit
    module load samtools


  # 1) Extract tuples from the linear BAM

  (
    # mate 1 (“fw”), primary only (-F 2304), MAPQ ≥ 60
    samtools view -f 64 -F 2304 -q 60 $LINEAR_BAM \
      | awk '{ print $1 "\tfw\t" $3 "\t" $4 }'
    # mate 2 (“rv”)
    samtools view -f 128 -F 2304 -q 60 $LINEAR_BAM \
      | awk '{ print $1 "\trv\t" $3 "\t" $4 }'
  ) | sort -k1,1 -k2,2 -k3,3 -k4,4 -u > linear.info
  # 2) Do the same for the graph BAM

    samtools view -F 2304 -q 60 $GRAPH_BAM \
      | awk '{ print $1 "\t" $3 "\t" $4 }'| sort  -k1,1 -k2,2 -k3,3 -u > graph.info


  # 3) Take the lines common to both
  join -t $'\t' \
    <( awk -F'\t' '{ key=$1 "|" $3; print key "\t" $0 }' linear.info | sort ) \
    <( awk -F'\t' '{ key=$1 "|" $2; print key "\t" $0 }' graph.info  | sort ) |
  awk -F'\t' '
      ($4 == $7) {                       # keep only same-chromosome matches
          # fields after join:
          #   $1 = "QNAME|CHR"
          #   $2 = QNAME      (linear)
          #   $3 = mate       ("fw"/"rv")
          #   $4 = CHR        (linear)
          #   $5 = POS_linear
          #   $6 = QNAME      (graph)
          #   $7 = CHR_graph
          #   $8 = POS_graph
          delta = $5 - $8
          print $2, $3, $4, $5, delta
      }
  ' OFS='\t' > common.info
  wc -l common.info

  awk -F'\t' '
  {
      abs = ($5 < 0 ? -$5 : $5)        # |Δ| of current line

      if (prev_q == $1 && prev_m == $2) {
          # We have a pair (two lines in a row with same read + mate)
          if (abs < prev_abs) {
              print $0                 # current one is better
          } else {
              print prev_line          # previous one is better
          }
          keep_prev = 0                # pair handled; don’t re-print
      } else {
          # New QNAME/mate; first write out the stored line (if any)
          if (keep_prev) print prev_line
          prev_line = $0
          prev_q    = $1
          prev_m    = $2
          prev_abs  = abs
          keep_prev = 1
      }
  }
  END {
      if (keep_prev) print prev_line   # flush last unpaired line
  }'  common.info > common.fixed

  awk '
  { cnt[$1]++ }                       # bump counter for each read name
  END {
      for (q in cnt) {
          if (cnt[q] == 1) n1++
          else if (cnt[q] == 2) n2++
          else if (cnt[q] == 3) n3++
          else if (cnt[q] == 4) n4++
      }
      printf "Read names appearing exactly 1× : %d\n", n1+0
      printf "Read names appearing exactly 2× : %d\n", n2+0
      printf "Read names appearing exactly 3× : %d\n", n3+0
      printf "Read names appearing exactly 4× : %d\n", n4+0
  }' common.fixed

  # Count how many 0 in column 5
  awk -F'\t' '$5 == 0 { n++ } END { print n }' common.fixed

  # Extract non-zero Δ rows and sort them descending by column 5
  awk -F'\t' '$5 != 0'  common.fixed |
      sort -t$'\t' -k5,5nr  > diff_loc.info
  awk -F'\t' '($5 >= -50 && $5 <= 50) { n++ } END { print n }' diff_loc.info

###########################################
# See how many hits per seed per read
# Cannot work this way, change strategy to vg autoindex to construct the graph again
  head -n400000 data/iPSC/Treatment1.trim.pair2.fastq.gz | gzip > data/iPSC/Treatment1.trim.pair2_100k.fastq.gz
  vg kmers -k 17  -t 12 $graph | sort -u > Pangenomes/vg_giraffe/vcfbub/chr6.42M_45M.k17.kmers.tsv

  # # 2) now run the same awk script, but point it at graph.k17.kmers.txt
  k=17
  kmers_file="Pangenomes/vg_giraffe/vcfbub/chr6.42M_45M.k17.kmers.tsv"
  fastq="data/iPSC/Treatment1.trim.pair1_100k.fastq.gz"
  out="part_results/current/iPSC_k${k}_seed_counts.txt"
  mkdir -p "$(dirname "$out")"

  # zcat into awk in C locale
  zcat "$fastq" \
    | LC_ALL=C awk -v k="$k" -v M="$kmers_file" '
      BEGIN {
        FS = "\t"
        # build lookup of only the kmer (field 1) from the .kmers file
        while ((getline < M) > 0) {
          km[$1] = 1
        }
        close(M)
        # for debug: print how many kmers we loaded
        print "loaded", length(km), "kmers from", M > "/dev/stderr"
        # and maybe print the first few:
        c = 0
        for (x in km) {
          print x > "/dev/stderr"
          if (++c == 10) break
        }
      }
      # seq lines in fastq are line 2 of every 4
      NR % 4 == 2 {
      seq = $0
      L = length(seq)
      # debug: show read length vs k
      # print "DEBUG read length=" L ", k=" k > "/dev/stderr"
      hits = 0

      # only iterate if L >= k
      for (i = 1; i <= L - k + 1; i++) {
        s = substr(seq, i, k)
        if (s in km) {
          hits++
          # debug first few matches
          # if (hits <= 3) {
          #   print "DEBUG match " hits ": " s " @pos " i > "/dev/stderr"
          # }
        }
      }

      # debug: total hits for this read
      print "DEBUG total hits=" hits > "/dev/stderr"

      # output the count for downstream
      print hits
    }
  ' > "$out"

###########################################
# Create small dataset (aligned fastq) in the chosen range
  bam="results/K27_FLU/linear_chm13/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.sorted.bam"
  region="chr6:42000000-45000000"
  samtools view -b "$bam" "$region" > "results/K27_FLU/linear_chm13/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.chr6.42M_45M.bam"
  # Grab just the reads in that region (names only)
  samtools view -b "$bam" "$region" \
    | samtools fastq \
        -1 data/cgroza_data/H3K27AC_FLU/treatment/Treatment.pair1.chr6.42M_45M.fastq.gz \
        -2 data/cgroza_data/H3K27AC_FLU/treatment/Treatment.pair2.chr6.42M_45M.fastq.gz \
        -0 /dev/null  \
        -n   # preserve original read names exactly
#
# Create small graph and different k and w min index files
  # Create range vcf file
  bcftools view -r chr6:42e6-45e6 $wd/Graph_genome_data/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz \
    -Oz -o chr6.42M_45M.vcf.gz
  # Create range fa file, change the header to "chr6" to match the vcf, then reindex
  # Check header by grep '^>' $wd/Graph_genome_data/chr6.42_45M.fa
  samtools faidx $wd/Graph_genome_data/chm13v2.0.fa chr6:42000000-45000000 > $wd/Graph_genome_data/chr6_42-45M.fa
  sed 's/^>chr6:.*/>chr6/' $wd/Graph_genome_data/chr6.42_45M.fa > $wd/Graph_genome_data/chr6.42_45M.fixed.fa
  samtools faidx $wd/Graph_genome_data/chr6.42M_45M.fixed.fa
  # vg index and create minimizer index
  vg autoindex --workflow giraffe -r $wd/Graph_genome_data/chr6.42_45M.fixed.fa -v $wd/Graph_genome_data/chr6.42M_45M.vcf.gz -p chr6.42M_45M --threads 64
  mv chr6.42M_45M.giraffe.gbz chr6.42M_45M.gbz
  vg minimizer -k 9 -w 4 -d chr6.42M_45M.dist chr6.42M_45M.gbz -o chr6.42M_45M.k9w4.min
  vg minimizer -k 11 -w 5 -d chr6.42M_45M.dist chr6.42M_45M.gbz -o chr6.42M_45M.k11w5.min
  vg minimizer -k 13 -w 6 -d chr6.42M_45M.dist chr6.42M_45M.gbz -o chr6.42M_45M.k13w6.min
  vg minimizer -k 15 -w 7 -d chr6.42M_45M.dist chr6.42M_45M.gbz -o chr6.42M_45M.k15w7.min
