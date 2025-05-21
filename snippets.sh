#!/bin/bash

# This file is to keep common snippets used but are not included in the pipeline

###################################################
# Access change in MAPQ=0 reads when moving from linear to graph (Excel)
# Input BAM files
  linear_bam="$wd/results/146/vg_giraffe_chm13/treatment_alignments.bam"
  graph_bam="$wd/results/146/vg_giraffe_vcfbub/treatment_alignments.bam"
  result_dir="$wd/part_results/current"
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
  bam="results/iPSC_K27/linear_chm13/alignment/iPSC_K27/iPSC_K27/iPSC_K27.iPSC_K27.sorted.dup.filtered.bam"
  region="chr6:42000000-45000000"
  fastq1="data/iPSC/Treatment1.trim.pair2.fastq.gz"
  out_fastq="data/iPSC/Treatment1.trim.pair2.chr6.42M_45M.fastq.gz"

  # 2) Grab just the reads in that region (names only)
  samtools view -@4 "$bam" "$region" \
    | cut -f1 \
    | sort -u \
    > data/iPSC/chr6.42M_45M_read_ids.txt

  # 3) Pull those IDs back out of the original FASTQ1
  #    Using seqtk (fast and memory‐efficient)
  module load seqtk
  seqtk subseq "$fastq1" data/iPSC/chr6.42M_45M_read_ids.txt \
    | gzip -c > "$out_fastq"

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