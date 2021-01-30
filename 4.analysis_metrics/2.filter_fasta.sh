python analysis/4.analysis_metrics/get_longest_transcript.py -i data/fasta/raw_fasta/B73v4.fa  -o data/fasta/filt_fasta/B73v4.fa -s "Zm00001d" -p "_T[0-9]+"
python analysis/4.analysis_metrics/get_longest_transcript.py -i data/fasta/raw_fasta/W22.fa    -o data/fasta/filt_fasta/W22.fa   -s "Zm00004b" -p "_P[0-9]+"
python analysis/4.analysis_metrics/get_longest_transcript.py -i data/fasta/raw_fasta/PH207.fa  -o data/fasta/filt_fasta/PH207.fa -s "Zm00008a" -p "_P[0-9]+"
python analysis/4.analysis_metrics/get_longest_transcript.py -i data/fasta/raw_fasta/Mo17.fa   -o data/fasta/filt_fasta/Mo17.fa -s "Zm00014a" -p "_P[0-9]+"
python analysis/4.analysis_metrics/get_longest_transcript.py -i data/fasta/raw_fasta/B73v3.fa  -o data/fasta/filt_fasta/B73v3_GRMZM.fa -s "^GRMZM" -p "_P[0-9]+"
python analysis/4.analysis_metrics/get_longest_transcript.py -i data/fasta/raw_fasta/B73v3.fa  -o data/fasta/filt_fasta/B73v3_AC.fa -s "^(?!GRMZM).*" -p "FGP" -r "FG"

cat data/fasta/filt_fasta/B73v3_GRMZM.fa data/fasta/filt_fasta/B73v3_AC.fa >  data/fasta/filt_fasta/B73v3.fa && \
rm data/fasta/filt_fasta/B73v3_GRMZM.fa data/fasta/filt_fasta/B73v3_AC.fa