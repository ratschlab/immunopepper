#!/usr/bin/env bash
# bash for use case run
immunopepper build --mutation_mode ref --samples ENCSR000BZG ERR2130621 --output_dir immunopepper_usecase_out --splice_path immunopepper_usecase.pickle --ann_path immunopepper_usecase.gtf --ref_path genome1.fa --kmer 9 --count_path immunopepper_usecase.count.hdf5
immunopepper build --mutation_mode germline --samples ENCSR000BZG ERR2130621 --output_dir immunopepper_usecase_out --splice_path immunopepper_usecase.pickle --ann_path immunopepper_usecase.gtf --ref_path genome1.fa --kmer 9 --count_path immunopepper_usecase.count.hdf5 --vcf_path immunopepper_usecase.vcf --maf_path immunopepper_usecase.maf
immunopepper build --mutation_mode somatic --samples ENCSR000BZG ERR2130621 --output_dir immunopepper_usecase_out --splice_path immunopepper_usecase.pickle --ann_path immunopepper_usecase.gtf --ref_path genome1.fa --kmer 9 --count_path immunopepper_usecase.count.hdf5 --vcf_path immunopepper_usecase.vcf --maf_path immunopepper_usecase.maf
immunopepper build --mutation_mode somatic_and_germline --samples ENCSR000BZG ERR2130621 --output_dir immunopepper_usecase_out --splice_path immunopepper_usecase.pickle --ann_path immunopepper_usecase.gtf --ref_path genome1.fa --kmer 9 --count_path immunopepper_usecase.count.hdf5 --vcf_path immunopepper_usecase.vcf --maf_path immunopepper_usecase.maf

# create background from ENCSR000BZG
immunopepper filter --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ENCSR000BZG/ref_mode_background_kmer.tsv --junction_kmer_tsv_path immunopepper_usecase_out/ENCSR000BZG/ref_junction_kmer.txt --junc_expr
immunopepper make_bg --kmer_files_list immunopepper_usecase_out/ENCSR000BZG/ref_mode_background_kmer.tsv --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/background_kmer.txt

# get uniq kmer (ref)
immunopepper diff --junction_kmer_file  immunopepper_usecase_out/ERR2130621/ref_junction_kmer.txt --bg_file_path immunopepper_usecase_out/background_kmer.txt --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/ref_junction_kmer_remove_bg.tsv --remove_bg
immunopepper filter --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/ref_junction_kmer_remove_bg_filter.tsv --junction_kmer_tsv_path immunopepper_usecase_out/ERR2130621/ref_junction_kmer_remove_bg.tsv --cross_junction --seg_expr --seg_expr_thre 2

immunopepper diff --junction_kmer_file  immunopepper_usecase_out/ERR2130621/germline_junction_kmer.txt --bg_file_path immunopepper_usecase_out/background_kmer.txt --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/germline_junction_kmer_remove_bg.tsv --remove_bg
immunopepper filter --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/germline_junction_kmer_remove_bg_filter.tsv --junction_kmer_tsv_path immunopepper_usecase_out/ERR2130621/germline_junction_kmer_remove_bg.tsv --cross_junction --seg_expr --seg_expr_thre 2

immunopepper diff --junction_kmer_file  immunopepper_usecase_out/ERR2130621/somatic_junction_kmer.txt --bg_file_path immunopepper_usecase_out/background_kmer.txt --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/somatic_junction_kmer_remove_bg.tsv --remove_bg
immunopepper filter --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/somatic_junction_kmer_remove_bg_filter.tsv --junction_kmer_tsv_path immunopepper_usecase_out/ERR2130621/somatic_junction_kmer_remove_bg.tsv --cross_junction --seg_expr --seg_expr_thre 2

immunopepper diff --junction_kmer_file  immunopepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer.txt --bg_file_path immunopepper_usecase_out/background_kmer.txt --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer_remove_bg.tsv --remove_bg
immunopepper filter --output_dir immunopepper_usecase_out --output_file_path immunopepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer_remove_bg_filter.tsv --junction_kmer_tsv_path immunopepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer_remove_bg.tsv --cross_junction --seg_expr --seg_expr_thre 2

tail -n +2 immunopepper_usecase_out/ERR2130621/*_junction_kmer_remove_bg_filter.tsv | cat | grep -v "==>" | cut -f1 | sort |uniq | grep . > neo_kmer.txt
