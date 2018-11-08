'''
In some cases, the output does not exist. So a special character is used in the output file.
Mainly on the flag field and expression counts field in the output .tsv file
For example,

- flag field IS_IN_JUNCTION_LIST if no external junction file is provided
- flag field JUNCTION_ANNOTATED if it is an isolated exon
- field VARIANT_COMBINE if no variant existing in the exon-pair
- expression field VARIANT_SEG_EXPR if no variant existing in the exon-pair
- expression field JUNCTION_EXPR if it is an isolated exon
'''

NOT_EXIST = '.'
