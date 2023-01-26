import sys
import os
import csv
import re

def parse_vcf(variants, reader, sample):
	for line in reader:
		if line[0].startswith('#') is False:
			variants.append(read_line_vcf(line, sample))
	
	return variants

def read_line_vcf(line, sample):
	variant = []

	# sample name
	variant.append(sample)
	# variant position
	variant.append(int(line[1]))
	# variant type
	variant.append(re.search('[A-Z]{3}', re.search('=[A-Z]{3}', re.search('SVTYPE=[A-Z]{3}', line[7]).group(0)).group(0)).group(0))
	# variant length
	variant.append(re.search('[0-9]+', re.search('SVLEN=-?[0-9]+', line[7]).group(0)).group(0))
	# variant sequence
	if variant[2] == "DEL":
		variant.append(line[3])
	elif variant[2] == "INS":
		variant.append(line[4])
	
	return variant

def main():
	extension = sys.argv[1]
	variants = []

	for root, directories, files in os.walk(sys.argv[2]):
		for fi in files:
			if fi.endswith(extension):
				# extract sample name in file name
				sample = str(fi.split('.')[0] + '.' + fi.split('.')[1])

				# open vcf file in csv reader
				vcfFile = open(os.path.join(root, fi))
				vcfReader = csv.reader(vcfFile, delimiter='\t')

				# extract SV from vcf file
				variants = parse_vcf(variants, vcfReader, sample)

				vcfFile.close()
	
	variants = sorted(variants, key=lambda x: (x[1], x[0]))
	insertions = [variant for variant in variants if variant[2] == "INS"]
	deletions = [variant for variant in variants if variant[2] == "DEL"]

        print("")

	print("insertions")
	for v in insertions:
		print(v[0], v[1], v[2], v[3], v[4], sep="\t", end="\n")
	print("deletion")
	for v in deletions:
		print(v[0], v[1], v[2], v[3], v[4], sep="\t", end="\n")

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: extract_sv.py <extension> <res_directory>")
		exit()
	main()
