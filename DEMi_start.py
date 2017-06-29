import sys
import os

'''
Take inputs: Fastq file and Barcodes File
Read files into usable data structures
'''
fastq_handle=sys.argv[1]
barcodes_handle=sys.argv[2]

output_directory=sys.argv[3]
if output_directory[-1]!="/":
    output_director+="/"

#Process Barcodes: Into  a List and Dictionary
barcodes_file=open(barcodes_handle)
barcodes_file_lines=barcodes_file.readlines()
barcodes_file.close()

barcode_list=[]
barcode_dicts={}


for barcode in barcodes_file_lines:
    barcode_list.append(barcode.strip("\n"))
    barcode_dicts[barcode.strip("\n")]=[]


fastq_file=open(fastq_handle)
fastq_lines=fastq_file.readlines()
fastq_file.close()

total_reads=len(fastq_lines)/4

matched_reads=[]
unmatched_reads=[]


read_starts=[]
cut_sites=[]

for index in range(total_reads):

    raw_read_line=fastq_lines[1+4*index]
    raw_read=raw_read_line.strip("\n")

    read_start=raw_read[0:6]
    if read_start in barcode_list:
        barcode_dicts[read_start].append(raw_read)
        matched_reads.append(raw_read)
    else: 
        unmatched_reads.append(raw_read)
    
    cut_site=raw_read[6:10]
    if cut_site not in cut_sites:
        cut_sites.append(cut_site)

for each_barcode in barcode_dicts:
    barcode_file=open(output_directory+each_barcode,"w")
    for assigned_reads in barcode_dicts[each_barcode]:
        barcode_file.write(assigned_reads+"\n")
    barcode_file.close()


barcode_file=open(output_directory+"unmatched","w")
for reads in unmatched_reads:
    barcode_file.write(reads+"\n")
barcode_file.close()

print cut_sites
