import sys
import os

'''
Take inputs: Fastq file and Barcodes File
Read files into usable data structures
'''
def arg_parse():
    fastq_handle=sys.argv[1]
    barcodes_handle=sys.argv[2]

    output_directory=sys.argv[3]
    if output_directory[-1]!="/":
        output_directory+="/"
    
    return(fastq_handle,barcodes_handle,output_directory)

def load_barcodes(barcodes):
    barcodes_file=open(barcodes)
    barcodes_file_lines=barcodes_file.readlines()
    barcodes_file.close()

    barcode_list=[]

    for barcode in barcodes_file_lines:
        barcode_list.append(barcode.strip("\n"))


    return barcode_list

def fastq_generator(fastq):
    fastq_file=open(fastq)
    read_set=[]
    for index, lines in enumerate(fastq_file):
        line=index+1
        read_set.append(lines)
        if line%4==0:
            yield read_set
            read_set=[]

def get_read_start(four_line_read):
    read_start=four_line_read[1][0:6]
    return read_start    

def check_for_barcode(read_start,barcode_list):

    if read_start in barcode_list:
        read_result=read_start
    else: 
        read_result="unmatched"
    return read_result

def write_to_file(out_dir,write_barcode,full_read_list):
    full_read=""
    for line in full_read_list:
        full_read+=line
    full_read=full_read.strip("\n")
    print full_read
    barcode_file=open(out_dir+write_barcode,"a")
    barcode_file.write(full_read)
    barcode_file.close()


def main():
    
    fastq_file, barcodes_file, out_dir=arg_parse()
    
    if not os.path.isdir(out_dir):
        os.system("mkdir "+out_dir)

    barcodes=load_barcodes(barcodes_file)
    
    fastq_items=fastq_generator(fastq_file)
   
    counter=0 

    for full_reads in fastq_items:
        full_read=next(fastq_items)
        read_start=get_read_start(full_read)
        barcode_match=check_for_barcode(read_start,barcodes)
        write_to_file(out_dir,barcode_match,full_read)

main()
