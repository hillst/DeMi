import collections
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

    argument_list=[fastq_handle,barcodes_handle,output_directory]
    
    #Set default barcode length, adjust if argument given    
    BC_length=6
    if "-barcode_length" in sys.argv:
        BC_option_index=sys.argv.index("-barcode_length")
        if len(sys.argv)<BC_option_index+2:
            print "Error: incorrect -barcode_length argument formatting"
            sys.exit()
        BC_length=sys.argv[BC_option_index+1]
        if not type(BC_length.isdigit):
            print "Error: incorrect -barcode_length argument formatting"
    
    BC_length=int(BC_length)

   
    #check for help request arguments 
    help_requests=["-h","--h","-help","--help"]
    help_check=set(argument_list).intersection(set(help_requests))
    if len(help_check)>0:
        print "help message"
        print "Arg1: fastq, Arg2: barcodes, Arg3: output directory"
        sys.exit()

    #Make sure all file/folder arguments are valid paths
    for argument in argument_list:
        if os.path.exists(argument)==False:
            print "error- file : ", argument
            print "no such path!!"
            sys.exit()


    if output_directory[-1]!="/":
        output_directory+="/"
    
    return(fastq_handle,barcodes_handle,output_directory,BC_length)

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
    read_description=collections.namedtuple("read_set","ID sequence plus quality")
    read_set=[]   
    
    #FIX THIS
    this_quality="start" 

    while this_quality:
        this_ID=fastq_file.readline()
        this_sequence=fastq_file.readline()
        this_plus=fastq_file.readline()
        this_quality=fastq_file.readline()

        read=read_description(ID=this_ID,sequence=this_sequence,plus=this_plus,quality=this_quality)
        yield read
         

def get_read_start(four_line_read,BClength):

    read_start=four_line_read.sequence[0:BClength]
    return read_start    

def check_for_barcode(read_start,barcode_list):

    if read_start in barcode_list:
        read_result=read_start
    else: 
        read_result="unmatched"
    return read_result

def write_to_file(out_dir,barcode_filename,full_read_list):
    full_read=""
    for line in full_read_list:
        full_read+=line

    barcode_file=open(out_dir+barcode_filename,"a")
    barcode_file.write(full_read)
    barcode_file.close()


def main():
    
    fastq_file, barcodes_file, out_dir, BC_length=arg_parse()
    
    if not os.path.isdir(out_dir):
        os.system("mkdir "+out_dir)
    barcodes=load_barcodes(barcodes_file)
    
    fastq_items=fastq_generator(fastq_file)

    for full_read in fastq_items:
        read_start=get_read_start(full_read,BC_length)
        barcode_match=check_for_barcode(read_start,barcodes)
        write_to_file(out_dir,barcode_match,full_read)


if __name__=="__main__":
    main()
