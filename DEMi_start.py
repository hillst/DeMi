import collections
import sys
import os

'''
Take inputs: Fastq file and Barcodes File
Read files into usable data structures
'''
def arg_parse():


    #check for help request arguments 
    help_requests=["-h","--h","-help","--help"]

    help_check=set(sys.argv).intersection(set(help_requests))
    if len(help_check)>0:
        print "help message"
        print "Arg1: fastq, Arg2: barcodes, Arg3: output directory"
        sys.exit()

    if len(sys.argv)<4:
        print "Error: Not enough arguments"
        sys.exit()

    R1fastq_handle=sys.argv[1]
    R2fastq_handle=sys.argv[2]
    barcodes_handle=sys.argv[3]
    output_directory=sys.argv[4]

       
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

    argument_file_list=sys.argv[1:4]

    #Make sure all file/folder arguments are valid paths
    for argument in [R1fastq_handle,R2fastq_handle,barcodes_handle,output_directory]:
        if os.path.exists(argument)==False:
            print "error- file : ", argument
            print "no such path!!"
            sys.exit()


    if output_directory[-1]!="/":
        output_directory+="/"
    
    return([R1fastq_handle,R2fastq_handle],barcodes_handle,output_directory,BC_length)

def load_barcodes(barcodes):
    barcodes_file=open(barcodes)
    barcodes_file_lines=barcodes_file.readlines()
    barcodes_file.close()

    barcode_list=[]

    for barcode in barcodes_file_lines:
        barcode_list.append(barcode.strip("\n"))


    return barcode_list

def fastq_generator(R1_fastq,R2_fastq):
    R1_fastq_file=open(R1_fastq)
    R2_fastq_file=open(R2_fastq)

    read_description=collections.namedtuple("read_set","ID sequence plus quality")
    
    while True:

        R1_ID=R1_fastq_file.readline()
        R2_ID=R2_fastq_file.readline()    
            
        if len(R1_ID)==0:
            break
        
        R1_sequence=R1_fastq_file.readline()
        R2_sequence=R2_fastq_file.readline()
        R1_plus=R1_fastq_file.readline()
        R2_plus=R2_fastq_file.readline()
        R1_quality=R1_fastq_file.readline()
        R2_quality=R2_fastq_file.readline()

        R1_read=read_description(ID=R1_ID,sequence=R1_sequence,plus=R1_plus,quality=R1_quality)
        R2_read=read_description(ID=R2_ID,sequence=R2_sequence,plus=R2_plus,quality=R2_quality)

        if R1_ID.split(" ")[0]!=R2_ID.split(" ")[0]:
            raise Exception("R1 and R2 have non-matching IDs: \nR1: "+R1_ID+"\nR2: "+R2_ID)

        yield R1_read, R2_read

    R1_fastq_file.close()
    R2_fastq_file.close()         

def get_read_start(four_line_read,BClength):
    read_start=four_line_read.sequence[0:BClength]
    return read_start    

def check_for_barcode(read_start,barcode_list):

    if read_start in barcode_list:
        read_result=read_start
    else: 
        read_result="unmatched"
    return read_result



def barcodefile_dict_maker(barcodes,out_dir):
    R1_barcode_dict={}
    R2_barcode_dict={}

    for each_barcode in barcodes:
        R1_barcode_dict[each_barcode]=open(out_dir+"R1_"+each_barcode,"w")
        R2_barcode_dict[each_barcode]=open(out_dir+"R2_"+each_barcode,"w")
        
    
    R1_barcode_dict["unmatched"]=open(out_dir+"R1_unmatched","w")
    R2_barcode_dict["unmatched"]=open(out_dir+"R2_unmatched","w")

    return R1_barcode_dict,R2_barcode_dict

def main():
    
    fastq_files, barcodes_file, out_dir, BC_length=arg_parse()
    
    if not os.path.isdir(out_dir):
        os.system("mkdir "+out_dir)
    barcodes=load_barcodes(barcodes_file)
     
    R1_barcodefile_dict,R2_barcodefile_dict=barcodefile_dict_maker(barcodes,out_dir)
    

    for full_read in fastq_generator(fastq_files[0],fastq_files[1]):
        R1_set=full_read[0]
        R2_set=full_read[1]
        read_start=get_read_start(R1_set,BC_length)
        barcode_match=check_for_barcode(read_start,barcodes)
        R1_to_write="".join(R1_set)
        R2_to_write="".join(R2_set)
        R1_barcodefile_dict[barcode_match].write(R1_to_write)
        R2_barcodefile_dict[barcode_match].write(R1_to_write)

    for barcodefiles in R1_barcodefile_dict.keys():
        R1_barcodefile_dict[barcodefiles].close()

    for barcodefiles in R2_barcodefile_dict.keys():
        R2_barcodefile_dict[barcodefiles].close()

if __name__=="__main__":
    main()
