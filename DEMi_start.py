import argparse
import collections
import sys
import os
import operator

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

def fastq_generator(R1_fastq_file,R2_fastq_file):

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

    print "closing files..."
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

def get_cutsites(R1fastq,R2fastq,BC_len,cutsite_len):
    R1_cutsite=R1fastq.sequence[BC_len:BC_len+cutsite_len]
    R2_cutsite=R2fastq.sequence[BC_len:BC_len+cutsite_len]
    
    return R1_cutsite,R2_cutsite

def barcodefile_dict_maker(R1fastq,R2fastq,barcodes,out_dir):
    R1_barcode_dict={}
    R2_barcode_dict={}

    R1fastq_name=R1fastq.strip("R1.fastq")
    R2fastq_name=R2fastq.strip("R2.fastq")

    # Name demultiplexed cells according to original fastq, end in .fastq 
    for each_barcode in barcodes:
        R1_barcode_dict[each_barcode]=open(out_dir+R1fastq_name+"_"+each_barcode+"_R1.fastq","w")
        R2_barcode_dict[each_barcode]=open(out_dir+R2fastq_name+"_"+each_barcode+"_R2.fastq","w")
        
    
    R1_barcode_dict["unmatched"]=open(out_dir+R1fastq_name+"_unmatched_R1.fastq","w")
    R2_barcode_dict["unmatched"]=open(out_dir+R2fastq_name+"_unmatched_R2.fastq","w")

    return R1_barcode_dict,R2_barcode_dict

def generate_reports(out_dir,R1barcode_counter, R1cutsite_counter, R2cutsite_counter):

    #Barcode counting and site counting reporting. Is this awkward?
    R1barcode_final=sorted(R1barcode_counter.items(),key=operator.itemgetter(1),reverse=True)
    R1cutsite_final=sorted(R1cutsite_counter.items(),key=operator.itemgetter(1),reverse=True)
    R2cutsite_final=sorted(R2cutsite_counter.items(),key=operator.itemgetter(1),reverse=True)   

    if not os.path.exists(out_dir+"/reports"):
        os.mkdir(out_dir+"/reports")

    barcode_report=open(out_dir+"reports/barcode_report","w")
    barcode_report.write("barcode \t counts \n")
    for barcodes in R1barcode_final:
        barcode_report.write(barcodes[0]+"\t"+str(barcodes[1])+"\n")
    barcode_report.close()

    R1cutsite_report=open(out_dir+"reports/R1cutsite_report","w")
    R1cutsite_report.write("cutsite \t counts \n")
    for cutsites in R1cutsite_final:
        R1cutsite_report.write(cutsites[0]+"\t"+str(cutsites[1])+"\n")
    R1cutsite_report.close()


    R2cutsite_report=open(out_dir+"reports/R2cutsite_report","w")
    R2cutsite_report.write("cutsite \t counts \n")
    for cutsites in R2cutsite_final:
        R2cutsite_report.write(cutsites[0]+"\t"+str(cutsites[1])+"\n")
    R2cutsite_report.close()
    





def main():
    '''
    parser=argparse.ArgumentParser(description="Demultiplexing program")
    parser.add_argument("R1.fastq", type=str, help="R1 fastqs")
    parser.add_argument("R2.fastq", type=str, help="R2 fastqs")
    parser.add_argument("output_dir", type=str, help="output directory") 
    parser.add_argument("--barcodes", type=str, help="barcode file")
   
    args=parser.parse_args()

    print args
    ''' 
 
    fastq_files, barcodes_file, out_dir, BC_length=arg_parse()
    
    if not os.path.isdir(out_dir):
        os.system("mkdir "+out_dir)

    barcodes=load_barcodes(barcodes_file)
    R1_barcodefile_dict,R2_barcodefile_dict=barcodefile_dict_maker(fastq_files[0],fastq_files[1],barcodes,out_dir)

    R1barcode_counter={}
    R1cutsite_counter={}
    R2cutsite_counter={}

    
    R1_fastq_file=open(fastq_files[0])
    R2_fastq_file=open(fastq_files[1])

    # Switch to fastq_generator(*fastq_files)... ? Look into star operator. 
    # ... packing and unpacking collections
    # Populating Barcode and Cutsite dictionarys is awkward
    for full_read in fastq_generator(R1_fastq_file,R2_fastq_file):
        R1_set=full_read[0]
        R2_set=full_read[1]
        read_start=get_read_start(R1_set,BC_length)
        barcode_match=check_for_barcode(read_start,barcodes)
        if barcode_match in R1barcode_counter:
            R1barcode_counter[barcode_match]+=1
        else:
            R1barcode_counter[barcode_match]=1
        
        #Package this in its own function?    
        R1cutsite,R2cutsite=get_cutsites(R1_set,R2_set,BC_length,5)        
        if R1cutsite in R1cutsite_counter:
            R1cutsite_counter[R1cutsite]+=1
        else:
            R1cutsite_counter[R1cutsite]=1
        if R2cutsite in R2cutsite_counter:
            R2cutsite_counter[R2cutsite]+=1
        else:
            R2cutsite_counter[R2cutsite]=1

        R1_to_write="".join(R1_set)
        R2_to_write="".join(R2_set)
        R1_barcodefile_dict[barcode_match].write(R1_to_write)
        R2_barcodefile_dict[barcode_match].write(R1_to_write)

    R1_fastq_file.close()
    R2_fastq_file.close()

    for barcodefiles in R1_barcodefile_dict.keys():
        R1_barcodefile_dict[barcodefiles].close()

    for barcodefiles in R2_barcodefile_dict.keys():
        R2_barcodefile_dict[barcodefiles].close()

    generate_reports(out_dir,R1barcode_counter,R1cutsite_counter,R2cutsite_counter)



if __name__=="__main__":
    main()
