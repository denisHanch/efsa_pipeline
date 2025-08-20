#!/usr/bin/env python3
import validation.utils as u
import validation.format as r
from Bio import SeqIO


class Status:
    def __init__(self):
        self.value = "SUCCESS"
        self.msg = None
    
    def success(self) ->bool:
        return self.value == "SUCCESS"
     
    def get_status(self):
        if self.success():
            #   SUCCESS
            return self.value,None
        else:
            #   ERROR
            return self.value,self.msg
        
    def log(self, e:Exception):
        self.value = type(e).__name__
        self.msg = str(e)

    def __str__(self):
        return f"Status(value={self.value}, msg={self.msg})"



class Validator:

    def __init__(self, filepath:str, file_format:str, encode_with:str=None, outdir:str=None, outsuffix:str=None, edit_specification:dict=None):
        """To validate, edit and format input file.

        Args:
            filepath (str): Path to the file that is to be validated and formatted
            file_format (str): File format to be on output
            encode_with (str, optional): Coding type to be on output. Defaults to None.
            outdir (str, optional): Outdir may be specified. Default will be setup to the same directory as the input is.
            outsuffix (str, optional): To tag files with suffix like <filename>_<suffix>.<ext>. Defaults to None.
        """
        
        #   Input specifications
        ##   By user
        self.filepath = filepath
        # self.edit_specification = edit_specification    #   dictionary that setup which editation has to be done on FASTA and with what arguments
        ##   By system
        self._allowed_exts = ()
        self._normalized_ext = ()
        self._encoded_with = None     #   compression technique applied in input file
        self._file_format = None      #   file format of input file, expected - FASTA/FASTQ/BAM/GTF/GFF/...? - defined by subclass
        self._basedir = None          #   directory containing input file   
        self._filename = None         #   input file name
        self._work_file = None        #   current file to work with (initialy input file, but may change based on the encoding and moving the original file)
        
        #   Output specifications
        ##  By user
        self.encode_with = encode_with  #   Compression method expected gz,bz2
        self.outdir = outdir            #   directory to save in output file
        self.outsuffix = outsuffix      #   ref,mod or any other marker?
        self.file_format = file_format.upper()  #   output file format, expected - FASTA/FASTQ
        ##  By system
        self._ext = None             #   output file extension
        self._outpath = None         #   output file absolute path
        self._reformat = False       #   bool to reformat if required from _file_format to file_format
        self.status = Status()       #   Contains ERROR/SUCCESS status

    def _validate_input(self):
        try:
            print(f"    {self.__class__.__name__}:_validate_input: +INFO: Setting up the Validator.")
            
            u.file_exists(self.filepath)    #   raise exception

            self._work_file = self.filepath
            self._basedir, self._filename, _ , self._encoded_with = u.split_filepath(self._work_file) #  formats will be upper case

            if not self.outdir:
                self.outdir = u.create_tmpdir(self._basedir) #   raise exception

            self._ext = self.file_format.lower()

            #   check encoding options
            if self._encoded_with:
                u.check_coding_format(self._encoded_with) #   raise exception
            if self.encode_with:
                u.check_coding_format(self.encode_with) #   raise exception
                self._ext += self.encode_with.lower()

            #   check formatting options
            if self._file_format != self.file_format:
                self._reformat = True

            self._outpath = self.outdir + "/" + self._filename 
            self._outpath += f"_{self.outsuffix}" + f".{self._ext}" if self.outsuffix else f".{self._ext}" 

            print(f"    {self.__class__.__name__}:_validate_input: |--Output will be directed to {self._outpath}")
            # print(f"    {self.__class__.__name__}:_validate_input: |--DONE\n")

        except Exception as e:
            print(f"    {self.__class__.__name__}:_validate_input: |--ERROR: {e}")
            raise


    def _cleanup(self):
        pass
        # print(f"    {self.__class__.__name__}:_cleanup: +INFO: Cleaning output directory {self.outdir}")
        # u.remove_dir(self.outdir)
        # print(f"    {self.__class__.__name__}:_cleanup: |--DONE\n")


    def run(self):
        try:
            self._validate_input()  #   get required input/output paths and names

            #   resolve input file encoding
            if self._encoded_with:
                tmpfilepath = '.'.join(self._work_file.split('.')[:-1])
                r.remove_coding(self._encoded_with,filepath=self._work_file,outpath=tmpfilepath)
                self._work_file = tmpfilepath

            #   TODO - editation before or after formatting??
            #   validate file
            self.validate()

            #   format file
            if self._reformat:
                self.format(self._file_format, self.file_format)

            #   resolve output file encoding
            if self.encode_with:
                r.add_coding(self.encode_with,filepath=self._work_file,outpath=self._outpath)

            #   move/save new file
            if self._work_file != self.filepath:
                #   move new file
                u.move_file(self._work_file,self._outpath)
            else:
                #   create copy of the original file?? TODO: maybe it is better to move it - damage original file config information - useless when an error is raised - then everything should be reverted, but how?
                u.copy_file(self._work_file,self._outpath)

        except Exception as e:
            print(f"    {self.__class__.__name__}:run: |--ERROR: {e}")
            self.status.log(e)
            self._cleanup()
            raise
        finally:
            return self.status


    def validate(self):
        raise NotImplementedError(f"Function validate() Not implemented")


    def format(self, format_to):
        raise NotImplementedError("Function format() Not implemented")


    def __str__(self):
        """Human-friendly summary string."""
        return (
            f"Validator(\n"
            f"User entered\n"
            f"  filepath            = {self.filepath}\n"
            f"With characteristics\n"
            f"  encoded_with        = {self._encoded_with}\n"
            f"  basedir             = {self._basedir}\n"
            f"  file_format         = {self._file_format}\n"
            f"  filename            = {self._filename}\n"
            f"  work_file           = {self._work_file}\n"

            f"Output settings\n"
            f"  encode_with         = {self.encode_with}\n"
            f"  output_file_format  = {self.file_format}\n"
            f"  outsuffix           = {self.outsuffix}\n"
            f"  outdir              = {self.outdir}\n"
            f"  output_ext          = {self._ext}\n"
            f"  outpath             = {self._outpath}\n"
            f"Current Status\n"
            f"  status              = {self.status.get_status()[0]}\n"
            f")"
        )


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Subclassed Validators
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

class FASTA_Validator(Validator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._file_format = "FASTA"


    def validate(self):
        print(f"    {self.__class__.__name__}:validate: +INFO: Validating FASTA file")
        try:
            # Try to parse the file as FASTA using Bio.SeqIO
            with open(self._work_file, 'r') as fa:
                records = list(SeqIO.parse(fa, "fasta"))
                if not records:
                    raise ValueError("No sequences found in FASTA file.")
        except Exception as e:
            print(f"    {self.__class__.__name__}:validate: |--ERROR: {e}")
            self.status.log(e)
            raise
        # print(f"    {self.__class__.__name__}:validate: |--DONE\n")


    def edit(self):
        print(f"    {self.__class__.__name__}:edit: +INFO: Editing FASTA file")
        records = None
        try:
            with open(self._work_file,'r') as handle:
                records = list(SeqIO.parse(handle, "fasta"))

            if len(records) > 1:
                print(f"    {self.__class__.__name__}:edit: |--WARNING: more than one sequence found, first will be taken as reference, the rest will be suppressed")
                for i,record in enumerate(records):
                    if i == 0:
                        #   first seq = ref
                        with open(self._work_file,'w') as f:
                            f.write(">chr1 " + record.name+'\n')
                            f.write(str(record.seq)+'\n')
                    else:
                        p = self._work_file[:-6]+f"_{i}.fasta"
                        with open(p,'w') as f:
                            f.write(">"+record.name+'\n')
                            f.write(str(record.seq)+'\n')
        except Exception as e:
            print(f"    {self.__class__.__name__}:edit: |--ERROR: Edit FASTA: {e}")
            self.status.log(e)
        # print(f"    {self.__class__.__name__}:edit: |--DONE\n")


    def format(self,format_to):
        print(f"    {self.__class__.__name__}:format: +INFO: Formatting FASTA file to {format_to}")
        u.TODO(pref=f"    {self.__class__.__name__}:format: |--")
        # print(f"    {self.__class__.__name__}:format: |--DONE\n")



class GBK_Validator(Validator):
    def validate(self):
        print(f"    {self.__class__.__name__}:validate: +INFO: Validating Gene Bank file")
        u.TODO(pref=f"    {self.__class__.__name__}:validate: |--")
        # print(f"    {self.__class__.__name__}:validate: |--DONE\n")


    def format(self,format_to):
        print(f"    {self.__class__.__name__}:format: +INFO: Formatting Gene Bank file to {format_to}")
        u.TODO(pref=f"    {self.__class__.__name__}:format: |--")
        # print(f"    {self.__class__.__name__}:format: |--DONE\n")



class FASTQ_Validator(Validator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._file_format = "FASTQ"


    def validate(self):
        print(f"    {self.__class__.__name__}:validate: +INFO: Validating FASTQ file")
        try:
            # Try to parse the file as FASTQ using Bio.SeqIO
            with open(self._work_file, 'r') as fq:
                records = list(SeqIO.parse(fq, "fastq"))
                if not records:
                    raise ValueError("No sequences found in FASTQ file.")
        except Exception as e:
            print(f"    {self.__class__.__name__}:validate: |--ERROR: Not a valid FASTQ file")
            self.status.log(e)
        # print(f"    {self.__class__.__name__}:validate: |--DONE\n")


    def format(self,format_to):
        print(f"    {self.__class__.__name__}:format: +INFO: Formatting FASTQ file to {format_to}")
        u.TODO(pref=f"    {self.__class__.__name__}:format: |--")
        # print(f"    {self.__class__.__name__}:format: |--DONE\n")



class BAM_Validator(Validator):
    def validate(self):
        print(f"    {self.__class__.__name__}:validate: +INFO: Validating BAM file")
        u.TODO(pref=f"    {self.__class__.__name__}:validate: |--")
        # print(f"    {self.__class__.__name__}:validate: |--DONE\n")

    def format(self,format_to):
        print(f"    {self.__class__.__name__}:format: +INFO: Formatting BAM file to {format_to}")
        u.TODO(pref=f"    {self.__class__.__name__}:format: |--")
        # print(f"    {self.__class__.__name__}:format: |--DONE\n")



class GTF_Validator(Validator):
    def validate(self):
        print(f"    {self.__class__.__name__}:validate: +INFO: Validating GTF file")
        u.TODO(pref=f"    {self.__class__.__name__}:validate: |--")
        # print(f"    {self.__class__.__name__}:validate: |--DONE\n")

    def format(self,format_to):
        print(f"    {self.__class__.__name__}:format: +INFO: Formatting GTF file to {format_to}")
        u.TODO(pref=f"    {self.__class__.__name__}:format: |--")
        # print(f"    {self.__class__.__name__}:format: |--DONE\n")



class GFF_Validator(Validator):
    def validate(self):
        print(f"    {self.__class__.__name__}:validate: +INFO: Validating GFF file")
        u.TODO(pref=f"    {self.__class__.__name__}:validate: |--")
        # print(f"    {self.__class__.__name__}:validate: |--DONE\n")

    def format(self,format_to):
        print(f"    {self.__class__.__name__}:format: +INFO: Formatting GFF file to {format_to}")
        u.TODO(pref=f"    {self.__class__.__name__}:format: |--")
        # print(f"    {self.__class__.__name__}:format: |--DONE\n")
