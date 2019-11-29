# Benchmark Tool created by Nicholas Ho 
# November 2019


import sys
from random import seed
from random import randint
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import os



arguments = len(sys.argv) - 1

arguments_received = False
argument_values = {}

#taking in the arguments.
position = 1
while (arguments >= position):
    # print ("parameter %i: %s" % (position, sys.argv[position]))

    if(sys.argv[position] == "--input"):
        argument_values["input"] = sys.argv[position + 1]

    if(sys.argv[position] == "--segSize"):
        argument_values["segSize"] = sys.argv[position + 1]

    if(sys.argv[position] == "--numPiece"):
        argument_values["numPiece"] = sys.argv[position + 1]

    if(sys.argv[position] == "--outputdir"):
        argument_values["outputdir"] = sys.argv[position + 1]

    position = position + 1



if "input" in argument_values.keys():
    print("Input received...", argument_values['input'])
    arguments_received = True
else:
    print ("Please Specify Input: --input file.txt")
    

if "segSize" in argument_values.keys():
    print("segSize received...", argument_values['segSize'])
else:
    print ("segSize not specified. using default 50")
    argument_values["segSize"] = 50

if "numPiece" in argument_values.keys():
    print("numPiece received...", argument_values['numPiece'])
else:
    print ("numPiece not specified. using default 50")
    argument_values["numPiece"] = 50

if "outputdir" in argument_values.keys():
    print("outputdir received...", argument_values['outputdir'])
else:
    print ("No outputdir specified...putting file in current directory.")
    argument_values["outputdir"] = "benchmark.fasta"



list_of_genome_paths = []
with open(argument_values["input"], 'r') as f:
    for line in f:
        list_of_genome_paths.append(line.rstrip())

# print("Genome Path Lists ",list_of_genome_paths)


def dna_parse(list, segment_size, num_pieces, strain_id, cds_name):
    strain_payload = []
    for i in range(num_pieces):
        posval = randint(0, len(list) - segment_size)
        newDescription = str(segment_size) + "-" + str(num_pieces)
        dna_segment = list[posval: (posval+segment_size)]
        newId = str(strain_id) + "?" + str(cds_name)

        record = SeqRecord((dna_segment),
                   id = newId,
                   description = newDescription)
        strain_payload.append(record)
    return strain_payload

try:
    if(arguments_received):
        outputdir = argument_values["outputdir"]

        #removing the old file if it exists
        with open(outputdir, "w") as output_handle:
            print("cleared")

        count = 0
        for genome_path in list_of_genome_paths:

            records = list(SeqIO.parse(genome_path, "fasta"))
            
            seg_size = argument_values["segSize"]
            num_piece = argument_values["numPiece"]
            thing = os.path.splitext(genome_path)[0]
            genome_name = thing.split("/")[-1]

            print(genome_name)
            for record in records:
                newList = dna_parse(list = record.seq, segment_size = seg_size, num_pieces = num_piece, strain_id = genome_name, cds_name = record.id)
                
                if(count % 50 == 0):
                    print(count, " Genes have been processed...", count/(len(records) * len(list_of_genome_paths)) * 100)
                count += 1
                
                with open(outputdir, "a") as output_handle:
                    SeqIO.write(newList, output_handle, "fasta")

except Exception:
    print("Something went wrong. Make sure the input file exists")