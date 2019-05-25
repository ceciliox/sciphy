#! /usr/bin/env python

# TITULO             : Mafft: Alinhamento de sequencias
# AUTOR              : Kary Soriano
# DATA               : 11/01/2008
# DIFICULDADE        : 1
# ==============================================================================
# Objetivo do script: Executado do myscript_nucl.py
#                     Executa Mafft 
# ==============================================================================
# Data da ultima alteracao do script: 30/01/2008
#                                   : 15/01/2008
# ==============================================================================
#-------------------------------------------------------------------------------
# declarando os modulos a usar 
#-------------------------------------------------------------------------------
import os, sys, subprocess, re, pickle, string 
#-------------------------------------------------------------------------------
# modulos do DfAnalyzer
#-------------------------------------
from dfa_lib_python.dataflow import Dataflow
from dfa_lib_python.transformation import Transformation
from dfa_lib_python.attribute import Attribute
from dfa_lib_python.attribute_type import AttributeType
from dfa_lib_python.set import Set
from dfa_lib_python.set_type import SetType
from dfa_lib_python.task import Task
from dfa_lib_python.dataset import DataSet
from dfa_lib_python.element import Element
from dfa_lib_python.dependency import Dependency
#-------------------------------------
#dirin_do_ficheiro = sys.argv[0]
#dirin_arg_pas = sys.argv[0:]
###print "O nome do diretorio de entrada do ficheiro e: " + dirin_do_ficheiro 
###print "E os argumentos passados sao: " + str(dirin_arg_pas)


############################
#PROVENIÊNCIA
############################

dataflow_tag = "mafft-df"
df = Dataflow(dataflow_tag)

##PROVENIÊNCIA PROSPECTIVA
#Transformação para extrair nome dos arquivos: ExtrairNome
tf1 = Transformation("ExtrairNome")
tf1_input = Set("iExtrairNome", SetType.INPUT,
  [Attribute("DIRIN_FILE", AttributeType.FILE)])
tf1_output = Set("oExtrairNome", SetType.OUTPUT,
  [Attribute("FASTA_FILE", AttributeType.FILE)])
tf1.set_sets([tf1_input, tf1_output])
df.add_transformation(tf1)

#Transformação para ler o arquivo e contar o numero de sequencias: ContarSequencias
tf2 = Transformation("ContarSequencias")
tf2_input = Set("iContarSequencias", SetType.INPUT,
  [Attribute("FASTA_FILE", AttributeType.FILE)])
tf2_output = Set("oContarSequencias", SetType.OUTPUT,
  [Attribute("NUMERO_SEQUENCIAS", AttributeType.NUMERIC)])

tf1_output.set_type(SetType.INPUT)
tf1_output.dependency=tf1._tag

tf2.set_sets([tf1_output, tf2_output])
df.add_transformation(tf2)

#Transormação para criar um alinhamento mafft: CriarAlinhamento
tf3 = Transformation("CriarAlinhamento")
tf3_input = Set("iCriarAlinhamento", SetType.INPUT,
  [Attribute("NUMERO_SEQUENCIAS", AttributeType.NUMERIC)])
tf3_output = Set("oCriarAlinhamento", SetType.OUTPUT,
  [Attribute("ALINHAMENTO", AttributeType.TEXT)])

tf2_output.set_type(SetType.INPUT)
tf2_output.dependency=tf2._tag

tf3.set_sets([tf2_output, tf3_output])
df.add_transformation(tf3)
df.save()

t1 = Task(1, dataflow_tag, "ExtrairNome")
t1_input = DataSet("iExtrairNome", [Element(["dirin"])])
t1.add_dataset(t1_input)
t1.begin()

dirin = "/home/linux/"
# Abrindo o diretorio....
sequence_count = '0';  #contar as sequencias pro ajuste do MAFFT
# For file in os.listdir (dirin_do_ficheiro):
for file in os.listdir (dirin):
  if re.search ('.fasta$', file) is not None:
    #--- Mount directory and separate filename
    diretorio = file.split("/")
    nome = diretorio[diretorio.count('/')-1].split(".")  
    #--- Mount name of fasta and mafft files
    fasta = os.path.join (dirin, file)
    mafft = os.path.join (dirin, nome[0]+ ".mafft")
    print ("Aligning: ", fasta, " ...")

    t1_output = DataSet("oExtrairNome", [Element([fasta])])
    t1.add_dataset(t1_output)
    t1.end()      
                
# Lendo o arquivo e contando o numero de sequencias para usar a melhor opcao oferecida pelo MAFFT 
    text = open(fasta).read()

    t2 = Task(2, dataflow_tag, "ContarSequencias")
    t2_input = DataSet("iContarSequencias", [Element([fasta])])#text-fasta
    t2.add_dataset(t2_input)
    t2.begin()

    sequence_count = str.count(text, '>')

    t2_output = DataSet("oContarSequencias", [Element([sequence_count])])
    t2.add_dataset(t2_output)
    t2.end()

    t3 = Task(3, dataflow_tag, "CriarAlinhamento")
    t3_input = DataSet("iCriarAlinhamento", [Element([sequence_count])])
    t3.add_dataset(t3_input)
    t3.begin()

  # Criar um alinhamento com mafft dependendo do numero de seqs, buscando maior precisao e eficiencia
    if sequence_count < 200:
      print ("\tCreating a multiple alignment for {} using MAFFT (L-INS-i < 200 seqs)...".format(fasta))
            # Foram usados para todos os casos, a linha de comando em extenso, pode ser usado os alias fornecidos com o programa
              #1 Builds the command line with a program name and the arguments.
              #2 Runs the command and stores a handle in the handle variable. A handle for a command is the same kind of objects as a file handle: you open it (with the popen command, read from it, and close it.
              #3 Reads all the lines from the handle, and prints the joint result.

      cmd = subprocess.call(["mafft --quiet --localpair --maxiterate 1000 {} > {}".format(fasta, mafft)], shell=True)      

      t3_output = DataSet("oCriarAlinhamento", [Element([mafft])])
      t3.add_dataset(t3_output)  
    #Limpando o contador
      sequence_count = '0';
    elif sequence_count > 2000:
      print ("\tCreating a multiple alignment for {} using MAFFT (FFT-NS-1 > 2000 seqs)...".format(fasta)) 
      cmd = subprocess.call(["mafft --quiet --retree 1 --maxiterate 0 {} > {}".format(fasta, mafft)], shell=True)

      t3_output = DataSet("oCriarAlinhamento", [Element([mafft])])
      t3.add_dataset(t3_output)     
      sequence_count = '0';
    else:
      print ("\tCreating a multiple alignment for {} using MAFFT (FFT-NS-i)...".format(fasta))       
      cmd = subprocess.call(["mafft --quiet --retree 2 --maxiterate 1000 {} > {}".format(fasta, mafft)], shell=True)

      t3_output = DataSet("oCriarAlinhamento", [Element([mafft])])
      t3.add_dataset(t3_output)    
      sequence_count = '0';
    t3.end()

