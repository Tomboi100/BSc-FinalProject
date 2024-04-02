# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:01:23 2023

@author: Tommy
"""
import gzip
from glob import glob
import matplotlib.pyplot as plt
import os

def openfile(filePath):  # openening and reading the fastafile and gfffile
    if filePath.endswith(".gz"):
        with gzip.open(filePath, 'rt') as f:
            return [l.strip() for l in f.readlines()]
    else:
        with open(filePath, 'r') as f:
            return [l.strip() for l in f.readlines()]

def GcContent(seq):  # calculating the gc content percentage
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100, 6)

def parse_fasta(filename):
    FastaFile = openfile(filename)
    FDict = {}
    chromosomeID = ""
    nucleotide_sequence = ""
    for line in FastaFile:
        line = line.strip()
        if line.startswith(">"):
            if chromosomeID != "":
                FDict[chromosomeID] = nucleotide_sequence
                nucleotide_sequence = ""
            chromosomeID = line.strip().split(' ')[2].split(':')[2]
        else:
            nucleotide_sequence += line
    # adds the last sequence to the dictionary because there isn't a > after the last line
    FDict[chromosomeID] = nucleotide_sequence
    return FDict

def parse_gff(gff_file):
    gff_dict = {}
    f = openfile(gff_file)
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] == 'gene':
            gene_id = fields[8].split(';')[0]
            start = int(fields[3])
            stop = int(fields[4])
            chromosome = fields[0]
            gff_dict[gene_id] = (start, stop, chromosome) # creates gff dict
    return gff_dict  # Press Ctrl+F8 to toggle the breakpoint.

def get_End_Seq(gff_dict, FDict, gene):
    stopRegion = gff_dict[gene][1]          #stop coordinates of gene
    chromosome = gff_dict[gene][2]
    endSeq = (FDict[chromosome][stopRegion + 2:stopRegion + 14])
    return endSeq

def get_Initial_StopCodon(gff_dict, FDict, gene):
    stopRegion = gff_dict[gene][1]
    chromosome = gff_dict[gene][2]
    initialStopCodon = (FDict[chromosome][stopRegion - 1:stopRegion+2])# stop codons ????
    return initialStopCodon

def get_start_Seq(gff_dict, FDict, gene):
    startRegion = gff_dict[gene][0]
    chromosome = gff_dict[gene][2]
    startSeq = (FDict[chromosome][startRegion - 13:startRegion - 1]) # front end
    return startSeq

def get_Initial_StartCodon(gff_dict, FDict, gene):
    startRegion = gff_dict[gene][0]
    chromosome = gff_dict[gene][2]
    initialStartCodon = (FDict[chromosome][startRegion - 1:startRegion + 2]) # start codons ???
    return initialStartCodon

def tandemStopCount(endSeq, initialStopCodon, total_stop_codons, total_stop_codons_dict):
    stop_codons = ["TAA", "TAG", "TGA"] # list of possible stop codons
    initialStopCodon_count = 0
    tan_TAA_count = 0
    tan_TAG_count = 0
    tan_TGA_count = 0
    freq = {codon: 0 for codon in stop_codons}
    indexes = {codon: [] for codon in stop_codons}
    for i in range(0, len(endSeq) - 2, 3):  # loops over all stop codons in the endseq
        if len(endSeq) - i >= 3:  # Check if there are at least 3 characters left to form a codon
            codon = endSeq[i:i + 3]
            if codon in stop_codons:
                if (i + 3) % 3 == 0:
                    freq[codon] += 1
                    indexes[codon].append(i)
                    if codon == 'TAA':
                        tan_TAA_count += 1 # iterates to count the condon usage
                    if codon == 'TAG':
                        tan_TAG_count += 1
                    if codon == 'TGA':
                        tan_TGA_count += 1
                total_stop_codons += 1
    if initialStopCodon in stop_codons:
        initialStopCodon_count += 1
        total_stop_codons_dict[initialStopCodon] = (
            total_stop_codons_dict[initialStopCodon][0] + initialStopCodon_count, # returns count to the dict to count frequecy on multiple genes
            total_stop_codons_dict[initialStopCodon][1] + tan_TAA_count,
            total_stop_codons_dict[initialStopCodon][2] + tan_TAG_count,
            total_stop_codons_dict[initialStopCodon][3] + tan_TGA_count
        )
    tandemStopCountData = [freq, indexes, total_stop_codons, total_stop_codons_dict]
    return tandemStopCountData

def tandemStartCount(startSeq, initialStartCodon, total_start_codons, total_start_codons_dict):
    start_codons = ["ATG", "GTG", "TTG"] # list of possible start codons
    initialStartCodon_count = 0
    tan_ATG_count = 0
    tan_GTG_count = 0
    tan_TTG_count = 0
    freq = {codon: 0 for codon in start_codons}
    indexes = {codon: [] for codon in start_codons}
    for i in range(0, len(startSeq) - 2, 3):  # loops over all stop codons in the startseq
        codon = startSeq[i:i + 3]
        if codon in start_codons:
            if (i + 3) % 3 == 0:
                freq[codon] += 1
                indexes[codon].append(i)
                if codon == 'ATG':
                    tan_ATG_count += 1 # iterates to count the condon usage
                if codon == 'GTG':
                    tan_GTG_count += 1
                if codon == 'TTG':
                    tan_TTG_count += 1
            total_start_codons += 1
    if initialStartCodon in start_codons:
        initialStartCodon_count += 1
        total_start_codons_dict[initialStartCodon] = (
            total_start_codons_dict[initialStartCodon][0] + initialStartCodon_count, # returns count to the dict to count frequecy on multiple genes
            total_start_codons_dict[initialStartCodon][1] + tan_ATG_count,
            total_start_codons_dict[initialStartCodon][2] + tan_GTG_count,
            total_start_codons_dict[initialStartCodon][3] + tan_TTG_count
        )
    tandemStartCountData = [freq, indexes, total_start_codons, total_start_codons_dict]
    return tandemStartCountData

def displayTandemStopCodon(stopFreq, stopIndexes, gene, initialStopCodon, f):
    # print("Stop codon frequencies:", gene + ":")# for terminal output
    # print("The initial stop codon was: " + initialStopCodon)
    f.write("Stop codon frequencies: " + gene + ":\n")
    f.write("The initial stop codon was: " + initialStopCodon + "\n")
    for codon, count in stopFreq.items():
        # print(codon, ":", count, "at indexes", stopIndexes[codon])
        f.write(codon + " : " + str(count) + " at indexes " + str(stopIndexes[codon]) + "\n\n")
    return f

def displayTandemStartCodon(startFreq, startIndexes, gene, initialStartCodon, f):
    # print("Start codon frequencies:", gene + ":")# for terminal output
    # print("The initial start codon was: " + initialStartCodon)
    f.write("Start codon frequencies: " + gene + ":\n")
    f.write("The initial start codon was: " + initialStartCodon + "\n")
    for codon, count in startFreq.items():
        #print(codon, ":", count, "at indexes", startIndexes[codon])
        f.write(codon + " : " + str(count) + " at indexes " + str(startIndexes[codon]) + "\n\n")
    return f

def displayFinalinfo(total_stop_codons, total_start_codons,GcResults, element, f):
    # print("Total tandem stop codons found: ", total_stop_codons, "\n")# for terminal output
    # print("Total tandem start codons found: ", total_start_codons, "\n")
    # print("Total gc content for " + element + " = " + (str(sum(GcResults.values()) / len(GcResults))) + "\n\n")
    f.write("Total tandem stop codons found: " + str(total_stop_codons))
    f.write("Total tandem start codons found: " + str(total_start_codons))
    f.write("Total gc content for " + element + " = " + (str(sum(GcResults.values()) / len(GcResults))) + "\n\n")
    return f

def barchartCreate(element, total_stop_codons_dict, total_start_codons_dict):
    initial_stop_codons = ['TAA', 'TAG', 'TGA']
    #actual data ??
    initial_counts = [total_stop_codons_dict['TAA'][0], total_stop_codons_dict['TAG'][0], total_stop_codons_dict['TGA'][0]]
    tandem_TAA_counts = [total_stop_codons_dict['TAA'][1], total_stop_codons_dict['TAG'][1], total_stop_codons_dict['TGA'][1]]
    tandem_TAG_counts = [total_stop_codons_dict['TAA'][2], total_stop_codons_dict['TAG'][2], total_stop_codons_dict['TGA'][2]]
    tandem_TGA_counts = [total_stop_codons_dict['TAA'][3], total_stop_codons_dict['TAG'][3], total_stop_codons_dict['TGA'][3]]

    initial_start_codons = ['ATG', 'GTG', 'TTG']
    initial_start_counts = [total_start_codons_dict['ATG'][0], total_start_codons_dict['GTG'][0], total_start_codons_dict['TTG'][0]]
    tandem_ATG_counts = [total_start_codons_dict['ATG'][1], total_start_codons_dict['GTG'][1], total_start_codons_dict['TTG'][1]]
    tandem_GTG_counts = [total_start_codons_dict['ATG'][2], total_start_codons_dict['GTG'][2], total_start_codons_dict['TTG'][2]]
    tandem_TTG_counts = [total_start_codons_dict['ATG'][3], total_start_codons_dict['GTG'][3], total_start_codons_dict['TTG'][3]]

    # Bar chart setup
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))  # Create subplots for stop codon and start codon bar charts
    width = 0.15  # the width of the bars
    x = [i for i in range(len(initial_stop_codons))]  # the x locations for the groups of bars

    # Stop codon bar chart
    # Plot the initial stop codon counts TAA
    bar1 = ax1.bar([i - 1.5 * width for i in x], initial_counts, width, label='Initial Stop')
    # Plot the tandem stop codon counts TAG
    bar2 = ax1.bar([i - 0.5 * width for i in x], tandem_TAA_counts, width, label='TAA')
    # Plot the tandem stop codon counts TGA
    bar3 = ax1.bar([i + 0.5 * width for i in x], tandem_TAG_counts, width, label='TAG')
    # Plot the tandem stop codon counts TGA
    bar4 = ax1.bar([i + 1.5 * width for i in x], tandem_TGA_counts, width, label='TGA')
    ax1.set_xticks(x)
    ax1.set_xticklabels(initial_stop_codons)
    ax1.legend()
    ax1.set_ylabel('Frequency')
    ax1.set_title('Initial Stop Codon vs Tandem Stop Codons')
    # Start codon bar chart
    # Plot the initial start codon counts ATG
    bar1 = ax2.bar([i - 1.5 * width for i in x], initial_start_counts, width, label='Initial Start')
    # Plot the tandem start codon counts ATG
    bar2 = ax2.bar([i - 0.5 * width for i in x], tandem_ATG_counts, width, label='ATG')
    # Plot the tandem start codon counts GTG
    bar3 = ax2.bar([i + 0.5 * width for i in x], tandem_GTG_counts, width, label='GTG')
    # Plot the tandem start codon counts TTG
    bar4 = ax2.bar([i + 1.5 * width for i in x], tandem_TTG_counts, width, label='TTG')
    ax2.set_xticks(x)
    ax2.set_xticklabels(initial_start_codons)
    ax2.legend()
    ax2.set_ylabel('Frequency')
    ax2.set_title('Initial Start Codon vs Tandem Start Codons')
    element = element.strip().split('.')# name formating
    element = element[0] + '.' + element[1]
    ax1.text(0.5, 1.15, element, transform=ax1.transAxes,
             ha='center', va='top', fontsize=12)
    # save the barchart
    plt.tight_layout()
    outfileName = ('bar_chart_' + element + '.png')
    #directory_path = "C:/Users/Tommy/PycharmProjects/CompleteFinalProject/Output/graphs/" # eukaryotic output
    directory_path = "C:/Users/Tommy/PycharmProjects/CompleteFinalProject/Output/pro_graphs/" # prokaryotic output
    out_filepath = os.path.join(directory_path, outfileName)
    plt.savefig(out_filepath)
    plt.cla()
    plt.close(fig)
    print(element)
    # Show the chart
    #plt.show()

def stopVariableInitialisation():
    # total stop codon counting dict
    total_stop_codons_dict = {}
    initialStopCodon_count = 0
    tan_TAA_count = 0
    tan_TAG_count = 0
    tan_TGA_count = 0
    total_stop_codons_dict["TAA"] = (initialStopCodon_count, tan_TAA_count, tan_TAG_count, tan_TGA_count)
    total_stop_codons_dict["TAG"] = (initialStopCodon_count, tan_TAA_count, tan_TAG_count, tan_TGA_count)
    total_stop_codons_dict["TGA"] = (initialStopCodon_count, tan_TAA_count, tan_TAG_count, tan_TGA_count)
    return total_stop_codons_dict

def startVariableInitialisation():
    # total start codon counting dict
    total_start_codons_dict = {}
    initialStartCodon_count = 0
    tan_ATG_count = 0
    tan_GTG_count = 0
    tan_TTG_count = 0
    total_start_codons_dict["ATG"] = (initialStartCodon_count, tan_ATG_count, tan_GTG_count, tan_TTG_count)
    total_start_codons_dict["GTG"] = (initialStartCodon_count, tan_ATG_count, tan_GTG_count, tan_TTG_count)
    total_start_codons_dict["TTG"] = (initialStartCodon_count, tan_ATG_count, tan_GTG_count, tan_TTG_count)
    return total_start_codons_dict

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # main start
    # FFiles = glob('C:/Users/Tommy/PycharmProjects/ensembl_fungi_genomes/ensembl_fungi_genomes/' + '*.fa.gz')
    FFiles = glob('C:/Users/Tommy/PycharmProjects/Ensembl_bacteria_genomes/' + '*.fa.gz')

    for element in FFiles:
        fastaName = os.path.basename(element)
        #FDict = parse_fasta(element) moved to after gff is checked to increased speed

        gffFilename = os.path.basename(element)
        gffFilename = gffFilename.strip().split('.')
        # gffFilename = gffFilename[0] + "." + gffFilename[1] + ".52.gff3.gz" # for the eukaryotic genome version
        gffFilename = gffFilename[0] + "." + gffFilename[1] + ".46.gff3.gz" # for the prokaryotic genome version

        # gffPath = os.path.join('C:/Users/Tommy/PycharmProjects/ensembl_fungi_gff/ensembl_fungi_gff/', gffFilename)
        gffPath = os.path.join('C:/Users/Tommy/PycharmProjects/Ensembl_bacteria_gff/', gffFilename)

        if os.path.exists(gffPath):
            FDict = parse_fasta(element) # put fasta file into a dictionary
            gff_dict = parse_gff(gffPath)# puts gff into a dictionary
            #print(element)
            #print(gffPath) for testing
        else:
            print("match not found")
            print(gffPath)
            continue

        # GC content calculation
        GcResults = {key: GcContent(value) for (key, value) in FDict.items()}

        # variable initialisation
        # Total start/stop codons
        total_stop_codons = 0
        total_start_codons = 0
        # total stop codon counting dict
        total_stop_codons_dict = stopVariableInitialisation()
        total_start_codons_dict = startVariableInitialisation()
        # file writing
        filename = fastaName
        if filename.endswith(".gz"):
            filename = filename[:-3]
        outfileName = ("outfile_" + filename + ".txt")
        #directory_path = "C:/Users/Tommy/PycharmProjects/CompleteFinalProject/Output/mainTextOutput/" # eukaryotic output
        directory_path = "C:/Users/Tommy/PycharmProjects/CompleteFinalProject/Output/pro_mainTextOutput/" # prokaryotic output
        out_filepath = os.path.join(directory_path, outfileName)

        with open(out_filepath, "w") as f:
            for gene in gff_dict:
                #function calling below to obtain the sequences first then the data on the
                endSeq = get_End_Seq(gff_dict, FDict, gene)
                #print(endSeq) # for testing
                initialStopCodon = get_Initial_StopCodon(gff_dict, FDict, gene)
                #print(initialStopCodon)# for testing
                startSeq = get_start_Seq(gff_dict, FDict, gene)
                # print(startSeq)# for testing
                initialStartCodon = get_Initial_StartCodon(gff_dict, FDict, gene)
                #print(initialStartCodon)
                tandemStopCountData = tandemStopCount(endSeq, initialStopCodon, total_stop_codons, total_stop_codons_dict)
                # tandemstopcount == [0] = freq, [1] = indexes, [2] = total_stop_codons, [3] = total_stop_codons_dict
                tandemStartCountData = tandemStartCount(startSeq, initialStartCodon, total_start_codons, total_start_codons_dict)
                # tandemstartcount == [0] = freq, [1] = indexes, [2] = total_start_codons, [3] = total_start_codons_dict

                f = displayTandemStopCodon(tandemStopCountData[0], tandemStopCountData[1], gene, initialStopCodon, f)
                f = displayTandemStartCodon(tandemStartCountData[0], tandemStartCountData[1], gene, initialStartCodon, f)

                total_stop_codons = tandemStopCountData[2]
                total_start_codons = tandemStartCountData[2]
                total_stop_codons_dict = tandemStopCountData[3]
                # total_stop_codons_dict == [0] = initialStopCodon_count, [1] = tan_TAA_count, [2] = tan_TAG_count, [3] = tan_TGA_count
                total_start_codons_dict = tandemStartCountData[3]
                # total_start_codons_dict == [0] = initialStartCodon_count, [1] = tan_ATG_count, [2] = tan_GTG_count, [3] = tan_TTG_count

            f = displayFinalinfo(total_stop_codons, total_start_codons,GcResults, fastaName, f)
            f.close()
            barchartCreate(filename, total_stop_codons_dict, total_start_codons_dict)

