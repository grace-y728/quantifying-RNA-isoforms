import argparse
import zipfile
import copy
import scipy as sp
from sklearn.linear_model import ElasticNet


def parse_annotation_file(annotation_fn):
  """
  :param annotation_fn: the gene annotations file
  :return: outputs a list of tuples, "genes". Every tuple represents one gene, the first element of the tuple is the
          list of exon index ranges for that gene, the second element of the tuple is the list of isoforms that exist
          for that gene. For example, genes[0][0][0] references the first exon range of the first gene (in a tuple),
          genes[2][1][3] references the fourth isoform of the third gene.
  """
  with open(annotation_fn, 'r') as aFile:
    N = int(aFile.readline().strip())
    genes = [None]*N
    for i in range(N):
      numExons = int(aFile.readline().strip())
      exons = [None]*numExons
      starts = [int(x) for x in aFile.readline().strip().split(' ')]
      ends = [int(x) for x in aFile.readline().strip().split(' ')]
      for j in range(numExons):
        exons[j] = (starts[j], ends[j])
      numIsoforms = int(aFile.readline().strip())
      isoforms = [None]*numIsoforms
      for j in range(numIsoforms):
        isoforms[j] = [int(x) for x in aFile.readline().strip().split(' ')]
        genes[i] = (exons, isoforms)
  return genes


def parse_genome_file(genome_fn):
  """
  :param genome_fn: the full genome file
  :return: the string containing the genome
  """
  with open(genome_fn, 'r') as gFile:
    return gFile.readline().strip()


def parse_reads_file(reads_fn):
  """
  :param reads_fn: the file of shuffled reads
  :return: a list containing all of the shuffled reads
  """
  out_reads = []
  with open(reads_fn, 'r') as rFile:
    for line in rFile:
      out_reads.append(line.strip())
  return out_reads


def generate_junctions(genes):
  # for each isoform of each gene parse each gene's junction list and add each number as a string to an array
  # for each of those arrays, iterate through the string and add each consecutive pair of 2 to another array
  # this time keep only 1 array for each gene
  # remove duplicates
  # sort by alphabetical order using sort()
  # return array containing the junctions for each gene
  gene_junctions = []
  for gene in genes:
    iso_list = []
    for isoform in gene[1]:
      exon_list = []
      for exon in isoform:
        exon_list.append(exon)
      iso_list.append(exon_list)
    arr = []
    for isoform in iso_list:
      for i in range(len(isoform)-1):
        arr.append(str(isoform[i]) + str(isoform[i+1]))
    arr2 = list(dict.fromkeys(arr))
    arr2.sort()
    gene_junctions.append(arr2)
  return gene_junctions


def match_read_exon(exon, read):
  '''helper function that returns whether a read maps to a given exon'''
  return read in exon


def construct_iso_transcripts(genes, genome):
  '''Returns a list of isoform transcripts for each gene'''
  all_gene_transcripts = []
  for gene in genes:
    gene_transcripts = []
    isoforms = gene[1]
    for isoform in isoforms:
      isoform_transcript = ""
      for exon in isoform:
        isoform_transcript += genome[gene[0][exon][0]:gene[0][exon][1]+1]
      gene_transcripts.append(isoform_transcript)
    all_gene_transcripts.append(gene_transcripts)
  return all_gene_transcripts


def construct_isoforms(genes, master_elist):
  all_isoforms = []
  for gene in genes:
    gene_isoforms = []
    for isoform in gene[1]:
      exon_list = []
      count = 0
      for exon in isoform:
        exon_list.append(master_elist[genes.index(gene)][exon])
        count += len(master_elist[genes.index(gene)][exon])
      for i in range(len(isoform)-1):
        left = master_elist[genes.index(gene)][isoform[i]][-49:]
        right = master_elist[genes.index(gene)][isoform[i+1]][0:49]
        exon_list.append(left + right)
      gene_isoforms.append(exon_list)
    all_isoforms.append(gene_isoforms)
  return all_isoforms


# return an array for each gene, kinda in a column of all exons followed by all junctions
def get_exon_counts(all_genes, reads):
  all_counts = []
  for gene in all_genes:
    gene_counts = []
    for _ in gene:
      gene_counts.append(0)
    all_counts.append(gene_counts)
  for read in reads:
    for gene in all_genes:
      for segment in gene:
        if (match_read_exon(segment, read)):
          all_counts[all_genes.index(gene)][gene.index(segment)] += 1
  return all_counts


def curr_junctions_helper(gene_num, genes):
  ''' Small helper function that maps a gene to its junctions '''
  all_junctions = generate_junctions(genes)
  return all_junctions[gene_num-1]


def construct_array(genes, gene_exon_junctions, gene_isoforms, gene_num):
  '''
  Constructs the input A array for Ax=b computation

  The array has a row for each isoform and is populated with frequencies
  to the corresponding exons and junctions the isoform has
  '''
  gene_array = []
  for _ in gene_isoforms:
    row = []
    for _ in gene_exon_junctions:
      row.append(0)
    gene_array.append(row)

  curr_gene = genes[gene_num-1]
  curr_junctions = curr_junctions_helper(gene_num, genes)

  for isoform in curr_gene[1]:
    ind = curr_gene[1].index(isoform)
    junction_list = []
    for i in range(len(isoform)-1):
      junction_list.append(str(isoform[i])+str(isoform[i+1]))

    for exon in isoform:
      gene_array[ind][exon] = len(gene_exon_junctions[exon])/50.0

    curr_index = len(curr_gene[0])
    for junction in curr_junctions:
      if junction in junction_list:
        j_index = curr_junctions.index(junction)
        gene_array[ind][curr_index + j_index] = len(gene_exon_junctions[curr_index + j_index])/50.0
  return gene_array


# master list of all exons
def construct_gene_exons(genes, genome):
  all_exon_transcripts = []
  for gene in genes:
    exon_transcripts = []
    exons = gene[0]
    for exon in exons:
      exon_transcripts.append(genome[exon[0]:(exon[1]+1)])
    all_exon_transcripts.append(exon_transcripts)
  return all_exon_transcripts


def add_junctions(exons, genes):
  '''master list of all possible exons and junctions'''
  gene_junctions = generate_junctions(genes)
  copy_exons = copy.deepcopy(exons)
  index = 0
  for gene_exons in exons:
    # for each gene
    for i in gene_junctions[exons.index(gene_exons)]:
      first_exon = int(i[0])
      second_exon = int(i[1])
      thing1 = gene_exons[first_exon][-49:]
      thing2 = gene_exons[second_exon][0:49]
      concat_junction = thing1 + thing2
      copy_exons[index].append(concat_junction)
    index += 1
  return copy_exons


def quantify_isoforms(genes, genome, reads):
    """
    :param genes: the list of gene tuples generated by the parser
    :param genome_fn: the full genome file
    :param reads_fn: the file of shuffled reads
    :return: a list of tuples, where the first element of the tuple is the transcript sequence (the isoform in terms of
            the exon sequences that form it in the genome), and the second element of the tuple is the abundance of that
            specific isoform
    """
    all_transcripts = construct_iso_transcripts(genes, genome)
    print("got transcripts")
    all_exons = construct_gene_exons(genes, genome)
    print("got gene exons")
    all_exon_junctions = add_junctions(all_exons, genes)
    print("got junctions added")
    all_isoforms = construct_isoforms(genes, all_exons)
    print("got isoforms")
    all_exon_frequencies = get_exon_counts(all_exon_junctions, reads)
    print("got frequencies")
    gene_arrays = []
    for gene_isoform in all_isoforms:
        ind = all_isoforms.index(gene_isoform)
        g_arr = construct_array(genes, all_exon_junctions[ind], gene_isoform, ind+1)
        gene_arrays.append(g_arr)

    all_gene_array_transpose = []
    for ga in gene_arrays:
        ga_transpose = []
        for _ in ga[0]:
            ga_transpose.append([])
        for i in range(len(ga)):
            for j in range(len(ga[i])):
                ga_transpose[j].append(ga[i][j])
        all_gene_array_transpose.append(ga_transpose)

    all_abundances = []
    for gene_array in all_gene_array_transpose:
        g_index = all_gene_array_transpose.index(gene_array)
        a_matrix = sp.array(gene_array)
        b_matrix = sp.array(all_exon_frequencies[g_index])
        coverage = ElasticNet(alpha=1.5, positive=True) # do not allow negative coverage
        coverage.fit(a_matrix, b_matrix)
        coverage_sum = float(sum(coverage.coef_))
        abundance = [x / coverage_sum for x in coverage.coef_]
        all_abundances.append(abundance)

    final_array = []
    for gene_transcripts in all_transcripts:
        ind = all_transcripts.index(gene_transcripts)
        for transcript in gene_transcripts:
            ind2 = gene_transcripts.index(transcript)
            tup = (transcript, all_abundances[ind][ind2])
            final_array.append(tup)
    return final_array


if __name__ == "__main__":
  """
  Usage: python proj4.py -g full_genome.txt -r shuffled_reads.txt -a DATA_PA_1100_0 -o test.out -t hw4_r_4_chr_1
  """
  parser = argparse.ArgumentParser(description='For now this starter code helps parse the files given, but leaves\n'
                                                'the actual function that must be implemented empty')
  parser.add_argument('-g', '--genome', required=True, dest='genome_file', help='File containing the full genome')
  parser.add_argument('-r', '--reads', required=True, dest='read_file', help='File containing the shuffled reads')
  parser.add_argument('-a', '--annotation', required=True, dest='annotation_file', help='File containing gene '
                                                                                        'annotations')
  parser.add_argument('-o', '--outputFile', required=True, dest='output_file', help='Output file name')
  parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                      help='String that needs to be output on the first line of the output file so that the online\n'
                            'submission system recognizes which leaderboard this file should be submitted to. For\n'
                            'hw4, this will be hw4_r_4_chr_1')

  args = parser.parse_args()
  genome_fn = args.genome_file
  reads_fn = args.read_file
  annotation_fn = args.annotation_file
  output_fn = args.output_file

  genes = parse_annotation_file(annotation_fn)
  genome = parse_genome_file(genome_fn)
  reads = parse_reads_file(reads_fn)

  output = quantify_isoforms(genes, genome, reads)
  with open(output_fn, 'w') as oFile:
    oFile.write('>' + args.output_header + '\n')
    oFile.write('>RNA\n')
    for isoform in output:
      out_str = '{} {}\n'.format(isoform[0], isoform[1])
      oFile.write(out_str)

  zip_fn = output_fn + '.zip'
  with zipfile.ZipFile(zip_fn, 'w') as zFile:
    zFile.write(output_fn)
