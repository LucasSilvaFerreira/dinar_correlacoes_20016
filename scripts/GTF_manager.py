__author__ = 'lucas'
# -*- coding: utf-8 -*-
import sys
import re
from tqdm import tqdm
from collections import namedtuple



class Transcript():
    def __init__(self, transcript_name, exons_array, attr_list):
        self.transcript_name = transcript_name
        self.exons = exons_array
        self.__attr_list = attr_list
        self.exon_count = len(exons_array)
        self.transcript_size = self.__get_transcrip_size()
        #print self.transcript_name, self.transcript_size
        self.__named_tuple_attrs = namedtuple('attr', ' '.join(self.__attr_list) )# Creating named tuple
        self.attrs = self.__parse_attrs() # aqui eu deixarei todos os atributos dos genes


    def exon_size(self, one_exon_in_array):
        return int(one_exon_in_array[4])-int(one_exon_in_array[3])


    def __get_transcrip_size(self):
        return sum([self.exon_size(one_exon_in_array=exon) for exon in self.exons])

    def __parse_attrs(self):
        print self.exons[-1][-1]
        fields_to_parse = {field.strip(' ').split(' ')[0]:field.strip(' ').split(' ')[1]
                           for field in self.exons[-1][-1].split(';')} # get the last exon to parse attrs
        hash_to_load_namedtuple = {}
        for attr in self.__attr_list:
            if attr in fields_to_parse:
                hash_to_load_namedtuple[attr] = fields_to_parse[attr]
            else:
                hash_to_load_namedtuple[attr] = None

        return self.__named_tuple_attrs(**hash_to_load_namedtuple)



class Gene_content():
    def __init__(self, gene_id,
                 strand,
                 transcripts_id_hash,
                 info_gene_str,
                 attr_transcript_field
                 ):

        self.gene_id = gene_id
        self.gene_info = info_gene_str
        self.strand = strand  # checar se strand existe na posicao correta else: return error
        self.transcripts_ids = transcripts_id_hash #sorted_dict
        self.__child_transcript_possible_fields = attr_transcript_field #hash com os valores que os transcritos vao possuir
        self.transcripts_list = self.__parse_transcripts()



    def __parse_transcripts(self):
        hash_transcripts_to_parse = {}
        return {Transcript(transcript_name=ts_key,
                           exons_array=ts_values_array,
                           attr_list=self.__child_transcript_possible_fields) for ts_key, ts_values_array in self.transcripts_ids.iteritems()}





class GTF_manager():
    def __init__(self, gtf_file):
        self.file_name = gtf_file
        self.gtf_file_open = [line_gtf.split('\t') for line_gtf in open(gtf_file, 'r').read().split('\n') if len(line_gtf)> 0 ]
        self.__attr_list = {}
        self.genes_hash = self.__parse_lines()

    def __parse_lines(self):
        hash_gtf_file = {}
        print 'Parsing genes...'

        for line in tqdm(self.gtf_file_open):
            if line:
                # Populating attrs list
                for attr_line in line[-1].split(';'):
                    self.__attr_list[attr_line.strip().split(" ")[0]] = ''
                # Parsing lines
                try:
                    gtf_gene_id = re.search('gene_id "(\S*)"', line[-1]).group(1)
                    if gtf_gene_id in hash_gtf_file:
                        hash_gtf_file[gtf_gene_id].append(line)
                    else:
                        hash_gtf_file[gtf_gene_id] = []
                        hash_gtf_file[gtf_gene_id].append(line)
                except:
                    raise IOError("\n The Line: \n{line}\n don't have a regular gene_id field!".format(line=line))

        gene_list_to_return = {}
        print 'Converting gene to proper format'
        for key_gene, gene_value in hash_gtf_file.iteritems():

            gene_id = key_gene
            gene_info = gene_value[0][1]
            strand = gene_value[0][6]
            transcript_hash = {}
            for transcript in gene_value:  # Loading  each exon in each transcript
                id_line = re.search('transcript_id "(\S*)"', transcript[-1]).group(1)

                if id_line in transcript_hash:
                    transcript_hash[id_line].append(transcript)
                else:
                    transcript_hash[id_line] = []
                    transcript_hash[id_line].append(transcript)
            #Gene content class creation
            gene_object = Gene_content(gene_id= gene_id,
                         strand=strand,
                         transcripts_id_hash=transcript_hash,
                         info_gene_str=gene_info,
                         attr_transcript_field=self.__attr_list.iterkeys())

            if gene_id in gene_list_to_return:
                raise IOError('\n Duplicated gene_id\n {gene_id}\n'.format(gene_id=gene_id))
            else:
                gene_list_to_return[gene_id]= gene_object


        return gene_list_to_return



        print '\n Parsing finished...'



    def print_attrs_fields(self):
        for attr_field_print in self.__attr_list.iterkeys():
            print attr_field_print

# checar tamanho total do transcrito


def main():
    '''Main Function Description'''

    gtf_object = GTF_manager(gtf_file='/home/lucas/PycharmProjects/dinar_correlacoes_20016/dados/results/genes-2.gtf')
    gtf_object.attrs_fields()

if __name__ == '__main__':
    sys.exit(main())
