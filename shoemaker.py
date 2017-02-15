


class PairsHeader:

    # This keys cannot be used as a key in additional header lines
    reserved_header_keys = ['shape','sorted','chromosomes','genome_assembly','columns']

    def __init__(self, chromosomes, shape = 'U', sorted = 'chr1-chr2-pos1-pos2', genome_assembly = None, additional_columns = [], additional_lines = dict()):
        """Shape = 'U' (upper triangle) or 'L' (lower triangle)
           Chr_list is the list of all chromosomes. 
           Chromosomes of each entry should be in this list, but not every chromosome in this list need to have an entry.
           Chr_list will be sorted internally for consistency with column sorting (to ensure upper or lower triangle).
           Column sorting is alpha-numeric (ordering example ['1', '11', '2', '3', '6', 'X', 'X12', 'XXX', 'Y', '_'] )
           Field 'sorted' is for line sorting not column sorting. For now only 'chr1-chr2-pos1-pos2' or 'chr1-pos1-chr2-pos2' are allowed.
           Genome assembly can be an artibrary string. e.g. 'hg19.female', 'GRCh38.mainonly', 'mm9', ... 
           Additional_columns is a list of column names to add other than the 7 reserved columns.
           Additional_lines is a dictionary of additional information. Each item will be displayed as '#key : value'
        """
        assert chromosomes is not None
        assert shape == 'U' or shape == 'L'
        assert sorted == 'chr1-chr2-pos1-pos2' or 'chr1-chr2-pos1-pos2'
        
        self.shape = shape
        self.sorted = sorted
        self.chromosomes = chromosomes
        self.chromosomes.sort() # alpha-numeric sorting. 
        self.genome_assembly = genome_assembly # for humans not computers

        # These first 7 columns are reserved. One can add more columns.
        self.columns = ['readID', 'chr1', 'position1', 'chr2', 'position2', 'strand1', 'strand2'] 
        # Add more user-defined columns
        if len(additional_columns)>0:
            for c in additional_columns:
                assert c not in self.columns
                self.columns.append(c)

        # Additional lines
        self.additional_lines = dict()
        assert isinstance(additional_lines, dict)
        for k, v in additional_lines.items():
            assert k not in self.reserved_header_keys
            self.additional_lines.update({k: v})

    def format_additional_lines(self):
        str = ''
        for k, v in self.additional_lines.items():
            str += "#{} : {}\n".format(k, v)
        return(str)

    def print_header(self, fout):
        if self.shape == 'U':
            shape = 'upper triangle'
        else:
            shape = 'lower triangle'
        fout.write( "#shape : {}\n".format(shape) )
        fout.write( "#sorted : {}\n".format(self.sorted) )
        fout.write( "#chromosomes : {}\n".format(' '.join(self.chromosomes)) )
        if self.genome_assembly is not None:
            fout.write( "#genome_assembly : {}\n".format(self.genome_assembly) )
        fout.write( "#columns : {}\n".format(' '.join(self.columns)) )
        additional_lines_formatted = self.format_additional_lines()
        if(additional_lines_formatted is not ''):
            fout.write(additional_lines_formatted)
        
class PairsLine:
    def __init__(self, header):
        self.header = header

    def add(self, chr1, pos1, chr2, pos2, strand1='.', strand2='.', readID='.', optional=dict()):
        """Add entries, making sure it's either upper or lower triangle according to the predefined shape, based on alpha-numeric sort.
           chr1,pos1,chr2,pos2 cannot be NA ('.').
        """
        assert chr1 in self.header.chromosomes
        assert chr2 in self.header.chromosomes
        assert isinstance(pos1, int) and pos1 > 0
        assert isinstance(pos2, int) and pos2 > 0

        # read ID
        self.readID = readID

        # optional columns
        self.optional = dict()
        if len(optional)>0:
            self.add_optional(optional)

        # main columns, swap two mates to ensure either Upper Triangle or Lower Triangle.
        if self.header.shape == 'U' and chr1 > chr2 or (chr1 == chr2 and pos1 > pos2):
            self.chr1 = chr2
            self.chr2 = chr1
            self.pos1 = pos2
            self.pos2 = pos1
            self.strand1 = strand2
            self.strand2 = strand1
        elif self.header.shape == 'L' and chr2 > chr1 or (chr1 == chr2 and pos1 < pos2):
            self.chr1 = chr2
            self.chr2 = chr1
            self.pos1 = pos2
            self.pos2 = pos1
            self.strand1 = strand2
            self.strand2 = strand1
        else:
            self.chr1 = chr1
            self.chr2 = chr2
            self.pos1 = pos1
            self.pos2 = pos2
            self.strand1 = strand1
            self.strand2 = strand2

    def add_optional(self, keyvaluepairs):
        """Add optional columns (it can be used as an exported function, but it is called by the add method.
           keyvaluepairs is a dictionary with columns names as keys and values as values
        """
        for k,v in keyvaluepairs.items():
            assert k in self.header.columns
            self.optional.update({k: v}) 
        
    def printline(self, fout):
        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.readID, self.chr1, self.pos1, self.chr2, self.pos2, self.strand1, self.strand2))
        if len(self.header.columns) > 7:
            for opt in self.optional.keys():
                fout.write("\t{}".format(self.optional[opt]))
        fout.write("\n")


if __name__ == '__main__':
    p = PairsLine(PairsHeader(['chr1','chr2'], genome_assembly='hg19', additional_columns = ['mismatch'], additional_lines = {'creator': 'haha'}))
    with open('lala','w') as f:
       p.header.print_header(f)
       p.add('chr1',1000,'chr2',3000)
       p.add_optional({'mismatch': 'chr1:1070:A>T'})
       p.printline(f) 
       p.add('chr2',4000,'chr2',5000, '+', '-', optional = {'mismatch': 'chr2:4916:G>C, chr2:5106:C>T'})
       p.printline(f) 
       p.add('chr1',1000,'chr2',3000)
       p.printline(f) 



